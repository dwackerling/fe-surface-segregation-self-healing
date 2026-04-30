# -*- coding: utf-8 -*-
"""
Analytical segregation model (multilayer) + Monte Carlo minimisation.

Modifications included (per your requests):
1) Plane-dependent elastic screening factors F_i:
   - L atomic layers (L=4 by default), layer i=1 at the surface -> d_i = (i-1)*d_hkl
   - F_i = Aexp(-d_i / lambda), with lambda = a (Fe-BCC lattice parameter, fixed)
   - d_hkl = a / sqrt(h^2 + k^2 + l^2)

2) Save energy term decomposition at the optimum (and optionally at MC trace points):
   - term1..term5, Eel, total (J/mol) and divided-by-Omega (J/m^2)

3) Save a separate MC-trace dataset (one dataset for run evaluation):
   - energy terms + yB vector sampled every MC_LOG_EVERY steps
   - includes delta, acceptance stats, adsorption

4) Elastic E(2) units:
   - Using SI inputs (Pa, m), the Eshelby prefactor gives energy per inclusion (J).
   - Converted to J/mol by multiplying by N_A.

Notes:
- Sign convention for elastic term is kept as in your current code: total = ... - Eel.
  (If the paper uses opposite sign, this is the only place to change.)
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from pycalphad import Database, equilibrium, variables as v

# ============================================================
# CONFIGURATION
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

TDB_PATH = PROJECT_ROOT / "data" / "thermodynamic" / "steel_database_fix.tdb"
SURFACES_JSON_PATH = PROJECT_ROOT / "data" / "surface_energy" / "surfaces.json"

PHASE = ["BCC_A2"]
P = 101325.0  # Pa

# Sweep
T_C_LIST = list(np.arange(400, 701, 25, dtype=float))
PLANES = ["100", "110", "111"]

# MC
L = 4
MC_STEPS = 200_000
DELTA_MAX = 1.5
SEED = 0

# Adaptation
ADAPT_STEPS = 50_000
ADAPT_EVERY = 2_000
TARGET_ACC = 0.55
ADAPT_RATE = 0.15
ALPHA = 1.0  # kT_eff = alpha * R * T

# Clamp for delta during adaptation (your latest setting)
DELTA_CLIP_MIN = 1e-5
DELTA_CLIP_MAX = 0.2

# Elastic E(2)
INCLUDE_ELASTIC = True

# Fe-BCC lattice parameter (fixed, not T-dependent)
A_FE_BCC_M = 2.866e-10  # m  (you can edit)

# MC trace logging (separate dataset)
SAVE_MC_TRACE = True
MC_LOG_EVERY = 500       # sample every N MC steps (per run)
MC_TRACE_MAX_ROWS = None # set to int to cap memory (None = no cap)

# Output
OUT_WIDE = "results_wide.csv"
OUT_LONG = "results_long.csv"
OUT_TRACE = "mc_trace.csv"

# ============================================================
# CONSTANTS
# ============================================================
R = 8.314  # J/mol/K
Z_BCC = 8  # nearest neighbours in BCC (1st shell)
N_A = 6.02214076e23  # 1/mol

# Geometry BCC per plane (Omega in m^2/mol, zl/zv coordination terms)
GEOM_BCC = {
    "100": {"Omega": 4.947e4, "zl": 0, "zv": 4},
    "110": {"Omega": 3.498e4, "zl": 4, "zv": 2},
    "111": {"Omega": 8.568e4, "zl": 0, "zv": 4},
}

# gamma0 of Fe per plane (J/m^2)
GAMMA0_FE_BY_PLANE = {"110": 2.45, "100": 2.50, "111": 2.73}

# Atomic masses (g/mol)
ATOMIC_MASS = {"FE": 55.845, "MO": 95.95, "W": 183.84, "CU": 63.546, "AU": 196.966569}

# Elastic parameters (SI): radii in m, moduli in Pa
FE_PARAMS = {"rA": 1.24e-10, "GA": 82e9}  # Fe metallic (BCC)
SOLUTE_PARAMS = {
    "MO": {"rB": 1.39e-10, "KB": 230e9},
    "W":  {"rB": 1.41e-10, "KB": 310e9},
    "CU": {"rB": 1.28e-10, "KB": 140e9},
    "AU": {"rB": 1.44e-10, "KB": 180e9},
}

# ============================================================
# UTILITIES
# ============================================================
def wt_to_atfrac_binary_feX(wt_solute: float, solute_symbol: str) -> float:
    """Convert Fe-(wt_solute)%X to atomic fraction x_X (binary Fe-X)."""
    sol = solute_symbol.upper()
    if sol not in ATOMIC_MASS:
        raise KeyError(f"Missing atomic mass for {sol}")
    wt_fe = 100.0 - float(wt_solute)
    n_fe = wt_fe / ATOMIC_MASS["FE"]
    n_x = float(wt_solute) / ATOMIC_MASS[sol]
    x_x = n_x / (n_fe + n_x)
    return float(x_x)

def load_gamma0_weighted(path_json: str) -> dict:
    """Load isotropic gamma0 per element from surfaces.json: weighted_surface_energy (J/m^2)."""
    with open(path_json, "r", encoding="utf-8") as f:
        data = json.load(f)
    df = pd.json_normalize(data)
    if "pretty_formula" not in df.columns or "weighted_surface_energy" not in df.columns:
        raise ValueError("surfaces.json must contain 'pretty_formula' and 'weighted_surface_energy'.")

    gamma0 = {}
    for el, sub in df.groupby("pretty_formula"):
        sub2 = sub.sort_values("e_above_hull", ascending=True) if "e_above_hull" in sub.columns else sub
        row = sub2.iloc[0]
        gamma0[el.upper()] = float(row["weighted_surface_energy"])
    return gamma0

def plane_spacing_bcc_m(plane: str, a_m: float) -> float:
    """Interplanar spacing d_hkl = a / sqrt(h^2+k^2+l^2)."""
    h, k, l = (int(plane[0]), int(plane[1]), int(plane[2]))
    return float(a_m / np.sqrt(h*h + k*k + l*l))

def build_F_vector(plane: str, L_layers: int, a_m: float, A_scale: float) -> np.ndarray:
    d_hkl = plane_spacing_bcc_m(plane, a_m)
    lam = float(a_m)
    idx = np.arange(L_layers, dtype=float)
    d_i = idx * d_hkl

    F0 = np.exp(-d_i / lam)               # perfil base
    F = np.clip(A_scale * F0, 0.0, 1.0)   # escala física

    return F.astype(float)

# ============================================================
# Fe-Au: omega(T) from literature L0 (J/mol)
# ============================================================
def omega_FeAu_L0_Jmol(T_K: float) -> float:
    """
    Fe-Au (BCC): w = L0 = 34122.3 − 12.36*T  (T in K), J/mol.
    """
    T_K = float(T_K)
    return float(34122.3 - 12.36 * T_K)

# ============================================================
# CALPHAD: ΔHmix(x) and omega LSQ (J/mol)
# ============================================================
db = Database(TDB_PATH)

def delta_H_mix_Jmol(elem1: str, x1: float, elem2: str, T_K: float) -> float:
    """
    ΔH_mix = H(x) - x H(pure1) - (1-x) H(pure2) in PHASE.
    Returns J/mol.
    """
    elem1 = elem1.upper()
    elem2 = elem2.upper()
    x1 = float(x1)
    x2 = 1.0 - x1
    comps = [elem1, elem2, "VA"]

    H_sys = equilibrium(db, comps, PHASE, {v.T: T_K, v.P: P, v.X(elem1): x1}, output="HM").HM.values.flatten()[0]
    H_1   = equilibrium(db, comps, PHASE, {v.T: T_K, v.P: P, v.X(elem1): 0.99999}, output="HM").HM.values.flatten()[0]
    H_2   = equilibrium(db, comps, PHASE, {v.T: T_K, v.P: P, v.X(elem1): 0.00001}, output="HM").HM.values.flatten()[0]
    return float(H_sys - (x1 * H_1 + x2 * H_2))

def omega_avg_LSQ_Jmol(elem1: str, elem2: str, T_K: float, z: int = Z_BCC,
                       x_min: float = 0.05, x_max: float = 0.95, npts: int = 19) -> float:
    """
    Fit constant omega (J/mol) minimising Σ(ΔH - z x(1-x) omega)^2 over x-range.
    """
    elem1 = elem1.upper()
    elem2 = elem2.upper()
    xs = np.linspace(float(x_min), float(x_max), int(npts))
    xs = np.clip(xs, 1e-3, 1.0 - 1e-3)
    dH = np.array([delta_H_mix_Jmol(elem1, x, elem2, T_K) for x in xs], dtype=float)
    f = float(z) * xs * (1.0 - xs)
    omega = float(np.dot(dH, f) / np.dot(f, f))
    return omega

# ============================================================
# Elastic term E(2): J/mol (SI -> J per inclusion -> *N_A)
# ============================================================
def E2_elastic_Jmol(yB: np.ndarray, xB: float, F: np.ndarray,
                    rA: float, rB: float, GA: float, KB: float) -> float:
    """
    Eshelby-like prefactor in SI gives energy per inclusion (J). Convert to J/mol via N_A.
    Your form: pref * Σ (yB - xB) * F

    Units check (SI):
      pref ~ Pa*m^3 = J  (per solute atom/inclusion)
      summ is dimensionless
      => J, then * N_A => J/mol
    """
    yB = np.asarray(yB, dtype=float)
    F = np.asarray(F, dtype=float)
    xB = float(xB)

    num = 24.0 * np.pi * KB * GA * rA * rB * (rA - rB) ** 2
    den = 3.0 * KB * rB + 4.0 * GA * rA
    pref_J_per_inclusion = num / den
    summ = float(np.sum((yB - xB) * F))
    return float(pref_J_per_inclusion * summ * N_A)

# ============================================================
# gamma*Omega (J/mol) and term breakdown
# ============================================================
def gammaOmega_terms_Jmol(yB: np.ndarray, xB: float,
                          gamma0_A: float, gamma0_B: float,
                          omega_AB: float, Omega: float, T_K: float,
                          zl: int, zv: int,
                          F: np.ndarray,
                          rA: float, rB: float, GA: float, KB: float,
                          include_elastic: bool = False) -> dict:
    """
    Return dict of individual terms (J/mol) and total (J/mol).
    """
    yB = np.asarray(yB, dtype=float)
    xB = float(xB)
    yA = 1.0 - yB
    xA = 1.0 - xB

    eps = 1e-15
    yA = np.clip(yA, eps, 1.0 - eps)
    yB = np.clip(yB, eps, 1.0 - eps)

    term1 = (yA[0] * gamma0_A + yB[0] * gamma0_B) * float(Omega)  # J/mol
    term2 = R * float(T_K) * float(np.sum(yA * np.log(yA / xA) + yB * np.log(yB / xB)))  # J/mol

    term3 = int(zv) * float(omega_AB) * (xA * xB - xA * yB[0] - xB * yA[0])  # J/mol

    yA_next = np.r_[yA[1:], xA]
    yB_next = np.r_[yB[1:], xB]
    term4 = int(zv) * float(omega_AB) * float(np.sum(
        (yB - xB) * (yA_next - xA) +
        (yA - xA) * (yB_next - xB)
    ))

    term5 = int(zl) * float(omega_AB) * float(np.sum((yA - xA) * (yB - xB)))

    Eel = 0.0
    if include_elastic:
        Eel = E2_elastic_Jmol(yB=yB, xB=xB, F=F, rA=rA, rB=rB, GA=GA, KB=KB)

    total = float(term1 + term2 + term3 + term4 + term5 - Eel)  # keep your sign convention

    return {
        "term1_Jmol": float(term1),
        "term2_Jmol": float(term2),
        "term3_Jmol": float(term3),
        "term4_Jmol": float(term4),
        "term5_Jmol": float(term5),
        "Eel_Jmol": float(Eel),
        "total_Jmol": float(total),
    }

def gammaOmega_Jmol(*args, **kwargs) -> float:
    """Compatibility wrapper returning only total (J/mol)."""
    return float(gammaOmega_terms_Jmol(*args, **kwargs)["total_Jmol"])

def adsorption(yB: np.ndarray, xB: float) -> float:
    return float(np.sum(np.asarray(yB, dtype=float) - float(xB)))

# ============================================================
# MC minimisation + optional trace logging
# ============================================================
def monte_carlo_minimise(yB_init: np.ndarray, xB: float,
                         gamma0_A: float, gamma0_B: float,
                         omega_AB: float, Omega: float, T_K: float,
                         zl: int, zv: int,
                         F: np.ndarray,
                         rA: float, rB: float, GA: float, KB: float,
                         steps: int = MC_STEPS,
                         delta_max: float = DELTA_MAX,
                         rng=None,
                         include_elastic: bool = False,
                         alpha: float = ALPHA,
                         adapt_steps: int = ADAPT_STEPS,
                         adapt_every: int = ADAPT_EVERY,
                         target_acc: float = TARGET_ACC,
                         adapt_rate: float = ADAPT_RATE,
                         trace_rows: list | None = None,
                         trace_every: int = MC_LOG_EVERY,
                         trace_max_rows: int | None = MC_TRACE_MAX_ROWS,
                         run_meta: dict | None = None) -> tuple:
    """
    Returns: best_yB, best_gamma(J/m^2), acc_total, delta_final, best_terms(dict)
    Optionally appends trace rows to trace_rows.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    yB = np.array(yB_init, dtype=float)
    Lloc = len(yB)

    cur_terms = gammaOmega_terms_Jmol(
        yB=yB, xB=xB,
        gamma0_A=gamma0_A, gamma0_B=gamma0_B,
        omega_AB=omega_AB,
        Omega=Omega, T_K=T_K,
        zl=zl, zv=zv,
        F=F, rA=rA, rB=rB, GA=GA, KB=KB,
        include_elastic=include_elastic
    )
    Ecur = float(cur_terms["total_Jmol"])

    best_yB = yB.copy()
    best_terms = dict(cur_terms)
    best_E = float(Ecur)

    accepted_total = 0
    accepted_block = 0
    proposed_block = 0

    delta = float(delta_max)
    kT_eff = float(alpha) * R * float(T_K)

    # initial trace point
    if trace_rows is not None and trace_every and (trace_max_rows is None or len(trace_rows) < trace_max_rows):
        meta = run_meta or {}
        row = {
            **meta,
            "step": 0,
            "delta": float(delta),
            "accepted_total": int(accepted_total),
            "proposed_total": 0,
            "acc_total": 0.0,
            "ads": adsorption(yB, xB),
            **cur_terms,
        }
        for i in range(Lloc):
            row[f"yB_L{i+1}"] = float(yB[i])
        trace_rows.append(row)

    for step in range(int(steps)):
        i = int(rng.integers(0, Lloc))
        d = float(rng.uniform(-delta, delta))

        y_trial = yB.copy()
        y_trial[i] = np.clip(y_trial[i] + d, 1e-12, 1.0 - 1e-12)

        tr_terms = gammaOmega_terms_Jmol(
            yB=y_trial, xB=xB,
            gamma0_A=gamma0_A, gamma0_B=gamma0_B,
            omega_AB=omega_AB,
            Omega=Omega, T_K=T_K,
            zl=zl, zv=zv,
            F=F, rA=rA, rB=rB, GA=GA, KB=KB,
            include_elastic=include_elastic
        )
        Etr = float(tr_terms["total_Jmol"])
        dE = float(Etr - Ecur)
        accept = (dE <= 0.0) or (float(rng.random()) < np.exp(-dE / kT_eff))

        proposed_block += 1
        if accept:
            yB = y_trial
            Ecur = float(Etr)
            cur_terms = tr_terms
            accepted_total += 1
            accepted_block += 1
            if Ecur < best_E:
                best_E = float(Ecur)
                best_yB = yB.copy()
                best_terms = dict(cur_terms)

        # adaptation
        if (step + 1) <= int(adapt_steps) and (step + 1) % int(adapt_every) == 0:
            acc_block = accepted_block / max(1, proposed_block)
            if acc_block > target_acc:
                delta *= (1.0 + adapt_rate)
            else:
                delta *= (1.0 - adapt_rate)
            delta = float(np.clip(delta, DELTA_CLIP_MIN, DELTA_CLIP_MAX))
            accepted_block = 0
            proposed_block = 0

        # trace logging
        if trace_rows is not None and trace_every and ((step + 1) % int(trace_every) == 0):
            if trace_max_rows is None or len(trace_rows) < trace_max_rows:
                meta = run_meta or {}
                row = {
                    **meta,
                    "step": int(step + 1),
                    "delta": float(delta),
                    "accepted_total": int(accepted_total),
                    "proposed_total": int(step + 1),
                    "acc_total": float(accepted_total / float(step + 1)),
                    "ads": adsorption(yB, xB),
                    **cur_terms,
                }
                for j in range(Lloc):
                    row[f"yB_L{j+1}"] = float(yB[j])
                trace_rows.append(row)

    acc_total = accepted_total / float(steps)
    best_gamma = float(best_E / float(Omega))  # J/m^2
    delta_final = float(delta)
    return best_yB, best_gamma, float(acc_total), float(delta_final), best_terms

def polish_coordinate_descent(yB: np.ndarray, xB: float,
                              gamma0_A: float, gamma0_B: float,
                              omega_AB: float, Omega: float, T_K: float,
                              zl: int, zv: int,
                              F: np.ndarray,
                              rA: float, rB: float, GA: float, KB: float,
                              include_elastic: bool,
                              n_sweeps: int = 20,
                              step0: float = 0.02) -> tuple:
    """
    Deterministic local polish (coordinate descent).
    Returns yB_best, gamma_best(J/m^2), terms_best(dict).
    """
    y = np.array(yB, dtype=float)
    step = float(step0)

    def terms(yvec):
        return gammaOmega_terms_Jmol(
            yB=yvec, xB=xB,
            gamma0_A=gamma0_A, gamma0_B=gamma0_B,
            omega_AB=omega_AB,
            Omega=Omega, T_K=T_K,
            zl=zl, zv=zv,
            F=F, rA=rA, rB=rB, GA=GA, KB=KB,
            include_elastic=include_elastic
        )

    cur_terms = terms(y)
    Ecur = float(cur_terms["total_Jmol"])

    for _ in range(int(n_sweeps)):
        improved = False
        for i in range(len(y)):
            for sgn in (-1.0, +1.0):
                yt = y.copy()
                yt[i] = np.clip(yt[i] + sgn * step, 1e-12, 1.0 - 1e-12)
                t = terms(yt)
                Et = float(t["total_Jmol"])
                if Et < Ecur:
                    y, Ecur, cur_terms = yt, Et, t
                    improved = True
        if not improved:
            step *= 0.5
            if step < 1e-5:
                break

    gamma_best = float(Ecur / float(Omega))
    return y, gamma_best, cur_terms

# ============================================================
# SWEEP + EXPORT
# ============================================================
def run_sweep_and_export(A_scale: float, include_elastic: bool, out_prefix: str):
    gamma0_weighted = load_gamma0_weighted(SURFACES_JSON_PATH)
    gamma0_weighted = {k.upper(): float(v) for k, v in gamma0_weighted.items()}

    INCLUDE_ELASTIC = bool(include_elastic)

    alloys = [
        {"name": "Fe-6.207Mo_wt", "solute": "MO", "wt_solute": 6.207},
        {"name": "Fe-3.8W_wt",    "solute": "W",  "wt_solute": 3.8},
        {"name": "Fe-1.15Cu_wt",  "solute": "CU", "wt_solute": 1.15},
        {"name": "Fe-Au_wt",      "solute": "AU", "wt_solute": 2.87},
    ]
    for a in alloys:
        a["xB"] = wt_to_atfrac_binary_feX(a["wt_solute"], a["solute"])

    rng = np.random.default_rng(SEED)

    rows_wide = []
    rows_long = []
    trace_rows = [] if SAVE_MC_TRACE else None

    rA, GA = float(FE_PARAMS["rA"]), float(FE_PARAMS["GA"])

    for alloy in alloys:
        sol = alloy["solute"].upper()
        xB = float(alloy["xB"])

        if sol not in gamma0_weighted:
            raise KeyError(f"No weighted gamma0 for {sol} in surfaces.json.")
        gamma0_B = float(gamma0_weighted[sol])

        if INCLUDE_ELASTIC:
            if sol not in SOLUTE_PARAMS:
                raise KeyError(f"Missing elastic params for {sol}.")
            rB, KB = float(SOLUTE_PARAMS[sol]["rB"]), float(SOLUTE_PARAMS[sol]["KB"])
        else:
            rB, KB = 0.0, 0.0

        for plane in PLANES:
            if plane not in GEOM_BCC:
                raise KeyError(f"Plane {plane} not in GEOM_BCC.")
            if plane not in GAMMA0_FE_BY_PLANE:
                raise KeyError(f"Missing Fe gamma0 for plane {plane}.")

            geom = GEOM_BCC[plane]
            Omega_plane = float(geom["Omega"])
            zl = int(geom["zl"])
            zv = int(geom["zv"])

            gamma0_Fe_plane = float(GAMMA0_FE_BY_PLANE[plane])

            F_plane = build_F_vector(
                    plane=plane,
                    L_layers=L,
                    a_m=A_FE_BCC_M,
                    A_scale=float(A_scale)
                )
            
            for T_C in T_C_LIST:
                T_K = float(T_C) + 273.15

                # omega: Au via L0(T); others via CALPHAD LSQ
                if sol == "AU":
                    omega_AB = omega_FeAu_L0_Jmol(T_K) / Z_BCC
                else:
                    omega_AB = omega_avg_LSQ_Jmol(sol, "FE", T_K, z=Z_BCC, x_min=0.05, x_max=0.95, npts=19)

                yB_init = np.full(L, xB, dtype=float)

                run_meta = {
                    "alloy": alloy["name"],
                    "solute": sol,
                    "plane": plane,
                    "T_C": float(T_C),
                    "T_K": float(T_K),
                    "xB_bulk": float(xB),
                    "omega_Jmol": float(omega_AB),
                    "gamma0_Fe_plane_Jm2": float(gamma0_Fe_plane),
                    "gamma0_solute_weighted_Jm2": float(gamma0_B),
                    "Omega_m2mol": float(Omega_plane),
                    "zl": int(zl),
                    "zv": int(zv),
                    "include_elastic": bool(INCLUDE_ELASTIC),
                    "omega_source": ("L0_literature" if sol == "AU" else "CALPHAD_LSQ"),
                }
                for i in range(L):
                    run_meta[f"F_L{i+1}"] = float(F_plane[i])

                best_yB, best_gamma, acc, delta_final, best_terms_mc = monte_carlo_minimise(
                    yB_init=yB_init, xB=xB,
                    gamma0_A=gamma0_Fe_plane, gamma0_B=gamma0_B,
                    omega_AB=omega_AB,
                    Omega=Omega_plane, T_K=T_K,
                    zl=zl, zv=zv,
                    F=F_plane,
                    rA=rA, rB=rB, GA=GA, KB=KB,
                    steps=MC_STEPS, delta_max=DELTA_MAX,
                    rng=rng,
                    include_elastic=INCLUDE_ELASTIC,
                    alpha=ALPHA,
                    adapt_steps=ADAPT_STEPS,
                    adapt_every=ADAPT_EVERY,
                    target_acc=TARGET_ACC,
                    adapt_rate=ADAPT_RATE,
                    trace_rows=trace_rows if SAVE_MC_TRACE else None,
                    trace_every=MC_LOG_EVERY,
                    trace_max_rows=MC_TRACE_MAX_ROWS,
                    run_meta=run_meta if SAVE_MC_TRACE else None
                )

                # polish
                best_yB, best_gamma, best_terms = polish_coordinate_descent(
                    yB=best_yB, xB=xB,
                    gamma0_A=gamma0_Fe_plane, gamma0_B=gamma0_B,
                    omega_AB=omega_AB,
                    Omega=Omega_plane, T_K=T_K,
                    zl=zl, zv=zv,
                    F=F_plane,
                    rA=rA, rB=rB, GA=GA, KB=KB,
                    include_elastic=INCLUDE_ELASTIC,
                    n_sweeps=20,
                    step0=0.02
                )

                ads = adsorption(best_yB, xB)

                # store optimum term breakdown (both J/mol and J/m^2)
                term_cols = dict(best_terms)
                term_cols_m2 = {k.replace("_Jmol", "_Jm2"): float(v) / float(Omega_plane) for k, v in term_cols.items()}

                row = {
                    "alloy": alloy["name"],
                    "solute": sol,
                    "plane": plane,
                    "T_C": float(T_C),
                    "T_K": float(T_K),
                    "wt_solute": float(alloy["wt_solute"]),
                    "xB_bulk": float(xB),
                    "ads_monolayers": float(ads),
                    "gamma_Jm2": float(best_gamma),
                    "omega_Jmol": float(omega_AB),
                    "gamma0_Fe_plane_Jm2": float(gamma0_Fe_plane),
                    "gamma0_solute_weighted_Jm2": float(gamma0_B),
                    "Omega_m2mol": float(Omega_plane),
                    "zl": int(zl),
                    "zv": int(zv),
                    "mc_steps": int(MC_STEPS),
                    "delta_max": float(DELTA_MAX),
                    "delta_final": float(delta_final),
                    "adapt_steps": int(ADAPT_STEPS),
                    "adapt_every": int(ADAPT_EVERY),
                    "target_acc": float(TARGET_ACC),
                    "adapt_rate": float(ADAPT_RATE),
                    "alpha": float(ALPHA),
                    "seed": int(SEED),
                    "acceptance_ratio": float(acc),
                    "include_elastic": bool(INCLUDE_ELASTIC),
                    "phase": PHASE[0],
                    "tdb_path": TDB_PATH,
                    "omega_source": ("L0_literature" if sol == "AU" else "CALPHAD_LSQ"),
                    "a_Fe_BCC_m": float(A_FE_BCC_M),
                }
                for i in range(L):
                    row[f"yB_L{i+1}"] = float(best_yB[i])
                    row[f"F_L{i+1}"] = float(F_plane[i])

                # term breakdown at optimum
                row.update(term_cols)
                row.update(term_cols_m2)

                rows_wide.append(row)

                for i in range(L):
                    rows_long.append({
                        "alloy": alloy["name"],
                        "solute": sol,
                        "plane": plane,
                        "T_C": float(T_C),
                        "layer": int(i + 1),
                        "yB": float(best_yB[i]),
                        "F_i": float(F_plane[i]),
                        "xB_bulk": float(xB),
                        "ads_monolayers": float(ads),
                        "gamma_Jm2": float(best_gamma),
                        "omega_Jmol": float(omega_AB),
                        "omega_source": ("L0_literature" if sol == "AU" else "CALPHAD_LSQ"),
                        "acceptance_ratio": float(acc),
                        "include_elastic": bool(INCLUDE_ELASTIC),
                    })

                print(
                    f"[OK] {alloy['name']} plane {plane} T={T_C:.0f}C  "
                    f"ads={ads:+.4e}  gamma={best_gamma:.6f} J/m2  "
                    f"omega={omega_AB:+.3e} J/mol  acc={acc:.3f}  "
                    f"delta={delta_final:+.3e}"
                )

    df_wide = pd.DataFrame(rows_wide)
    df_long = pd.DataFrame(rows_long)
    df_wide.to_csv(out_prefix + "__wide.csv", index=False)
    df_long.to_csv(out_prefix + "__long.csv", index=False)

    if trace_rows is not None:
        df_trace = pd.DataFrame(trace_rows)
        df_trace.to_csv(out_prefix + "__trace.csv", index=False)
    else:
        df_trace = None

    return df_wide, df_long, df_trace

