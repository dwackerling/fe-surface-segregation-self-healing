# -*- coding: utf-8 -*-
"""
single_mc.py

Analytical multilayer segregation model + Monte Carlo minimisation
for a single alloy / plane / temperature condition.

Key features
------------
- Same physical model as the original script
- Reusable run_single_case(...) entry point
- Configurable:
    * alloy
    * plane
    * temperature
    * A_scale
    * L
    * seed
    * omega_scale
    * include_elastic
- Optional MC trace output
- Returns a single results dictionary suitable for robustness studies

Notes
-----
- Elastic screening is defined as:
      F_i = A * exp(-d_i / lambda)
- Sign convention:
      total = term1 + term2 + term3 + term4 + term5 - Eel
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from pycalphad import Database, equilibrium, variables as v


# ============================================================
# GLOBAL CONFIGURATION
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

TDB_PATH = PROJECT_ROOT / "data" / "thermodynamic" / "steel_database_fix.tdb"
SURFACES_JSON_PATH = PROJECT_ROOT / "data" / "surface_energy" / "surfaces.json"

PHASE = ["BCC_A2"]
P = 101325.0  # Pa

# Default MC parameters
DEFAULT_MC_STEPS = 200_000
DEFAULT_DELTA_MAX = 1.5
DEFAULT_ALPHA = 1.0

# Adaptation
DEFAULT_ADAPT_STEPS = 50_000
DEFAULT_ADAPT_EVERY = 2_000
DEFAULT_TARGET_ACC = 0.55
DEFAULT_ADAPT_RATE = 0.15
DEFAULT_DELTA_CLIP_MIN = 1e-5
DEFAULT_DELTA_CLIP_MAX = 0.2

# Elastic screening
A_FE_BCC_M = 2.866e-10  # m

# Constants
R = 8.314  # J/mol/K
Z_BCC = 8
N_A = 6.02214076e23

# Geometry BCC per plane
GEOM_BCC = {
    "100": {"Omega": 4.947e4, "zl": 0, "zv": 4},
    "110": {"Omega": 3.498e4, "zl": 4, "zv": 2},
    "111": {"Omega": 8.568e4, "zl": 0, "zv": 4},
}

# Fe reference surface energies by plane
GAMMA0_FE_BY_PLANE = {"110": 2.45, "100": 2.50, "111": 2.73}

# Atomic masses (g/mol)
ATOMIC_MASS = {
    "FE": 55.845,
    "MO": 95.95,
    "W": 183.84,
    "CU": 63.546,
    "AU": 196.966569,
}

# Elastic parameters (SI)
FE_PARAMS = {"rA": 1.24e-10, "GA": 82e9}
SOLUTE_PARAMS = {
    "MO": {"rB": 1.39e-10, "KB": 230e9},
    "W":  {"rB": 1.41e-10, "KB": 310e9},
    "CU": {"rB": 1.28e-10, "KB": 140e9},
    "AU": {"rB": 1.44e-10, "KB": 180e9},
}

# Default alloy definitions
DEFAULT_ALLOYS = {
    "Fe-6.207Mo_wt": {"solute": "MO", "wt_solute": 6.207},
    "Fe-3.8W_wt":    {"solute": "W",  "wt_solute": 3.8},
    "Fe-1.15Cu_wt":  {"solute": "CU", "wt_solute": 1.15},
    "Fe-Au_wt":      {"solute": "AU", "wt_solute": 2.87},
}

# Load database once
DB = Database(TDB_PATH)

L_EXPORT_MAX = 10

# ============================================================
# UTILITIES
# ============================================================
def wt_to_atfrac_binary_feX(wt_solute: float, solute_symbol: str) -> float:
    """Convert Fe-(wt_solute)%X to atomic fraction x_X."""
    sol = solute_symbol.upper()
    if sol not in ATOMIC_MASS:
        raise KeyError(f"Missing atomic mass for {sol}")
    wt_fe = 100.0 - float(wt_solute)
    n_fe = wt_fe / ATOMIC_MASS["FE"]
    n_x = float(wt_solute) / ATOMIC_MASS[sol]
    return float(n_x / (n_fe + n_x))


def load_gamma0_weighted(path_json: str | Path) -> Dict[str, float]:
    """Load isotropic gamma0 per element from surfaces.json."""
    with open(path_json, "r", encoding="utf-8") as f:
        data = json.load(f)

    df = pd.json_normalize(data)
    if "pretty_formula" not in df.columns or "weighted_surface_energy" not in df.columns:
        raise ValueError("surfaces.json must contain 'pretty_formula' and 'weighted_surface_energy'.")

    gamma0: Dict[str, float] = {}
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
    """
    Build screening vector:
        F_i = A * exp(-d_i / lambda)
    with d_i = (i-1) d_hkl and lambda = a_m
    """
    d_hkl = plane_spacing_bcc_m(plane, a_m)
    lam = float(a_m)
    idx = np.arange(L_layers, dtype=float)
    d_i = idx * d_hkl
    F = np.clip(float(A_scale) * np.exp(-d_i / lam), 0.0, 1.0)
    return F.astype(float)


def omega_FeAu_L0_Jmol(T_K: float) -> float:
    """Fe-Au (BCC): L0 = 34122.3 − 12.36 T, J/mol."""
    return float(34122.3 - 12.36 * float(T_K))


def delta_H_mix_Jmol(elem1: str, x1: float, elem2: str, T_K: float) -> float:
    """Mixing enthalpy in BCC_A2 from CALPHAD."""
    elem1 = elem1.upper()
    elem2 = elem2.upper()
    x1 = float(x1)
    x2 = 1.0 - x1
    comps = [elem1, elem2, "VA"]

    H_sys = equilibrium(DB, comps, PHASE, {v.T: T_K, v.P: P, v.X(elem1): x1}, output="HM").HM.values.flatten()[0]
    H_1 = equilibrium(DB, comps, PHASE, {v.T: T_K, v.P: P, v.X(elem1): 0.99999}, output="HM").HM.values.flatten()[0]
    H_2 = equilibrium(DB, comps, PHASE, {v.T: T_K, v.P: P, v.X(elem1): 0.00001}, output="HM").HM.values.flatten()[0]

    return float(H_sys - (x1 * H_1 + x2 * H_2))


def omega_avg_LSQ_Jmol(
    elem1: str,
    elem2: str,
    T_K: float,
    z: int = Z_BCC,
    x_min: float = 0.05,
    x_max: float = 0.95,
    npts: int = 19,
) -> float:
    """
    Fit constant omega by least squares:
        ΔH_mix(x) ≈ z x(1-x) omega
    """
    xs = np.linspace(float(x_min), float(x_max), int(npts))
    xs = np.clip(xs, 1e-3, 1.0 - 1e-3)
    dH = np.array([delta_H_mix_Jmol(elem1, x, elem2, T_K) for x in xs], dtype=float)
    f = float(z) * xs * (1.0 - xs)
    return float(np.dot(dH, f) / np.dot(f, f))


def adsorption(yB: np.ndarray, xB: float) -> float:
    """Total adsorption in monolayer equivalents."""
    return float(np.sum(np.asarray(yB, dtype=float) - float(xB)))


def adsorption_first_n_layers(yB: np.ndarray, xB: float, n_layers: int = 2) -> float:
    """Adsorption restricted to first n layers."""
    yB = np.asarray(yB, dtype=float)
    n = min(int(n_layers), len(yB))
    return float(np.sum(yB[:n] - float(xB)))


# ============================================================
# ELASTIC TERM
# ============================================================
def E2_elastic_Jmol(
    yB: np.ndarray,
    xB: float,
    F: np.ndarray,
    rA: float,
    rB: float,
    GA: float,
    KB: float,
) -> float:
    """
    Eshelby-like elastic term, converted to J/mol.
    """
    yB = np.asarray(yB, dtype=float)
    F = np.asarray(F, dtype=float)

    num = 24.0 * np.pi * KB * GA * rA * rB * (rA - rB) ** 2
    den = 3.0 * KB * rB + 4.0 * GA * rA
    pref_J_per_inclusion = num / den
    summ = float(np.sum((yB - float(xB)) * F))
    return float(pref_J_per_inclusion * summ * N_A)


# ============================================================
# OBJECTIVE FUNCTION
# ============================================================
def gammaOmega_terms_Jmol(
    yB: np.ndarray,
    xB: float,
    gamma0_A: float,
    gamma0_B: float,
    omega_AB: float,
    Omega: float,
    T_K: float,
    zl: int,
    zv: int,
    F: np.ndarray,
    rA: float,
    rB: float,
    GA: float,
    KB: float,
    include_elastic: bool = False,
) -> Dict[str, float]:
    """
    Return term decomposition and total gamma*Omega in J/mol.
    """
    yB = np.asarray(yB, dtype=float)
    xB = float(xB)
    yA = 1.0 - yB
    xA = 1.0 - xB

    eps = 1e-15
    yA = np.clip(yA, eps, 1.0 - eps)
    yB = np.clip(yB, eps, 1.0 - eps)

    term1 = (yA[0] * gamma0_A + yB[0] * gamma0_B) * float(Omega)
    term2 = R * float(T_K) * float(np.sum(yA * np.log(yA / xA) + yB * np.log(yB / xB)))
    term3 = int(zv) * float(omega_AB) * (xA * xB - xA * yB[0] - xB * yA[0])

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

    total = float(term1 + term2 + term3 + term4 + term5 - Eel)

    return {
        "term1_Jmol": float(term1),
        "term2_Jmol": float(term2),
        "term3_Jmol": float(term3),
        "term4_Jmol": float(term4),
        "term5_Jmol": float(term5),
        "Eel_Jmol": float(Eel),
        "total_Jmol": float(total),
    }


# ============================================================
# MONTE CARLO
# ============================================================
def monte_carlo_minimise(
    yB_init: np.ndarray,
    xB: float,
    gamma0_A: float,
    gamma0_B: float,
    omega_AB: float,
    Omega: float,
    T_K: float,
    zl: int,
    zv: int,
    F: np.ndarray,
    rA: float,
    rB: float,
    GA: float,
    KB: float,
    steps: int = DEFAULT_MC_STEPS,
    delta_max: float = DEFAULT_DELTA_MAX,
    rng: Optional[np.random.Generator] = None,
    include_elastic: bool = False,
    alpha: float = DEFAULT_ALPHA,
    adapt_steps: int = DEFAULT_ADAPT_STEPS,
    adapt_every: int = DEFAULT_ADAPT_EVERY,
    target_acc: float = DEFAULT_TARGET_ACC,
    adapt_rate: float = DEFAULT_ADAPT_RATE,
    delta_clip_min: float = DEFAULT_DELTA_CLIP_MIN,
    delta_clip_max: float = DEFAULT_DELTA_CLIP_MAX,
    trace_every: Optional[int] = None,
    trace_max_rows: Optional[int] = None,
    run_meta: Optional[Dict[str, Any]] = None,
) -> Tuple[np.ndarray, float, float, float, Dict[str, float], List[Dict[str, Any]]]:
    """
    Returns
    -------
    best_yB, best_gamma_Jm2, acc_total, delta_final, best_terms, trace_rows
    """
    if rng is None:
        rng = np.random.default_rng(0)

    yB = np.array(yB_init, dtype=float)
    Lloc = len(yB)

    cur_terms = gammaOmega_terms_Jmol(
        yB=yB, xB=xB,
        gamma0_A=gamma0_A, gamma0_B=gamma0_B,
        omega_AB=omega_AB, Omega=Omega, T_K=T_K,
        zl=zl, zv=zv, F=F,
        rA=rA, rB=rB, GA=GA, KB=KB,
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

    trace_rows: List[Dict[str, Any]] = []

    if trace_every is not None:
        meta = run_meta or {}
        row = {
            **meta,
            "step": 0,
            "delta": float(delta),
            "accepted_total": 0,
            "proposed_total": 0,
            "acc_total": 0.0,
            "ads_total": adsorption(yB, xB),
            "ads_first2": adsorption_first_n_layers(yB, xB, n_layers=2),
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
            omega_AB=omega_AB, Omega=Omega, T_K=T_K,
            zl=zl, zv=zv, F=F,
            rA=rA, rB=rB, GA=GA, KB=KB,
            include_elastic=include_elastic
        )
        Etr = float(tr_terms["total_Jmol"])
        dE = float(Etr - Ecur)

        accept = (dE <= 0.0) or (float(rng.random()) < np.exp(-dE / kT_eff))

        proposed_block += 1
        if accept:
            yB = y_trial
            Ecur = Etr
            cur_terms = tr_terms
            accepted_total += 1
            accepted_block += 1

            if Ecur < best_E:
                best_E = float(Ecur)
                best_yB = yB.copy()
                best_terms = dict(cur_terms)

        # Adaptation
        if (step + 1) <= int(adapt_steps) and (step + 1) % int(adapt_every) == 0:
            acc_block = accepted_block / max(1, proposed_block)
            if acc_block > target_acc:
                delta *= (1.0 + adapt_rate)
            else:
                delta *= (1.0 - adapt_rate)
            delta = float(np.clip(delta, delta_clip_min, delta_clip_max))
            accepted_block = 0
            proposed_block = 0

        # Trace logging
        if trace_every is not None and (step + 1) % int(trace_every) == 0:
            if trace_max_rows is None or len(trace_rows) < trace_max_rows:
                meta = run_meta or {}
                row = {
                    **meta,
                    "step": int(step + 1),
                    "delta": float(delta),
                    "accepted_total": int(accepted_total),
                    "proposed_total": int(step + 1),
                    "acc_total": float(accepted_total / float(step + 1)),
                    "ads_total": adsorption(yB, xB),
                    "ads_first2": adsorption_first_n_layers(yB, xB, n_layers=2),
                    **cur_terms,
                }
                for j in range(Lloc):
                    row[f"yB_L{j+1}"] = float(yB[j])
                trace_rows.append(row)

    acc_total = accepted_total / float(steps)
    best_gamma = float(best_E / float(Omega))
    delta_final = float(delta)

    return best_yB, best_gamma, float(acc_total), float(delta_final), best_terms, trace_rows


def polish_coordinate_descent(
    yB: np.ndarray,
    xB: float,
    gamma0_A: float,
    gamma0_B: float,
    omega_AB: float,
    Omega: float,
    T_K: float,
    zl: int,
    zv: int,
    F: np.ndarray,
    rA: float,
    rB: float,
    GA: float,
    KB: float,
    include_elastic: bool,
    n_sweeps: int = 20,
    step0: float = 0.02,
) -> Tuple[np.ndarray, float, Dict[str, float]]:
    """
    Deterministic local polish.
    """
    y = np.array(yB, dtype=float)
    step = float(step0)

    def terms(yvec: np.ndarray) -> Dict[str, float]:
        return gammaOmega_terms_Jmol(
            yB=yvec, xB=xB,
            gamma0_A=gamma0_A, gamma0_B=gamma0_B,
            omega_AB=omega_AB, Omega=Omega, T_K=T_K,
            zl=zl, zv=zv, F=F,
            rA=rA, rB=rB, GA=GA, KB=KB,
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
# CASE SETUP
# ============================================================
def get_alloy_definition(alloy_name: str) -> Dict[str, Any]:
    if alloy_name not in DEFAULT_ALLOYS:
        raise KeyError(f"Unknown alloy '{alloy_name}'. Available: {list(DEFAULT_ALLOYS.keys())}")
    d = dict(DEFAULT_ALLOYS[alloy_name])
    d["name"] = alloy_name
    d["xB"] = wt_to_atfrac_binary_feX(d["wt_solute"], d["solute"])
    return d


def build_case_inputs(
    alloy_name: str,
    plane: str,
    T_C: float,
    A_scale: float,
    L: int,
    include_elastic: bool,
    omega_scale: float = 1.0,
) -> Dict[str, Any]:
    alloy = get_alloy_definition(alloy_name)
    sol = alloy["solute"].upper()
    xB = float(alloy["xB"])
    T_K = float(T_C) + 273.15

    gamma0_weighted = load_gamma0_weighted(SURFACES_JSON_PATH)
    gamma0_B = float(gamma0_weighted[sol])
    gamma0_A = float(GAMMA0_FE_BY_PLANE[plane])

    geom = GEOM_BCC[plane]
    Omega = float(geom["Omega"])
    zl = int(geom["zl"])
    zv = int(geom["zv"])

    if sol == "AU":
        omega_AB = omega_FeAu_L0_Jmol(T_K) / Z_BCC
        omega_source = "L0_literature"
    else:
        omega_AB = omega_avg_LSQ_Jmol(sol, "FE", T_K, z=Z_BCC, x_min=0.05, x_max=0.95, npts=19)
        omega_source = "CALPHAD_LSQ"

    omega_AB *= float(omega_scale)

    if include_elastic:
        rA, GA = float(FE_PARAMS["rA"]), float(FE_PARAMS["GA"])
        rB, KB = float(SOLUTE_PARAMS[sol]["rB"]), float(SOLUTE_PARAMS[sol]["KB"])
    else:
        rA, GA = float(FE_PARAMS["rA"]), float(FE_PARAMS["GA"])
        rB, KB = 0.0, 0.0

    F = build_F_vector(plane=plane, L_layers=L, a_m=A_FE_BCC_M, A_scale=A_scale)

    return {
        "alloy": alloy,
        "solute": sol,
        "xB": xB,
        "plane": plane,
        "T_C": float(T_C),
        "T_K": float(T_K),
        "L": int(L),
        "A_scale": float(A_scale),
        "include_elastic": bool(include_elastic),
        "omega_scale": float(omega_scale),
        "gamma0_A": gamma0_A,
        "gamma0_B": gamma0_B,
        "omega_AB": float(omega_AB),
        "omega_source": omega_source,
        "Omega": Omega,
        "zl": zl,
        "zv": zv,
        "F": F,
        "rA": rA,
        "rB": rB,
        "GA": GA,
        "KB": KB,
    }


# ============================================================
# PUBLIC ENTRY POINT
# ============================================================
def run_single_case(
    alloy_name: str,
    plane: str,
    T_C: float,
    A_scale: float = 0.0,
    L: int = 4,
    seed: int = 0,
    include_elastic: bool = False,
    omega_scale: float = 1.0,
    mc_steps: int = DEFAULT_MC_STEPS,
    delta_max: float = DEFAULT_DELTA_MAX,
    alpha: float = DEFAULT_ALPHA,
    adapt_steps: int = DEFAULT_ADAPT_STEPS,
    adapt_every: int = DEFAULT_ADAPT_EVERY,
    target_acc: float = DEFAULT_TARGET_ACC,
    adapt_rate: float = DEFAULT_ADAPT_RATE,
    delta_clip_min: float = DEFAULT_DELTA_CLIP_MIN,
    delta_clip_max: float = DEFAULT_DELTA_CLIP_MAX,
    save_trace: bool = False,
    trace_every: int = 500,
    trace_max_rows: Optional[int] = None,
    polish: bool = True,
) -> Tuple[Dict[str, Any], Optional[pd.DataFrame]]:
    """
    Run one alloy / plane / temperature case and return:
        result_dict, trace_df (or None)

    Important:
    ----------
    Results are exported with a fixed number of layer columns up to L_EXPORT_MAX,
    filling missing layers with NaN. This keeps CSV outputs consistent across runs.
    """
    inputs = build_case_inputs(
        alloy_name=alloy_name,
        plane=plane,
        T_C=T_C,
        A_scale=A_scale,
        L=L,
        include_elastic=include_elastic,
        omega_scale=omega_scale,
    )

    alloy = inputs["alloy"]
    xB = inputs["xB"]

    rng = np.random.default_rng(seed)
    yB_init = yB_init = np.array([0.9, 0.5, 0.2, xB])

    run_meta = {
        "alloy": alloy["name"],
        "solute": inputs["solute"],
        "plane": plane,
        "T_C": float(T_C),
        "T_K": inputs["T_K"],
        "xB_bulk": xB,
        "L": int(L),
        "seed": int(seed),
        "A_scale": float(A_scale),
        "omega_scale": float(omega_scale),
        "omega_Jmol": inputs["omega_AB"],
        "omega_source": inputs["omega_source"],
        "include_elastic": bool(include_elastic),
    }

    best_yB_mc, best_gamma_mc, acc, delta_final, best_terms_mc, trace_rows = monte_carlo_minimise(
        yB_init=yB_init,
        xB=xB,
        gamma0_A=inputs["gamma0_A"],
        gamma0_B=inputs["gamma0_B"],
        omega_AB=inputs["omega_AB"],
        Omega=inputs["Omega"],
        T_K=inputs["T_K"],
        zl=inputs["zl"],
        zv=inputs["zv"],
        F=inputs["F"],
        rA=inputs["rA"],
        rB=inputs["rB"],
        GA=inputs["GA"],
        KB=inputs["KB"],
        steps=mc_steps,
        delta_max=delta_max,
        rng=rng,
        include_elastic=include_elastic,
        alpha=alpha,
        adapt_steps=adapt_steps,
        adapt_every=adapt_every,
        target_acc=target_acc,
        adapt_rate=adapt_rate,
        delta_clip_min=delta_clip_min,
        delta_clip_max=delta_clip_max,
        trace_every=trace_every if save_trace else None,
        trace_max_rows=trace_max_rows,
        run_meta=run_meta,
    )

    gamma_before_polish = best_gamma_mc
    yB_before_polish = best_yB_mc.copy()

    if polish:
        best_yB, best_gamma, best_terms = polish_coordinate_descent(
            yB=best_yB_mc,
            xB=xB,
            gamma0_A=inputs["gamma0_A"],
            gamma0_B=inputs["gamma0_B"],
            omega_AB=inputs["omega_AB"],
            Omega=inputs["Omega"],
            T_K=inputs["T_K"],
            zl=inputs["zl"],
            zv=inputs["zv"],
            F=inputs["F"],
            rA=inputs["rA"],
            rB=inputs["rB"],
            GA=inputs["GA"],
            KB=inputs["KB"],
            include_elastic=include_elastic,
            n_sweeps=20,
            step0=0.02,
        )
    else:
        best_yB, best_gamma, best_terms = best_yB_mc, best_gamma_mc, best_terms_mc

    ads_total = adsorption(best_yB, xB)
    ads_first2 = adsorption_first_n_layers(best_yB, xB, n_layers=2)

    result: Dict[str, Any] = {
        "alloy": alloy["name"],
        "solute": inputs["solute"],
        "plane": plane,
        "T_C": float(T_C),
        "T_K": inputs["T_K"],
        "wt_solute": float(alloy["wt_solute"]),
        "xB_bulk": xB,
        "L": int(L),
        "seed": int(seed),
        "A_scale": float(A_scale),
        "omega_scale": float(omega_scale),
        "include_elastic": bool(include_elastic),
        "omega_Jmol": inputs["omega_AB"],
        "omega_source": inputs["omega_source"],
        "gamma0_Fe_plane_Jm2": inputs["gamma0_A"],
        "gamma0_solute_weighted_Jm2": inputs["gamma0_B"],
        "Omega_m2mol": inputs["Omega"],
        "zl": inputs["zl"],
        "zv": inputs["zv"],
        "mc_steps": int(mc_steps),
        "delta_max": float(delta_max),
        "delta_final": float(delta_final),
        "acceptance_ratio": float(acc),
        "gamma_mc_Jm2": float(gamma_before_polish),
        "gamma_Jm2": float(best_gamma),
        "gamma_polish_improvement_Jm2": float(gamma_before_polish - best_gamma),
        "ads_monolayers": float(ads_total),
        "ads_first2_monolayers": float(ads_first2),
        "phase": PHASE[0],
        "tdb_path": TDB_PATH,
        "a_Fe_BCC_m": float(A_FE_BCC_M),
    }

    # fixed-width export for layer-dependent outputs
    for i in range(L_EXPORT_MAX):
        idx = i + 1
        result[f"F_L{idx}"] = float(inputs["F"][i]) if i < len(inputs["F"]) else np.nan
        result[f"yB_mc_L{idx}"] = float(yB_before_polish[i]) if i < len(yB_before_polish) else np.nan
        result[f"yB_L{idx}"] = float(best_yB[i]) if i < len(best_yB) else np.nan

    result.update(best_terms)
    for k, v in best_terms.items():
        if k.endswith("_Jmol"):
            result[k.replace("_Jmol", "_Jm2")] = float(v) / float(inputs["Omega"])

    trace_df = pd.DataFrame(trace_rows) if save_trace else None
    return result, trace_df

# ============================================================
# EXAMPLE
# ============================================================
if __name__ == "__main__":
    result, trace_df = run_single_case(
        alloy_name="Fe-3.8W_wt",
        plane="110",
        T_C=550,
        A_scale=0.4,
        L=4,
        seed=0,
        include_elastic=True,
        omega_scale=1.0,
        save_trace=True,
        trace_every=1000,
        trace_max_rows=500,
    )

    print(pd.Series(result))
    if trace_df is not None:
        trace_df.to_csv("single_case_trace.csv", index=False)