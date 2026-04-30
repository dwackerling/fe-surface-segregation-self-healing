# -*- coding: utf-8 -*-
"""
run_robustness_tests.py

Driver for robustness studies using single_mc.run_single_case.

Implemented tests
-----------------
1) seed robustness
2) L robustness
3) omega robustness

Features
--------
- Incremental save after every completed case
- Automatic resume: already completed cases are skipped
- Progress bar with tqdm
- Optional trace saving for selected representative cases
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Set

import pandas as pd
from tqdm import tqdm
import numpy as np

from single_mc import run_single_case, DEFAULT_ALLOYS


# ============================================================
# USER SETTINGS
# ============================================================
OUTDIR = Path("robustness_outputs")
OUTDIR.mkdir(parents=True, exist_ok=True)

LONG_CSV = OUTDIR / "robustness_results__long.csv"
SUMMARY_CSV = OUTDIR / "robustness_results__summary.csv"

# Common conditions
PLANES = ["100", "110", "111"]
T_C_LIST = list(range(400, 701, 25))
A_SCALE_BASE = 0.0
INCLUDE_ELASTIC_BASE = False

# Test grids
SEEDS = [0, 1, 2]
L_VALUES = [3, 4, 6, 8, 10]
OMEGA_SCALES = np.arange(0.8, 1.21, 0.05).tolist()

# Optional trace saving for representative cases only
SAVE_TRACE_FOR_REP_CASES = True
TRACE_CASES = {
    ("Fe-6.207Mo_wt", "110", 550.0),
    ("Fe-Au_wt", "110", 550.0),
    ("Fe-3.8W_wt", "110", 550.0),
    ("Fe-1.15Cu_wt", "110", 550.0)
}
TRACE_EVERY = 1000
TRACE_MAX_ROWS = 500


# ============================================================
# HELPERS
# ============================================================
def make_case_id(
    test_type: str,
    alloy: str,
    plane: str,
    T_C: float,
    test_value: float | int,
) -> str:
    return f"{test_type}|{alloy}|{plane}|{float(T_C):.1f}|{test_value}"


def load_completed_case_ids(long_csv: Path) -> Set[str]:
    if not long_csv.exists():
        return set()

    try:
        df = pd.read_csv(long_csv)
    except Exception:
        return set()

    if "case_id" not in df.columns:
        return set()

    return set(df["case_id"].astype(str).unique())


def append_result_row(long_csv: Path, row: Dict[str, Any]) -> None:
    df = pd.DataFrame([row])
    header = not long_csv.exists()
    df.to_csv(long_csv, mode="a", header=header, index=False)


def save_trace_if_needed(
    trace_df: pd.DataFrame | None,
    test_type: str,
    alloy_name: str,
    plane: str,
    T_C: float,
    test_value: float | int,
) -> None:
    if trace_df is None or trace_df.empty:
        return

    safe_test_value = str(test_value).replace(".", "p")
    trace_name = (
        f"trace__{test_type}__{alloy_name}__{plane}__"
        f"{int(T_C)}C__{safe_test_value}.csv"
    )
    trace_df.to_csv(OUTDIR / trace_name, index=False)


def build_seed_tasks() -> List[Dict[str, Any]]:
    tasks: List[Dict[str, Any]] = []
    for alloy_name in DEFAULT_ALLOYS.keys():
        for plane in PLANES:
            for T_C in T_C_LIST:
                for seed in SEEDS:
                    tasks.append({
                        "test_type": "seed",
                        "alloy_name": alloy_name,
                        "plane": plane,
                        "T_C": float(T_C),
                        "seed": int(seed),
                        "L": 4,
                        "omega_scale": 1.0,
                        "A_scale": A_SCALE_BASE,
                        "include_elastic": INCLUDE_ELASTIC_BASE,
                        "test_value": int(seed),
                    })
    return tasks


def build_L_tasks() -> List[Dict[str, Any]]:
    tasks: List[Dict[str, Any]] = []
    for alloy_name in DEFAULT_ALLOYS.keys():
        for plane in PLANES:
            for T_C in T_C_LIST:
                for L in L_VALUES:
                    tasks.append({
                        "test_type": "L",
                        "alloy_name": alloy_name,
                        "plane": plane,
                        "T_C": float(T_C),
                        "seed": 0,
                        "L": int(L),
                        "omega_scale": 1.0,
                        "A_scale": A_SCALE_BASE,
                        "include_elastic": INCLUDE_ELASTIC_BASE,
                        "test_value": int(L),
                    })
    return tasks


def build_omega_tasks() -> List[Dict[str, Any]]:
    tasks: List[Dict[str, Any]] = []
    for alloy_name in DEFAULT_ALLOYS.keys():
        for plane in PLANES:
            for T_C in T_C_LIST:
                for omega_scale in OMEGA_SCALES:
                    tasks.append({
                        "test_type": "omega",
                        "alloy_name": alloy_name,
                        "plane": plane,
                        "T_C": float(T_C),
                        "seed": 0,
                        "L": 4,
                        "omega_scale": float(omega_scale),
                        "A_scale": A_SCALE_BASE,
                        "include_elastic": INCLUDE_ELASTIC_BASE,
                        "test_value": float(omega_scale),
                    })
    return tasks


def summarise_results(long_csv: Path, summary_csv: Path) -> None:
    if not long_csv.exists():
        return

    df_long = pd.read_csv(long_csv)
    if df_long.empty:
        return

    metrics = [
        "gamma_Jm2",
        "ads_monolayers",
        "ads_first2_monolayers",
        "acceptance_ratio",
        "delta_final",
        "term1_Jmol",
        "term2_Jmol",
        "term3_Jmol",
        "term4_Jmol",
        "term5_Jmol",
        "Eel_Jmol",
        "total_Jmol",
    ]
    existing_metrics = [m for m in metrics if m in df_long.columns]

    group_cols = ["test_type", "alloy", "solute", "plane", "T_C"]

    summary = (
        df_long
        .groupby(group_cols)[existing_metrics]
        .agg(["mean", "std", "min", "max"])
        .reset_index()
    )

    summary.columns = [
        "_".join(col).strip("_") if isinstance(col, tuple) else col
        for col in summary.columns
    ]

    summary.to_csv(summary_csv, index=False)


def run_tasks(
    tasks: Iterable[Dict[str, Any]],
    completed_ids: Set[str],
) -> None:
    tasks = list(tasks)
    remaining_tasks = []

    for task in tasks:
        case_id = make_case_id(
            test_type=task["test_type"],
            alloy=task["alloy_name"],
            plane=task["plane"],
            T_C=task["T_C"],
            test_value=task["test_value"],
        )
        task["case_id"] = case_id
        if case_id not in completed_ids:
            remaining_tasks.append(task)

    if not remaining_tasks:
        print("[INFO] No pending tasks in this block.")
        return

    for task in tqdm(remaining_tasks, desc="Running robustness tasks", unit="case"):
        alloy_name = task["alloy_name"]
        plane = task["plane"]
        T_C = task["T_C"]
        seed = task["seed"]
        L = task["L"]
        omega_scale = task["omega_scale"]
        A_scale = task["A_scale"]
        include_elastic = task["include_elastic"]
        test_type = task["test_type"]
        test_value = task["test_value"]
        case_id = task["case_id"]

        save_trace = (
            SAVE_TRACE_FOR_REP_CASES
            and (alloy_name, plane, T_C) in TRACE_CASES
            and test_type == "seed"
        )

        try:
            result, trace_df = run_single_case(
                alloy_name=alloy_name,
                plane=plane,
                T_C=T_C,
                A_scale=A_scale,
                L=L,
                seed=seed,
                include_elastic=include_elastic,
                omega_scale=omega_scale,
                save_trace=save_trace,
                trace_every=TRACE_EVERY,
                trace_max_rows=TRACE_MAX_ROWS,
            )

            result["test_type"] = test_type
            result["test_value"] = test_value
            result["case_id"] = case_id

            append_result_row(LONG_CSV, result)

            save_trace_if_needed(
                trace_df=trace_df,
                test_type=test_type,
                alloy_name=alloy_name,
                plane=plane,
                T_C=T_C,
                test_value=test_value,
            )

        except Exception as exc:
            print(
                f"\n[ERROR] Failed case: {case_id}\n"
                f"        alloy={alloy_name}, plane={plane}, T={T_C}, "
                f"seed={seed}, L={L}, omega_scale={omega_scale}\n"
                f"        error={exc}"
            )


# ============================================================
# MAIN
# ============================================================
def main() -> None:
    # remove old corrupted CSVs manually before rerun
    completed_ids = load_completed_case_ids(LONG_CSV)

    # print("[INFO] Running seed robustness...")
    run_tasks(build_seed_tasks(), completed_ids)
    completed_ids = load_completed_case_ids(LONG_CSV)

    print("[INFO] Running L robustness...")
    run_tasks(build_L_tasks(), completed_ids)
    completed_ids = load_completed_case_ids(LONG_CSV)

    print("[INFO] Running omega robustness...")
    run_tasks(build_omega_tasks(), completed_ids)

    # print("[INFO] Running seed robustness ONLY (trace generation)...")
    # run_tasks(build_seed_tasks(), completed_ids)

    summarise_results(LONG_CSV, SUMMARY_CSV)

    print(f"[OK] Long results saved to: {LONG_CSV}")
    print(f"[OK] Summary results saved to: {SUMMARY_CSV}")


if __name__ == "__main__":
    main()