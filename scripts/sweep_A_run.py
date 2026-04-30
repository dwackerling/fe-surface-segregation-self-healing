# -*- coding: utf-8 -*-
"""
Sweep of elastic scaling factor A in [A_MIN, A_MAX] + single no-elastic baseline.
Parallelises independent runs over (A, include_elastic) pairs.

Outputs per run:
  {OUT_DIR}/EL{0|1}__A{A}__wide.csv
  {OUT_DIR}/EL{0|1}__A{A}__long.csv
  {OUT_DIR}/EL{0|1}__A{A}__trace.csv   (if enabled in seg_model)
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Tuple
from tqdm import tqdm

# ---------------------------
# USER SETTINGS
# ---------------------------
A_MIN = 0.0
A_MAX = 1.0
N_A = 11                  # e.g. 0.0, 0.1, ..., 1.0
RUN_NO_ELASTIC = True    # single baseline without elastic term
N_WORKERS = 5            # <=2–3 recommended (pycalphad + SciPy are heavy on Windows)

OUT_DIR = "A_sweep_outputs"
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------
# IMPORT YOUR CORE MODEL
# ---------------------------
from seg_model import run_sweep_and_export

# ---------------------------
# JOB DEFINITION
# ---------------------------
@dataclass(frozen=True)
class Job:
    A: float
    include_elastic: bool

def _fmt_A(A: float) -> str:
    # stable filenames: A0p100 instead of A0.1
    return f"{A:.3f}".replace(".", "p")

def worker(job: Job) -> Tuple[Job, str]:
    A = float(job.A)
    inc = bool(job.include_elastic)

    tagA = _fmt_A(A)
    tagE = "EL1" if inc else "EL0"
    out_prefix = os.path.join(OUT_DIR, f"{tagE}__A{tagA}")

    # Call your model
    run_sweep_and_export(
        A_scale=A,
        include_elastic=inc,
        out_prefix=out_prefix
    )

    return job, out_prefix

# ---------------------------
# MAIN
# ---------------------------
def main():
    A_values = np.linspace(A_MIN, A_MAX, int(N_A)).astype(float).tolist()

    # Elastic runs for all A
    jobs = [Job(A=A, include_elastic=True) for A in A_values]

    # Single baseline run without elastic term (A irrelevant)
    if RUN_NO_ELASTIC:
        jobs.append(Job(A=A_values[0], include_elastic=False))

    print(f"Running {len(jobs)} jobs with {N_WORKERS} workers...")

    results = []
    with Pool(processes=N_WORKERS) as pool:
        for item in tqdm(
            pool.imap_unordered(worker, jobs),
            total=len(jobs),
            desc="A-sweep",
        ):
            results.append(item)

    # Manifest for post-processing
    manifest_rows = []
    for job, prefix in results:
        manifest_rows.append({
            "A": job.A,
            "include_elastic": job.include_elastic,
            "out_prefix": prefix,
            "wide_csv": prefix + "__wide.csv",
            "long_csv": prefix + "__long.csv",
            "trace_csv": prefix + "__trace.csv",
        })

    manifest = pd.DataFrame(manifest_rows).sort_values(["include_elastic", "A"])
    manifest_path = os.path.join(OUT_DIR, "manifest.csv")
    manifest.to_csv(manifest_path, index=False)

    print(f"\nDone. Manifest saved to: {manifest_path}")


if __name__ == "__main__":
    main()