#!/usr/bin/env python3
"""run_all.py - execute the whole analysis pipeline in order, then validate."""
import os, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STEPS = ["prepare_data.py", "exp1_pooling.py", "exp2_specificity.py", "exp3_attenuation.py",
         "exp4_dose_response.py", "exp5_robustness.py", "make_figures.py", "validate.py"]

t0 = time.time()
for s in STEPS:
    print(f"\n{'=' * 74}\n>>> {s}\n{'=' * 74}")
    rc = subprocess.run([sys.executable, os.path.join(ROOT, "src", s)], cwd=ROOT).returncode
    if rc != 0:
        sys.exit(f"FAILED at {s} (exit {rc})")
print(f"\nPipeline complete in {time.time() - t0:.1f}s")
