"""
config.py — single place for the machine-specific disk paths.

Every CdL script imports its paths from here, so moving the project to a new machine
means editing only the three roots below (or setting the matching environment variables) —
no path edits scattered across the scripts, and nothing to run (replaces port_paths.py).

    BASE         project root that holds CdL_model/ and CdL_pest/ (and, usually, code/)
    MODFLOW_DIR  folder holding the MODFLOW 6 / PEST++ / Triangle executables
    PYTHON_EXE   the conda-env python that PEST++ workers use to run forward_run.py

Everything else is derived from these. `CODE` self-locates to this file's folder.
"""
import os
from pathlib import Path

# ==== the only machine-specific settings (edit these, or set the env vars) ============
BASE        = Path(os.environ.get("CDL_BASE",        r"X:\3p1p1\DR3PDA1"))
MODFLOW_DIR = Path(os.environ.get("CDL_MODFLOW_DIR", r"C:\sw\MODFLOWandCo"))
PYTHON_EXE  =      os.environ.get("CDL_PYTHON_EXE",  r"C:\Users\su-alain.frances\AppData\Local\miniconda3\envs\mf6models\python.exe")

# ==== derived paths (usually no need to edit) ========================================
CODE      = Path(__file__).resolve().parent          # this folder (the repo's CdL/code)
MODEL     = BASE / "CdL_model"                        # MF6 workspace; gis/ forcing/ RS/ live under here
PEST      = BASE / "CdL_pest"                         # PEST / pyEMU workspace
MODIS_DIR = MODEL / "RS" / "MODIS4CDL"               # MODIS MOD16A2GF granules

MF6_EXE      = str(MODFLOW_DIR / "mf6.7.0_win64" / "bin" / "mf6.exe")
TRIANGLE_EXE = str(MODFLOW_DIR / "win64" / "triangle.exe")
PESTPP_IES   = str(MODFLOW_DIR / "pestpp-5.2.27-win" / "bin" / "pestpp-ies.exe")
