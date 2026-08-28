"""
obs_collapse.py — make the MF6 time-series OBS output deterministic under ATS.

WHY:  the SFR/HEAD continuous-OBS csvs get ONE row per MF6 time step.  With ATS
(adaptive time-stepping) a stress period is subdivided into extra sub-steps only
when a step fails to converge, so the NUMBER of rows varies per parameter set.
PstFrom's instruction files navigate POSITIONALLY (l1 = advance one line), so any
realization whose ATS sub-stepping differs from the reference run misaligns and
fails ("EOF encountered when executing line advance").

FIX:  keep only the rows at the STRESS-PERIOD END times (cumsum of perlen).  ATS
always lands exactly on a period boundary, so every run has all NPER period-end
rows; the intermediate sub-steps are dropped.  Result: a deterministic NPER-row
csv in EVERY run, so the positional .ins always aligns.  Applied identically to
the reference csvs (pest_prep -> the .ins is built for NPER rows) and to every
forward run (forward_run.py).  make_obs.py already weights only period-end obs and
drops sub-steps, so collapsing to period ends is consistent with the weighting.
"""
import os
import numpy as np
import pandas as pd

# obs csvs that carry per-time-step rows and are read by a positional .ins
TIMESERIES_OBS = ("cdl_gwf.obs.head.csv", "cdl_gwf.obs.sfr.csv")


def period_end_times(tdis_path):
    """Cumulative end-of-period times (days) parsed from an MF6 TDIS file.
    Handles inline perioddata and the OPEN/CLOSE form."""
    perlen, inblk = [], False
    with open(tdis_path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines):
        s = lines[i].strip(); lo = s.lower()
        if lo.startswith("begin perioddata"):
            inblk = True; i += 1; continue
        if inblk and lo.startswith("end"):
            break
        if inblk and s and not s.startswith("#"):
            if lo.startswith("open/close"):
                # perlen is the first column of the referenced file
                fn = s.split()[1].strip("'\"")
                ext = os.path.join(os.path.dirname(tdis_path), fn)
                for ln in open(ext):
                    t = ln.strip()
                    if t and not t.startswith("#"):
                        perlen.append(float(t.split()[0]))
                break
            perlen.append(float(s.split()[0]))
        i += 1
    if not perlen:
        raise ValueError(f"no perioddata parsed from {tdis_path}")
    return np.cumsum(np.asarray(perlen, dtype=float))


def collapse_obs_csv(csv_path, ends, atol=1e-2):
    """Rewrite csv_path keeping the header + only rows whose 'time' is a period end.
    Returns (n_before, n_after). No-op if the file is missing or has no 'time' col."""
    if not os.path.exists(csv_path):
        return (0, 0)
    d = pd.read_csv(csv_path)
    if "time" not in d.columns:
        return (len(d), len(d))
    t = d["time"].to_numpy(dtype=float)
    keep = np.isclose(t[:, None], ends[None, :], atol=atol).any(axis=1)
    out = d.loc[keep].drop_duplicates(subset="time", keep="last")
    out.to_csv(csv_path, index=False)
    return (len(d), len(out))


def collapse_run_obs(workdir=".", tdis="cdl_gwf.tdis", verbose=False):
    """Collapse every TIMESERIES_OBS csv in workdir to period-end rows using workdir/tdis."""
    ends = period_end_times(os.path.join(workdir, tdis))
    for name in TIMESERIES_OBS:
        nb, na = collapse_obs_csv(os.path.join(workdir, name), ends)
        if verbose and nb:
            print(f"   collapsed {name}: {nb} -> {na} rows (period ends; {nb - na} sub-steps dropped)")
    return len(ends)
