"""
Post-processing for the CdL pestpp-ies calibration.  Reads the master-dir outputs
(no model re-run) and writes figures + CSVs to CdL_pest/ies_output/:

  1. phi_progress.png   - phi (mean/min-max) vs IES iteration  (if phi.actual.csv)
  2. par_ensemble.png   - prior vs posterior parameters, by group  (ensemble runs only)
  3. head_fit.png       - simulated head vs observed, per obs point
  4. aet_fit.png        - simulated zonal AET vs MODIS, per zone
  5. phi_by_group.png   - phi split head / aet / sfr-inlet / sfr-outlet (last stage)

Adapts to the run mode:
  NOPTMAX 0  -> base run: uses cdl.base.obs.csv (single simulated series)
  NOPTMAX -1 -> prior MC: uses cdl.0.obs.csv (prior ensemble)
  NOPTMAX >0 -> IES: uses cdl.0.obs.csv (prior) + cdl.N.obs.csv (posterior)
Run:  in the mf6models env (e.g. from Spyder), or from a shell:
      conda run -p C:/Users/su-alain.frances/AppData/Local/miniconda3/envs/mf6models python postproc_ies.py
"""
import config
import glob, re
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import pyemu

PEST   = Path(str(config.PEST))
MASTER = PEST / "master"
OUT    = PEST / "ies_output"; OUT.mkdir(exist_ok=True)
CASE   = "cdl"
# SIM_START from the model's sidecar (single source of truth = cdl_gwf_model_fable_v2.py).
_SSF = Path((str(config.MODEL) + r"\last_sim_start.txt"))
SIM_START = pd.Timestamp(_SSF.read_text().strip()) if _SSF.exists() else pd.Timestamp(1981, 1, 1)
SPINUP_NPER = 12


def _iters():
    its = []
    for f in glob.glob(str(MASTER / f"{CASE}.*.obs.csv")):
        m = re.search(rf"{CASE}\.(\d+)\.obs\.csv$", Path(f).name)   # numbered ensembles only
        if m:
            its.append(int(m.group(1)))
    return sorted(its)


def _load(name):
    f = MASTER / name
    return pd.read_csv(f, index_col=0) if f.exists() else None


def _stages():
    """[(label, DataFrame[realizations x obsnme], colour)] — base run or prior/posterior."""
    its = _iters()
    if its:
        st = [("prior", _load(f"{CASE}.{its[0]}.obs.csv"), "0.7")]
        if its[-1] != its[0]:
            st.append(("posterior", _load(f"{CASE}.{its[-1]}.obs.csv"), "tab:blue"))
        return [(l, d, c) for l, d, c in st if d is not None]
    base = _load(f"{CASE}.base.obs.csv")                            # NOPTMAX=0 base run
    return [("base run", base, "tab:blue")] if base is not None else []


def main():
    pst = pyemu.Pst(str(MASTER / f"{CASE}.pst"))
    od = pst.observation_data.copy()
    its = _iters(); stages = _stages()
    print(f">> iterations: {its or 'none (base run)'}; stages: {[s[0] for s in stages]}")

    # ---- 1. phi progress ----
    phi = _load(f"{CASE}.phi.actual.csv")
    if phi is not None:
        phi = phi.reset_index()
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(phi["iteration"], phi["mean"], "-o", color="tab:blue", label="mean")
        if {"min", "max"}.issubset(phi.columns):
            ax.fill_between(phi["iteration"], phi["min"], phi["max"], alpha=0.2, color="tab:blue", label="min-max")
        ax.set_yscale("log"); ax.set_xlabel("IES iteration"); ax.set_ylabel("phi (actual)")
        ax.set_title("CdL IES - phi progress"); ax.grid(alpha=0.3); ax.legend()
        fig.tight_layout(); fig.savefig(OUT / "phi_progress.png", dpi=140, bbox_inches="tight"); plt.close(fig)
        print("   wrote phi_progress.png  (base/last mean phi = %.1f)" % phi["mean"].iloc[-1])

    # ---- 2. parameter ensembles (ensemble runs only) ----
    if its:
        pe0, peN = _load(f"{CASE}.{its[0]}.par.csv"), _load(f"{CASE}.{its[-1]}.par.csv")
        if pe0 is not None:
            pdta = pst.parameter_data; groups = list(pst.par_groups)
            fig, axes = plt.subplots(1, len(groups), figsize=(2.4 * len(groups) + 1, 4.5), squeeze=False)
            for ax, g in zip(axes.ravel(), groups):
                pn = pdta.index[pdta.pargp == g]; log = not g.startswith("ghb_head")
                f = (lambda v: np.log10(np.clip(v, 1e-6, None))) if log else (lambda v: v)
                data = [f(pe0[pn].values.ravel())] + ([f(peN[pn].values.ravel())] if peN is not None and its[-1] != its[0] else [])
                ax.boxplot(data, showfliers=False)                       # set labels separately (version-agnostic)
                ax.set_xticks(range(1, len(data) + 1)); ax.set_xticklabels(["prior", "post"][:len(data)])
                ax.set_title(g, fontsize=8); ax.tick_params(labelsize=7)
                ax.set_ylabel(("log10 " if log else "") + "value", fontsize=7); ax.grid(alpha=0.3, axis="y")
            fig.suptitle("CdL IES - parameter ensembles"); fig.tight_layout()
            fig.savefig(OUT / "par_ensemble.png", dpi=140, bbox_inches="tight"); plt.close(fig)
            print("   wrote par_ensemble.png")

    if not stages:
        print("   (no obs output found). Done."); return

    # ---- identify obs by OBSNME + reconstruct key/date (obgnme is unreliable: the
    #      first AET/head use-col was grouped under the bare obsgp name) ----
    aetmask = od.index.str.contains("aet_otype")
    hmask = od.index.str.contains("head_otype")
    od["key"] = ""; od["date_dt"] = pd.NaT
    am = od.index[aetmask]
    od.loc[am, "key"] = [re.search(r"usecol:(.+?)_date:", i).group(1) for i in am]
    od.loc[am, "date_dt"] = pd.to_datetime([re.search(r"_date:(.+)$", i).group(1) for i in am], errors="coerce")
    hm = od.index[hmask]
    # head usecol = piezo (p#_l#) OR the virtual-P4 outlets (outw_l1 / outs_l2) -> match ANY usecol
    od.loc[hm, "key"] = [re.search(r"usecol:(.+?)_time:", i).group(1) for i in hm]
    htime = np.array([float(re.search(r"_time:(.+)$", i).group(1)) for i in hm])
    real_times = sorted(np.unique(htime))[SPINUP_NPER:]
    t2d = dict(zip(real_times, pd.date_range(SIM_START, periods=len(real_times), freq="MS")))
    od.loc[hm, "date_dt"] = [t2d.get(t, pd.NaT) for t in htime]

    def fit_panels(mask, title, fname, ncols=3):
        sub = od[mask & (od.weight > 0)]
        if sub.empty:
            print(f"   ({fname}: no weighted obs, skipped)"); return
        keys = sorted(sub["key"].unique())
        nrows = int(np.ceil(len(keys) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 3.0 * nrows), squeeze=False)
        for ax, k in zip(axes.ravel(), keys):
            rows = sub[sub["key"] == k].sort_values("date_dt"); d = rows["date_dt"].values; nm = rows.index
            for lab, ens, col in stages:
                sim = ens[nm]
                if len(sim) > 1:
                    ax.fill_between(d, sim.quantile(0.05).values, sim.quantile(0.95).values, color=col, alpha=0.25)
                    ax.plot(d, sim.median().values, color=col, lw=1.0, label=lab)
                else:
                    ax.plot(d, sim.iloc[0].values, color=col, lw=1.1, label=lab)
            ax.plot(d, rows["obsval"].values, "o" if "head" in fname else "-", ms=3, color="k", label="observed")
            ax.set_title(str(k), fontsize=9); ax.tick_params(labelsize=7); ax.grid(alpha=0.3)
        for ax in axes.ravel()[len(keys):]:
            ax.axis("off")
        axes.ravel()[0].legend(fontsize=7, loc="best")
        fig.suptitle(title); fig.tight_layout(); fig.savefig(OUT / fname, dpi=130, bbox_inches="tight"); plt.close(fig)
        print(f"   wrote {fname} ({len(keys)} panels)")

    fit_panels(hmask, "CdL - head fit (simulated vs observed)", "head_fit.png", ncols=3)
    fit_panels(aetmask, "CdL - zonal AET fit (simulated vs MODIS)", "aet_fit.png", ncols=3)

    # ---- 5. phi by group (last stage median) — head / aet / sfr-inlet / sfr-outlet ----
    # SFR is split out (was silently lumped into "head"): the stream OUTLET is the baseflow /
    # groundwater-discharge target the DRN params calibrate to; the inlet is prescribed (~0).
    sim = stages[-1][1]
    med = sim.median() if len(sim) > 1 else sim.iloc[0]
    res = med.reindex(od.index) - od["obsval"]
    contrib = (od["weight"] ** 2 * res ** 2)

    def _grp(i):
        if "aet_otype" in i:      return "aet"
        if "head_otype" in i:     return "head"
        if "usecol:outlet_" in i: return "sfr-outlet"      # stream outlet (baseflow / GW target)
        if "usecol:inlet_" in i:  return "sfr-inlet"       # prescribed inflow (residual ~0)
        return "other"
    grp = contrib.groupby(od.index.map(_grp)).sum()
    _order = ["head", "aet", "sfr-inlet", "sfr-outlet", "other"]
    grp = grp.reindex([g for g in _order if g in grp.index])
    _cmap = {"head": "tab:brown", "aet": "tab:green", "sfr-inlet": "tab:cyan",
             "sfr-outlet": "tab:blue", "other": "0.6"}
    grp.to_csv(OUT / "phi_by_group.csv")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(grp.index, grp.values, color=[_cmap.get(g, "0.6") for g in grp.index])
    ax.set_ylabel("phi contribution"); ax.set_title(f"CdL - phi by group ({stages[-1][0]})")
    for i, v in enumerate(grp.values):
        ax.annotate(f"{v:,.0f}", (i, v), ha="center", va="bottom", fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "phi_by_group.png", dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"   wrote phi_by_group.png  (phi: {grp.to_dict()})")
    print(">> done. Outputs in", OUT)


if __name__ == "__main__":
    main()
