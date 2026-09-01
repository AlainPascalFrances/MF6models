"""
port_paths.py — bulk search-and-replace of the hardcoded path prefixes across every file in a
folder, to move the CdL model/PEST code to a new machine.

The search strings are the SIX prefixes documented in CdL/instructions/02_PATHS_TO_CHANGE.md
(reproduced below with their FULL values — the doc truncates prefix #3 with an ellipsis).

HOW TO USE
  1. Edit the NEW values in REPLACEMENTS below (leave "" to skip a prefix for now).
  2. Run once as-is (DRY_RUN = True): it prints every file/line that WOULD change, and — for any
     prefix whose NEW value is still "" — it just COUNTS where that old prefix occurs (scan mode),
     so you can see exactly what is present before deciding.
  3. When the preview looks right, set DRY_RUN = False and run again to write the changes
     (a .bak copy of each modified file is kept if MAKE_BACKUP = True).

Only text files with an extension in TEXT_EXT are touched; binaries (.npz/.pkl/.grb/.png/.tif/…)
and the .git folder are skipped. Safe to re-run.
"""
from pathlib import Path

# ------------------------------------------------------------------ config
# Folder to process (default: the code folder next to this script).
TARGET_DIR = Path(__file__).resolve().parent / "code"

# (label, OLD prefix, NEW prefix)  — fill the NEW column for your machine.
# OLD values are the six prefixes from 02_PATHS_TO_CHANGE.md (full, real strings).
REPLACEMENTS = [
    ("1 DATA_ROOT",  r"E:\00code_ws\DRYAD",                                                        r""),   # e.g. r"D:\DRYAD\work"
    ("2 CODE_DIR",   r"E:\00code\flopy\dryad_cdl",                                                 r""),   # e.g. r"D:\DRYAD\MF6models\CdL\code"
    ("3 GIS_ROOT",   r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD", r""),   # e.g. r"D:\DRYAD\gis"
    ("4 COS_TIF",    r"E:\ArcGis_Data\WorkSpace\COS\COSc2025\COSc_2025_N3_v0_TM06.tif",             r""),   # your COS raster path
    ("5 MODFLOW",    r"C:\00MODFLOW",                                                              r""),   # your MODFLOW / PEST++ root
    ("6 CONDA_ENV",  r"C:\miniconda3\envs\flopy",                                                  r""),   # your conda env prefix
]

DRY_RUN     = True      # True = preview only (no files written); False = apply the changes
MAKE_BACKUP = True      # keep a <file>.bak next to each modified file
ALSO_SLASH  = True      # also replace the forward-slash form of each Windows path (some appear as E:/...)
TEXT_EXT    = {".py", ".md", ".txt", ".json", ".yml", ".yaml", ".cfg", ".ini",
               ".toml", ".bat", ".sh", ".pst", ".nam"}
SKIP_DIRS   = {".git", "__pycache__", ".ipynb_checkpoints", "node_modules"}


# ------------------------------------------------------------------ build the effective pair list
def _pairs():
    """Return [(label, old, new), …] including forward-slash variants; new='' means SCAN only."""
    out = []
    for label, old, new in REPLACEMENTS:
        out.append((label, old, new))
        if ALSO_SLASH and "\\" in old:
            out.append((label + " (/)", old.replace("\\", "/"), new.replace("\\", "/") if new else ""))
    # apply the most specific (longest) old string first
    out.sort(key=lambda p: len(p[1]), reverse=True)
    return out


def main():
    pairs = _pairs()
    self_name = Path(__file__).name
    print(f">> scanning {TARGET_DIR}")
    print(f">> mode: {'DRY-RUN (no writes)' if DRY_RUN else 'APPLY'};  "
          f"{sum(1 for _, _, n in REPLACEMENTS if n.strip())}/{len(REPLACEMENTS)} prefixes have a NEW value\n")

    files = [p for p in TARGET_DIR.rglob("*")
             if p.is_file() and p.suffix.lower() in TEXT_EXT
             and not (set(p.parts) & SKIP_DIRS) and p.name != self_name]

    totals = {label: 0 for label, _, _ in REPLACEMENTS}
    changed_files = 0
    for f in sorted(files):
        try:
            text = f.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue                                    # skip binaries / unreadable
        new_text = text
        hits = []                                       # (label, count, replaced?)
        for label, old, new in pairs:
            n = new_text.count(old)
            if not n:
                continue
            base = label.split(" (")[0]
            totals[base] = totals.get(base, 0) + n
            if new:                                      # replace
                new_text = new_text.replace(old, new)
                hits.append((label, n, True))
            else:                                        # scan only
                hits.append((label, n, False))
        if hits:
            rel = f.relative_to(TARGET_DIR)
            tag = "SCAN " if all(not r for _, _, r in hits) else ("would edit" if DRY_RUN else "EDIT")
            print(f"  [{tag}] {rel}")
            for label, n, replaced in hits:
                print(f"         {label:16s} x{n}  {'-> replaced' if replaced else '(no NEW value set -> scan only)'}")
            if new_text != text and not DRY_RUN:
                if MAKE_BACKUP:
                    f.with_suffix(f.suffix + ".bak").write_text(text, encoding="utf-8")
                f.write_text(new_text, encoding="utf-8")
                changed_files += 1

    print("\n== summary (occurrences found) ==")
    for label, _, new in REPLACEMENTS:
        print(f"  {label:16s} {totals.get(label, 0):4d}   {'(NEW set)' if new.strip() else '(NEW empty -> scan only)'}")
    print(f"\n>> files scanned: {len(files)}", end="")
    if DRY_RUN:
        print("   (DRY-RUN: nothing written; set DRY_RUN = False to apply)")
    else:
        print(f"   files modified: {changed_files}"
              f"{'  (.bak backups kept)' if MAKE_BACKUP and changed_files else ''}")


if __name__ == "__main__":
    main()
