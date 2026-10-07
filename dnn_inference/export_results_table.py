#!/usr/bin/env python3
"""
export_results_table.py

Collects every result.mat into two tables:

  results_v25_long.csv   one row per (model, stimulus), every scalar field
  results_v25_wide.csv   one row per stimulus, with <model>_dir / _mag / _err
                         columns side by side -- the layout you want in Excel

and, if pandas + openpyxl are available, results_v25.xlsx with one sheet per
experiment family plus a summary sheet of mean |error| per model.

Angles are in the stINRF3 convention: 0 = rightward, +90 = upward.
Errors are wrapped to [-180, 180].

Usage:
    python export_results_table.py [results_root] [out_prefix]
"""
import csv
import glob
import os
import re
import sys

import numpy as np
import scipy.io as sio

root = sys.argv[1] if len(sys.argv) > 1 else "results_v25"
prefix = sys.argv[2] if len(sys.argv) > 2 else "results_v25"
STIM_ROOT = sys.argv[3] if len(sys.argv) > 3 else "stimuli_by_model"

# stimulus name -> experiment family
FAMILIES = [
    ("barTex", r"^barTex_"), ("gratingTex", r"^gratingTex_"),
    ("gratings", r"^gratings_"), ("plaidI", r"^plaidI_"),
    ("plaidII", r"^plaidII_"), ("DH", r"^DH_"),
    ("Contrast", r"^Contrast_"), ("Masking", r"^Masking_"),
    ("Missing", r"^Missing_"), ("Opponency", r"^Opponency_"),
    ("ReversePhi", r"^ReversePhi_"), ("SecondOrder", r"^SecondOrder_"),
    ("TimeReversal", r"^TR_"),
]

# stimulus parameters worth carrying into the table
PARAMS = ["gt_direction_deg", "speed_px_per_frame", "component_px_per_frame",
          "bar_width_px", "bar_width_deg", "duty_cycle", "inver",
          "texture_rms_contrast", "grating_contrast",
          "cycles_per_frame", "px_per_deg", "framerate", "dov",
          "sfm_cpd", "sfm_envelope_cpd", "f_spatial_cpd", "sigma_t_ms",
          "contrast_level", "mask_ratio", "bar_width", "grain_px",
          "ioc_dir_deg", "direction_deg", "rigid_dir_deg", "type",
          "condition", "n_pairs_processed", "n_windows_processed",
          "n_eff_px", "sigma_c_px"]


def load_meta(path):
    """Scalar/string fields of a frames.mat (v7 or v7.3), tensor skipped."""
    try:
        d = sio.loadmat(path)
        return {k: v for k, v in d.items() if not k.startswith("_")}
    except NotImplementedError:
        pass
    import h5py
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            if k.startswith("#") or k == "tensor":
                continue
            v = f[k][()]
            if isinstance(v, np.ndarray):
                if v.dtype.kind == "u" and \
                        f[k].attrs.get("MATLAB_class") == b"char":
                    v = "".join(chr(c) for c in v.ravel())
            out[k] = v
    return out


# stimulus metadata, read once per stimulus name from whichever canvas group
# has it. The result.mat files do not carry the ground truth, so it is joined
# in from the stimuli here.
STIM_META = {}
for _g in sorted(glob.glob(os.path.join(STIM_ROOT, "*", "v25", "*",
                                        "frames.mat"))):
    _name = os.path.basename(os.path.dirname(_g))
    if _name not in STIM_META:
        try:
            STIM_META[_name] = load_meta(_g)
        except Exception as _e:
            print(f"  aviso: no pude leer {_g}: {_e}")


def family(name):
    for fam, pat in FAMILIES:
        if re.match(pat, name):
            return fam
    return "otro"


def val(d, k):
    if k not in d:
        return None
    a = np.asarray(d[k])
    if a.dtype.kind in "US":
        return str(np.ravel(a)[0]) if a.size else ""
    if a.size == 1:
        return float(np.ravel(a)[0])
    return None


def wrap(x):
    return (x + 180.0) % 360.0 - 180.0


files = sorted(glob.glob(os.path.join(root, "*", "*", "result.mat")))
if not files:
    print(f"no result.mat under {root}/")
    sys.exit(1)

long_rows, long_keys = [], set()
cells = {}                      # (stimulus) -> {model: (dir, mag, err)}
meta = {}                       # (stimulus) -> parameter dict
models = set()

for f in files:
    parts = f.split(os.sep)
    model, stim = parts[-3], parts[-2]
    d = sio.loadmat(f)
    models.add(model)

    r = {"model": model, "experiment": family(stim), "stimulus": stim}
    for k, v in d.items():
        if k.startswith("_"):
            continue
        x = val(d, k)
        if x is not None:
            r[k] = x
    long_keys.update(r)
    long_rows.append(r)

    dirv, magv = val(d, "tb_dir"), val(d, "tb_mag")
    gt = val(d, "gt_direction_deg")
    if gt is None:
        gt = val(STIM_META.get(stim, {}), "gt_direction_deg")
    err = wrap(dirv - gt) if (dirv is not None and gt is not None) else None
    cells.setdefault(stim, {})[model] = (dirv, magv, err)
    if stim not in meta:
        m = {"experiment": family(stim), "stimulus": stim}
        src = dict(STIM_META.get(stim, {}))
        src.update({k: v for k, v in d.items() if not k.startswith("_")})
        for p in PARAMS:
            x = val(src, p)
            if x is not None:
                m[p] = x
        meta[stim] = m

models = sorted(models)

# ---------------- long ----------------
cols = ["model", "experiment", "stimulus"] + sorted(
    k for k in long_keys if k not in ("model", "experiment", "stimulus"))
with open(prefix + "_long.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, restval="")
    w.writeheader()
    w.writerows(long_rows)

# ---------------- wide ----------------
pcols = [p for p in PARAMS if any(p in m for m in meta.values())]
wcols = ["experiment", "stimulus"] + pcols
for m in models:
    wcols += [f"{m}_dir", f"{m}_mag", f"{m}_err"]

wide = []
for stim in sorted(meta, key=lambda s: (meta[s]["experiment"], s)):
    row = dict(meta[stim])
    for m in models:
        dirv, magv, err = cells[stim].get(m, (None, None, None))
        row[f"{m}_dir"] = dirv
        row[f"{m}_mag"] = magv
        row[f"{m}_err"] = err
    wide.append(row)

with open(prefix + "_wide.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=wcols, restval="")
    w.writeheader()
    w.writerows(wide)

print(f"{len(long_rows)} resultados, {len(models)} modelos, "
      f"{len(meta)} estimulos")
print(f"  {prefix}_long.csv   ({len(cols)} columnas)")
print(f"  {prefix}_wide.csv   ({len(wcols)} columnas)")

# ---------------- xlsx ----------------
try:
    import pandas as pd
    df = pd.DataFrame(wide, columns=wcols)
    with pd.ExcelWriter(prefix + ".xlsx", engine="openpyxl") as xl:
        summ = []
        for fam in sorted({r["experiment"] for r in wide}):
            sub = df[df["experiment"] == fam].dropna(axis=1, how="all")
            sub.to_excel(xl, sheet_name=fam[:31], index=False)
            row = {"experiment": fam, "n": len(sub)}
            for m in models:
                c = f"{m}_err"
                if c in sub:
                    row[m] = sub[c].abs().mean()
            summ.append(row)
        pd.DataFrame(summ).to_excel(xl, sheet_name="resumen", index=False)
    print(f"  {prefix}.xlsx        (una hoja por experimento + resumen)")
except ImportError as e:
    print(f"  (sin xlsx: falta {e.name}; usa el CSV ancho en Excel)")

# ---------------- summary to stdout ----------------
print()
print("|error| medio frente a la direccion real, en grados")
print("%-14s %5s %s" % ("experimento", "n",
                        "".join("%22s" % m[:20] for m in models)))
for fam in sorted({r["experiment"] for r in wide}):
    sub = [r for r in wide if r["experiment"] == fam]
    line = "%-14s %5d" % (fam, len(sub))
    for m in models:
        es = [abs(r[f"{m}_err"]) for r in sub
              if r.get(f"{m}_err") is not None]
        line += "%22s" % ("%.2f" % np.mean(es) if es else "-")
    print(line)
