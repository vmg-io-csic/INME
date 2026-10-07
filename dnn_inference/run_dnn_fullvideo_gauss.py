#!/usr/bin/env python3
"""
run_dnn_fullvideo_gauss.py

Drop-in replacement for run_dnn_fullvideo.py that reports direction the
same way stINRF3 does, and saves the flow maps needed to draw optical-flow
figures.

Differences from run_dnn_fullvideo.py
-------------------------------------
1. SPATIAL readout: Gaussian filter c_s (sigma_c, truncated to a square
   support) instead of the mean over the whole frame. Matches
   stINRF3_readout_gaussian.m.
2. TEMPORAL pooling: Tb (beta window, a=9, b=1) as the primary estimate.
   Uniform and contrast-weighted are still saved for comparison, now also
   read out through the Gaussian.
3. Saves the Tb-integrated flow map (flow_u_tb, flow_v_tb) so the optical
   flow can be plotted afterwards. Optionally the per-pair maps too.
4. Skipped pairs no longer bias the pooled vector: their Tb weight is set
   to zero and the rest renormalised.
5. Memory stays O(1) in the number of frames: the threshold pass runs
   first (no model needed), then the flow is accumulated incrementally.

Usage
-----
    python run_dnn_fullvideo_gauss.py --model raft --ckpt sintel \
        --stimuli_dir stimuli/ --output_dir results_gauss/ \
        --sigma_c 57.735 --support_half 100

    # keep the readout area fixed in DEGREES instead of pixels (needed
    # when the export uses a px_per_deg different from stINRF3's):
    python run_dnn_fullvideo_gauss.py ... --sigma_c_deg 0.5774 \
        --px_per_deg 180
"""
import argparse
import glob
import os
import time

import numpy as np
import scipy.io as sio
import torch
import ptlflow

from stinrf3_readout import (LEGACY_HALF_REGION, LEGACY_SIGMA_C, TB_A, TB_B,
                             gaussian_readout, gaussian_weights,
                             make_tb_estimates, renormalise)

DEFAULT_CKPTS = {"liteflownet2": "sintel"}


def load_mat(path):
    """Read a .mat file saved as either v7 or v7.3 (HDF5).

    scipy.io.loadmat cannot read v7.3, which is what MATLAB's
    save(..., '-v7.3') produces and what dnn_v25_save.m uses.

    MATLAB is column-major and HDF5 row-major, so a MATLAB array of size
    (rows, cols, frames) comes back from h5py with its axes reversed. Rather
    than trust that convention blindly, the tensor is checked against the
    N_rows / N_cols the file also stores, and the axes are permuted to match.
    If no permutation fits, it raises instead of silently returning a
    transposed stimulus.

    Returns a dict of numpy arrays, same shape convention as loadmat.
    """
    try:
        return sio.loadmat(path)
    except NotImplementedError:
        pass                                  # v7.3 -> HDF5
    import itertools
    import h5py
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            if k.startswith("#"):
                continue
            v = f[k][()]
            if isinstance(v, np.ndarray):
                if v.ndim > 1:
                    v = np.transpose(v, tuple(range(v.ndim))[::-1])
                if v.dtype.kind == "u" and \
                        f[k].attrs.get("MATLAB_class") == b"char":
                    v = "".join(chr(c) for c in v.ravel())
            out[k] = v

    t = out.get("tensor")
    if isinstance(t, np.ndarray) and t.ndim == 3 \
            and "N_rows" in out and "N_cols" in out:
        want = (int(np.ravel(out["N_rows"])[0]),
                int(np.ravel(out["N_cols"])[0]))
        if t.shape[:2] != want:
            for perm in itertools.permutations(range(3)):
                if np.transpose(t, perm).shape[:2] == want:
                    out["tensor"] = np.transpose(t, perm)
                    break
            else:
                raise ValueError(
                    f"{path}: tensor {t.shape} has no permutation matching "
                    f"(N_rows, N_cols) = {want}")
    return out


def pad_to_min(frame, min_size=128):
    """Pad a 2D frame to the minimum size required by the model."""
    h, w = frame.shape
    pad_h = max(0, min_size - h)
    pad_w = max(0, min_size - w)
    if pad_h == 0 and pad_w == 0:
        return frame, (0, 0, h, w)
    top = pad_h // 2
    left = pad_w // 2
    padded = np.pad(frame,
                    ((top, pad_h - top), (left, pad_w - left)),
                    mode="constant",
                    constant_values=frame.mean())
    return padded, (top, left, h, w)


def run_model_on_pair(model, frame1, frame2, device):
    """Run optical flow on a single frame pair. Returns (flow_u, flow_v)."""
    # ptlflow models expect [0, 1]; passing 0-255 breaks the flow scale
    f1, crop = pad_to_min(np.clip(frame1, 0, 255) / 255.0)
    f2, _ = pad_to_min(np.clip(frame2, 0, 255) / 255.0)
    rgb1 = np.stack([f1] * 3, axis=0)
    rgb2 = np.stack([f2] * 3, axis=0)
    images = torch.from_numpy(np.stack([rgb1, rgb2])[None]).float().to(device)
    with torch.no_grad():
        preds = model({"images": images})
    flow = preds["flows"].cpu().numpy()[0, 0]        # (2, H, W)
    top, left, h, w = crop
    return (flow[0, top:top + h, left:left + w],
            flow[1, top:top + h, left:left + w])


def process_full_video(model, tensor, device, W, sigma_c, support_half,
                       skip_threshold=0.0, save_pairs=False,
                       tb_grid="midpoints"):
    """Process every consecutive frame pair and read out as stINRF3 does.

    Parameters:
        tensor (ndarray, H x W x T) stimulus
        W (ndarray, H x W) precomputed Gaussian weight map
        save_pairs (bool) also return the per-pair flow maps (large)

    Returns:
        (result dict, O_H, O_V, pair_maps or None)
        O_H / O_V are the Tb-integrated channels in maths convention.
    """
    if tensor.ndim != 3:
        raise ValueError(f"Expected 3D tensor, got shape {tensor.shape}")
    H, Wd, T = tensor.shape
    n_pairs = T - 1
    if n_pairs < 1:
        raise ValueError(f"Need >= 2 frames, got {T}")

    # ---- pass 1: which pairs carry signal (no model needed) ------------
    keep = np.zeros(n_pairs, dtype=bool)
    contrast = np.zeros(n_pairs)
    for t in range(n_pairs):
        diff = np.abs(tensor[:, :, t + 1].astype(np.float32)
                      - tensor[:, :, t].astype(np.float32))
        contrast[t] = float(np.mean(diff))
        keep[t] = float(np.max(diff)) >= skip_threshold

    n_processed = int(keep.sum())
    if n_processed == 0:
        empty = {k: 0.0 for k in
                 ("tb_u", "tb_v", "tb_dir", "tb_mag",
                  "uniform_u", "uniform_v", "uniform_dir", "uniform_mag",
                  "weighted_u", "weighted_v", "weighted_dir", "weighted_mag",
                  "fullframe_dir", "fullframe_mag")}
        empty.update(n_pairs_processed=0, n_pairs_total=n_pairs)
        return empty, None, None, None

    # ---- weights -------------------------------------------------------
    tb = renormalise(make_tb_estimates(T, n_pairs, tb_grid, TB_A, TB_B),
                     keep)
    uni = keep.astype(np.float64) / n_processed
    cw = contrast * keep
    cw = cw / cw.sum() if cw.sum() > 0 else uni.copy()

    # ---- pass 2: accumulate the three weighted maps incrementally ------
    acc = {name: [np.zeros((H, Wd), np.float64), np.zeros((H, Wd), np.float64)]
           for name in ("tb", "uniform", "weighted")}
    weights = {"tb": tb, "uniform": uni, "weighted": cw}
    pair_maps = ([], []) if save_pairs else None
    ff_u = ff_v = 0.0                                # legacy full-frame mean

    for t in range(n_pairs):
        if not keep[t]:
            continue
        fu, fv = run_model_on_pair(model,
                                   tensor[:, :, t].astype(np.float32),
                                   tensor[:, :, t + 1].astype(np.float32),
                                   device)
        for name, w in weights.items():
            if w[t] != 0.0:
                acc[name][0] += w[t] * fu
                acc[name][1] += w[t] * fv
        ff_u += uni[t] * float(np.mean(fu))
        ff_v += uni[t] * float(np.mean(fv))
        if save_pairs:
            pair_maps[0].append(fu.astype(np.float32))
            pair_maps[1].append(fv.astype(np.float32))

    # ---- Gaussian readout on each pooled map ---------------------------
    result = {}
    for name in ("tb", "uniform", "weighted"):
        mu, mv = acc[name]
        r = gaussian_readout(mu, -mv, sigma_c, support_half, W=W)
        result[f"{name}_u"] = r["mean_OH"]          # image convention, right+
        result[f"{name}_v"] = -r["mean_OV"]         # image convention, down+
        result[f"{name}_dir"] = r["direction_deg"]
        result[f"{name}_mag"] = r["magnitude"]
    result["n_eff_px"] = float(1.0 / np.sum(W ** 2))

    # legacy full-frame mean, for continuity with earlier runs
    result["fullframe_dir"] = float(np.degrees(np.arctan2(-ff_v, ff_u)))
    result["fullframe_mag"] = float(np.hypot(ff_u, ff_v))

    result["n_pairs_processed"] = n_processed
    result["n_pairs_total"] = n_pairs
    result["tb_weights"] = tb
    result["contrast_per_pair"] = contrast
    result["pair_kept"] = keep

    O_H = acc["tb"][0]
    O_V = acc["tb"][1]                               # image convention (down+)
    if save_pairs:
        pair_maps = (np.stack(pair_maps[0]), np.stack(pair_maps[1]))
    return result, O_H, O_V, pair_maps


def main():
    ap = argparse.ArgumentParser(
        description="Optical flow on full-video stimuli, stINRF3 readout")
    ap.add_argument("--model", required=True)
    ap.add_argument("--ckpt", default=None)
    ap.add_argument("--stimuli_dir", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--skip_threshold", type=float, default=0.0,
                    help="Skip pairs with max|f2-f1| < this. DEFAULT 0 = "
                         "process every pair, which is what stINRF3 does "
                         "(it integrates Tb over all frames, including the "
                         "tails of the temporal envelope). Only raise this "
                         "to save compute, and say so in the methods.")
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--sigma_c", type=float, default=LEGACY_SIGMA_C,
                    help=f"Gaussian sigma [px] (default {LEGACY_SIGMA_C:.3f}, "
                         "the stINRF3 legacy value)")
    ap.add_argument("--sigma_c_deg", type=float, default=None,
                    help="Gaussian sigma in DEGREES; overrides --sigma_c. "
                         "Requires --px_per_deg. Use this when the export "
                         "uses a px_per_deg different from stINRF3's, so "
                         "both read out the same area of the visual field.")
    ap.add_argument("--px_per_deg", type=float, default=None)
    ap.add_argument("--support_half", type=float, default=LEGACY_HALF_REGION,
                    help="Half-side of the square truncation [px]")
    ap.add_argument("--save_flow", choices=["none", "tb", "all"], default="tb",
                    help="none: scalars only. tb: the Tb-integrated map "
                         "(~2-4 MB/stimulus). all: every pair (VERY large).")
    ap.add_argument("--flow_dtype", choices=["float32", "float16"],
                    default="float32")
    ap.add_argument("--tb_grid", choices=["midpoints", "pairs"],
                    default="midpoints",
                    help="How Tb is carried from the n_frames time base to "
                         "the n_frames-1 flow estimates. midpoints (default) "
                         "evaluates Tb at the instant each pair refers to.")
    args = ap.parse_args()

    sigma_c = args.sigma_c
    if args.sigma_c_deg is not None:
        if args.px_per_deg is None:
            ap.error("--sigma_c_deg requires --px_per_deg")
        sigma_c = args.sigma_c_deg * args.px_per_deg

    ckpt = args.ckpt or DEFAULT_CKPTS.get(args.model, "things")
    model_tag = f"{args.model}_{ckpt}"
    device = torch.device(f"cuda:{args.gpu}"
                          if torch.cuda.is_available() else "cpu")

    print(f"Loading {args.model} (checkpoint: {ckpt}) on {device}...")
    model = ptlflow.get_model(args.model, ckpt).to(device).eval()
    print(f"Readout: Gaussian sigma_c={sigma_c:.3f} px, "
          f"support half={args.support_half} px, Tb(a={TB_A}, b={TB_B})")
    print(f"Flow maps: {args.save_flow} ({args.flow_dtype})\n")

    mat_files = sorted(glob.glob(
        os.path.join(args.stimuli_dir, "**", "frames.mat"), recursive=True))
    if not mat_files:
        print(f"No frames.mat found under {args.stimuli_dir}")
        return
    print(f"Found {len(mat_files)} stimulus files.\n")

    out_root = os.path.join(args.output_dir, model_tag)
    os.makedirs(out_root, exist_ok=True)

    W_cache = {}
    n_ok = n_fail = 0
    t_start = time.time()

    for i, mat_path in enumerate(mat_files):
        rel = os.path.relpath(os.path.dirname(mat_path), args.stimuli_dir)
        out_dir = os.path.join(out_root, rel)
        os.makedirs(out_dir, exist_ok=True)

        try:
            data = load_mat(mat_path)
            if "tensor" not in data:
                keys = [k for k in data if not k.startswith("_")]
                print(f"  [{i+1}/{len(mat_files)}] [SKIP] {rel}: keys={keys}")
                n_fail += 1
                continue

            tensor = data["tensor"]
            shape = tensor.shape[:2]
            if shape not in W_cache:
                W_cache[shape] = gaussian_weights(shape, sigma_c,
                                                  args.support_half)
            t0 = time.time()
            result, O_H, O_V, pair_maps = process_full_video(
                model, tensor, device, W_cache[shape], sigma_c,
                args.support_half, args.skip_threshold,
                save_pairs=(args.save_flow == "all"),
                tb_grid=args.tb_grid)
            dt = time.time() - t0

            # readout provenance, so a figure can always be reproduced
            result.update(sigma_c_px=sigma_c, support_half_px=args.support_half,
                          tb_a=TB_A, tb_b=TB_B, readout="gaussian",
                          tb_grid=args.tb_grid,
                          skip_threshold=args.skip_threshold,
                          flow_convention="u right+, v down+")
            for k in ("px_per_deg", "framerate", "dov", "speed_dps",
                      "sfm_cpd", "experiment"):
                if k in data:
                    result[k] = data[k]

            if args.save_flow != "none" and O_H is not None:
                dt_np = np.float16 if args.flow_dtype == "float16" else np.float32
                result["flow_u_tb"] = O_H.astype(dt_np)
                result["flow_v_tb"] = O_V.astype(dt_np)
                if args.save_flow == "all" and pair_maps is not None:
                    result["flow_u_pairs"] = pair_maps[0].astype(dt_np)
                    result["flow_v_pairs"] = pair_maps[1].astype(dt_np)

            sio.savemat(os.path.join(out_dir, "result.mat"), result,
                        do_compression=True)
            print(f"  [{i+1}/{len(mat_files)}] {rel:44s} "
                  f"tb={result['tb_dir']:7.1f}deg  "
                  f"uni={result['uniform_dir']:7.1f}deg  "
                  f"({result['n_pairs_processed']}/{result['n_pairs_total']} "
                  f"pairs) [{dt:.1f}s]")
            n_ok += 1

        except Exception as e:
            print(f"  [{i+1}/{len(mat_files)}] [FAIL] {rel}: {e}")
            n_fail += 1

    print(f"\nDone. {n_ok} OK, {n_fail} failed. "
          f"Total time: {(time.time()-t_start)/60:.1f} min.")


if __name__ == "__main__":
    main()
