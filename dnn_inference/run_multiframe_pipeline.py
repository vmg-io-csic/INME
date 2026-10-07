#!/usr/bin/env python3
"""
run_multiframe_pipeline.py

Run ME-Attention or Dual model on full video stimuli using the stINRF3
readout pipeline: sliding window → ROI spatial averaging → Tb filtering.

ME-Attention: 11 input frames, images normalised to [0, 1]
Dual:         15 input frames, images in [0, 255]

Usage:
    # ME-Attention
    python run_multiframe_pipeline.py --model me_attention \
        --model_dir me_attention_dual/ \
        --checkpoint me_attention_dual/Model_example.pth.tar \
        --stimuli_dir stimuli/ --output_dir results_pipeline/

    # Dual
    python run_multiframe_pipeline.py --model dual \
        --model_dir dual_model/ \
        --checkpoint me_attention_dual/final_sintel_kitti.pth \
        --stimuli_dir stimuli/ --output_dir results_pipeline/
"""
import argparse
import os
import sys
import glob
import time
import math
import numpy as np
import scipy.io as sio
import torch
import torch.nn as nn
import torch.nn.functional as F


# Default number of input frames per model
DEFAULT_N_FRAMES = {
    'me_attention': 11,
    'dual': 15,
}


# ============================================================
# Tb filter (identical to stINRF3 makeTb)
# ============================================================
def make_tb(n_frames, a=9, b=1):
    """Beta-shaped temporal pooling weights.

    MATLAB grid, from stINRF3_compute.m -> makeTb():
        t  = 0 : n_frames-1;
        Tb = (t/n_frames).^lpa .* (1 - t/n_frames).^lpb;
    (The previous version used a midpoint grid t + 0.5/n.)
    """
    t = np.arange(n_frames, dtype=np.float64) / float(n_frames)
    tb = (t ** a) * ((1 - t) ** b)
    s = np.sum(tb)
    if s > 0:
        tb = tb / s
    return tb


# Gaussian readout settings, filled in from the CLI in main().
_READOUT = {"sigma_c": None, "support_half": None, "_W": {}}


def _gaussian_W(shape, sigma_c, support_half):
    """Weight map of stINRF3_readout_gaussian.m, cached per canvas size.

    MATLAB is 1-based with the centre at round(cols/2)+0.5; the same
    sub-pixel centre in 0-based indexing is floor(cols/2 + 0.5) - 0.5.
    """
    key = (shape, sigma_c, support_half)
    if key in _READOUT["_W"]:
        return _READOUT["_W"][key]
    rows, cols = shape
    cx = math.floor(cols / 2.0 + 0.5) - 0.5
    cy = math.floor(rows / 2.0 + 0.5) - 0.5
    yy, xx = np.mgrid[0:rows, 0:cols]
    dx = xx - cx
    dy = yy - cy
    W = np.exp(-(dx * dx + dy * dy) / (2.0 * sigma_c ** 2)) \
        / (2.0 * np.pi * sigma_c ** 2)
    if support_half is not None:
        W[(np.abs(dx) > support_half) | (np.abs(dy) > support_half)] = 0.0
    s = W.sum()
    if s <= 0:
        raise ValueError("gaussian readout: empty support")
    W = W / s
    _READOUT["_W"][key] = W
    return W




def extract_roi_mean(flow_u, flow_v, roi_half):
    """Spatial readout with the Gaussian filter c_s.

    Port of stINRF3_readout_gaussian.m. Replaces the flat square ROI so the
    DNN direction is computed exactly as stINRF3 computes it. Returns the
    same (mu, mv) in image convention (v positive downward) as before, so
    all call sites are unchanged.

    sigma_c defaults to 2*roi_half/sqrt(12) = 57.735 px for roi_half = 100.
    """
    sigma_c = _READOUT["sigma_c"]
    support_half = _READOUT["support_half"]
    if support_half is None:
        support_half = roi_half
    if sigma_c is None:
        if roi_half is None:
            return float(np.mean(flow_u)), float(np.mean(flow_v))
        sigma_c = 2.0 * roi_half / math.sqrt(12.0)
    W = _gaussian_W(flow_u.shape, sigma_c, support_half)
    return float(np.sum(W * flow_u)), float(np.sum(W * flow_v))


# ============================================================
# Model loading
# ============================================================
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


def load_me_attention(model_dir, checkpoint, device):
    """Load ME-Attention model (FFV1DNN from me_attention_dual/)."""
    sys.path.insert(0, os.path.abspath(model_dir))
    from FFV1MT_MS import FFV1DNN

    model = FFV1DNN(
        num_scales=8,
        num_cells=256,
        upsample_factor=8,
        feature_channels=256,
        scale_factor=16,
        num_layers=6,
    )

    ckpt = torch.load(checkpoint, map_location='cpu')
    # Handle different checkpoint formats
    if isinstance(ckpt, dict) and 'model_state_dict' in ckpt:
        state = ckpt['model_state_dict']
    elif isinstance(ckpt, dict) and 'state_dict' in ckpt:
        state = ckpt['state_dict']
    else:
        state = ckpt

    # Remove 'module.' prefix if wrapped in DataParallel
    if any(k.startswith('module.') for k in state.keys()):
        state = {k.replace('module.', ''): v for k, v in state.items()}

    model.load_state_dict(state, strict=False)
    model = model.to(device).eval()
    print(f"  ME-Attention loaded ({sum(p.numel() for p in model.parameters()):,} params)")
    return model


def load_dual(model_dir, checkpoint, device):
    """Load Dual model (FFV1DNNV2 from dual_model/)."""
    sys.path.insert(0, os.path.abspath(model_dir))
    from model.nmi6.FFV1MT_MS import FFV1DNNV2

    model = FFV1DNNV2(
        num_scales=8,
        upsample_factor=8,
        scale_factor=16,
        num_layers=6,
    )

    model = nn.DataParallel(model)
    model_dict = torch.load(checkpoint, map_location='cpu')
    model.load_state_dict(model_dict, strict=True)
    model = model.module
    model = model.to(device).eval()
    print(f"  Dual loaded ({sum(p.numel() for p in model.parameters()):,} params)")
    return model


# ============================================================
# Frame preprocessing
# ============================================================
def prepare_frames(tensor, t_start, n_frames, device, model_type):
    """Extract n_frames from tensor and format for the model.

    ME-Attention: [1, 3, H, W] in [0, 1]
    Dual:         [1, 3, H, W] in [0, 255], height/width multiples of 8
    """
    H, W, T = tensor.shape
    frames = []
    for t in range(t_start, min(t_start + n_frames, T)):
        frame = tensor[:, :, t].astype(np.float32)

        if model_type == 'me_attention':
            # Normalise to [0, 1]
            fmax = max(frame.max(), 1.0)
            frame_norm = frame / fmax
        else:
            # Dual: keep in [0, 255] range
            frame_norm = np.clip(frame, 0, 255)

        # Grayscale → RGB (3 channels)
        rgb = np.stack([frame_norm, frame_norm, frame_norm], axis=0)  # (3, H, W)
        t_frame = torch.from_numpy(rgb).unsqueeze(0).float().to(device)  # (1, 3, H, W)
        frames.append(t_frame)

    # Dual: resize to multiples of 8
    if model_type == 'dual' and len(frames) > 0:
        _, _, fH, fW = frames[0].shape
        H8 = int(np.round(fH / 8)) * 8
        W8 = int(np.round(fW / 8)) * 8
        if H8 != fH or W8 != fW:
            frames = [F.interpolate(f, size=(H8, W8), mode='bicubic',
                                    align_corners=True).clamp_(0, 255)
                      for f in frames]

    return frames


# ============================================================
# Forward pass
# ============================================================
def run_model(model, frames, model_type, orig_H, orig_W):
    """Run forward pass and extract flow field (u, v) as numpy arrays (H, W)."""
    with torch.no_grad():
        if model_type == 'me_attention':
            # FFV1DNN.forward returns {'flow_seq': [...]}, the raw flow fields.
            # forward_viz was used before and returns result['flow'] as an
            # RGBA render of the flow (800x1000x4 from matplotlib), not a
            # velocity field, which is why every window was rejected with
            # "Unexpected flow shape".
            result = model.forward(frames, layer=6)
            if isinstance(result, dict):
                flow = result['flow_seq'][-1]
            else:
                flow = result

        elif model_type == 'dual':
            # Dual: forward returns dict with 'flow_seq' list
            result = model.forward(frames, layer=6)
            # result['flow_seq'][-1] is the final flow prediction
            flow = result['flow_seq'][-1]

        # Convert to numpy (H, W, 2) or (2, H, W)
        if isinstance(flow, torch.Tensor):
            flow = flow.detach().cpu().squeeze().numpy()

        if isinstance(flow, np.ndarray):
            if flow.ndim == 3:
                if flow.shape[0] == 2:
                    # (2, H, W) → transpose
                    flow_u = flow[0]
                    flow_v = flow[1]
                elif flow.shape[2] == 2:
                    # (H, W, 2)
                    flow_u = flow[:, :, 0]
                    flow_v = flow[:, :, 1]
                else:
                    raise ValueError(f"Unexpected flow shape: {flow.shape}")
            else:
                raise ValueError(f"Unexpected flow ndim: {flow.ndim}, shape: {flow.shape}")
        else:
            raise ValueError(f"Unexpected flow type: {type(flow)}")

        # Resize flow back to original resolution if needed
        fH, fW = flow_u.shape
        if fH != orig_H or fW != orig_W:
            # Scale flow values proportionally
            scale_h = orig_H / fH
            scale_w = orig_W / fW
            flow_u_t = torch.from_numpy(flow_u[None, None]).float()
            flow_v_t = torch.from_numpy(flow_v[None, None]).float()
            flow_u = F.interpolate(flow_u_t, size=(orig_H, orig_W),
                                   mode='bilinear', align_corners=True).squeeze().numpy() * scale_w
            flow_v = F.interpolate(flow_v_t, size=(orig_H, orig_W),
                                   mode='bilinear', align_corners=True).squeeze().numpy() * scale_h

        return flow_u, flow_v


# ============================================================
# Main pipeline
# ============================================================
def process_stimulus(model, model_type, tensor, device, n_input_frames,
                     roi_half, tb_a, tb_b, skip_threshold=0.0):
    H, W, T = tensor.shape
    n_windows = T - n_input_frames + 1

    if n_windows < 1:
        return None

    u_series = np.zeros(n_windows)
    v_series = np.zeros(n_windows)
    processed = np.zeros(n_windows, dtype=bool)

    for t in range(n_windows):
        # Check signal in this window
        f1 = tensor[:, :, t].astype(np.float32)
        f2 = tensor[:, :, t + 1].astype(np.float32)
        if np.max(np.abs(f2 - f1)) < skip_threshold:
            continue

        frames = prepare_frames(tensor, t, n_input_frames, device, model_type)

        try:
            flow_u, flow_v = run_model(model, frames, model_type, H, W)
            mu, mv = extract_roi_mean(flow_u, flow_v, roi_half)
            u_series[t] = mu
            v_series[t] = mv
            processed[t] = True
        except Exception as e:
            print(f"    [window {t}] error: {e}")
            continue

    n_proc = int(np.sum(processed))
    # (The Tb weights are renormalised over the processed windows further
    # down, where tb is actually built; the patch inserted a duplicate of
    # that here which referenced tb before assignment and OH_map/OV_map,
    # which this function does not have.)
    if n_proc == 0:
        return {
            "tb_dir": 0.0, "tb_mag": 0.0,
            "uniform_dir": 0.0, "uniform_mag": 0.0,
            "weighted_dir": 0.0, "weighted_mag": 0.0,
            "n_windows_processed": 0, "n_windows_total": n_windows,
        }

    # --- Tb filtering (weights renormalised over processed windows) ---
    tb = make_tb(n_windows, a=tb_a, b=tb_b) * processed
    _s = tb.sum()
    tb = tb / _s if _s > 0 else processed / max(n_proc, 1)
    tb_u = float(np.sum(tb * u_series))
    tb_v = float(np.sum(tb * v_series))
    tb_dir = float(np.degrees(np.arctan2(-tb_v, tb_u)))
    tb_mag = float(np.sqrt(tb_u**2 + tb_v**2))

    # --- Uniform average ---
    uni_u = float(np.mean(u_series[processed]))
    uni_v = float(np.mean(v_series[processed]))
    uni_dir = float(np.degrees(np.arctan2(-uni_v, uni_u)))
    uni_mag = float(np.sqrt(uni_u**2 + uni_v**2))

    # --- Contrast-weighted ---
    weights = np.zeros(n_windows)
    for t in range(n_windows):
        if processed[t]:
            f1 = tensor[:, :, t].astype(np.float32)
            f2 = tensor[:, :, t + 1].astype(np.float32)
            weights[t] = float(np.mean(np.abs(f2 - f1)))
    w_sum = np.sum(weights)
    if w_sum > 0:
        wn = weights / w_sum
        cw_u = float(np.sum(wn * u_series))
        cw_v = float(np.sum(wn * v_series))
    else:
        cw_u = cw_v = 0.0
    cw_dir = float(np.degrees(np.arctan2(-cw_v, cw_u)))
    cw_mag = float(np.sqrt(cw_u**2 + cw_v**2))

    return {
        "tb_u": tb_u, "tb_v": tb_v,
        "tb_dir": tb_dir, "tb_mag": tb_mag,
        "uniform_u": uni_u, "uniform_v": uni_v,
        "uniform_dir": uni_dir, "uniform_mag": uni_mag,
        "weighted_u": cw_u, "weighted_v": cw_v,
        "weighted_dir": cw_dir, "weighted_mag": cw_mag,
        "n_windows_processed": n_proc,
        "n_windows_total": n_windows,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Multi-frame model inference with stINRF3 pipeline")
    parser.add_argument("--model", required=True,
        choices=["me_attention", "dual"])
    parser.add_argument("--model_dir", required=True,
        help="me_attention_dual/ or dual_model/")
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--stimuli_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--n_input_frames", type=int, default=None,
        help="Frames per window (default: 11 for ME-Attention, 15 for Dual)")
    parser.add_argument("--roi_half", type=int, default=100)
    parser.add_argument("--tb_a", type=int, default=9)
    parser.add_argument("--tb_b", type=int, default=1)
    parser.add_argument("--skip_threshold", type=float, default=0.0,
        help="0 = process every window, as stINRF3 does")
    parser.add_argument("--sigma_c", type=float, default=None,
        help="Gaussian readout sigma [px]; default 2*roi_half/sqrt(12)")
    parser.add_argument("--sigma_c_deg", type=float, default=None,
        help="Gaussian sigma in DEGREES; overrides --sigma_c "
             "(requires --px_per_deg)")
    parser.add_argument("--px_per_deg", type=float, default=None)
    parser.add_argument("--gpu", type=int, default=0)
    args = parser.parse_args()

    n_frames = args.n_input_frames or DEFAULT_N_FRAMES[args.model]

    device = torch.device(
        f"cuda:{args.gpu}" if torch.cuda.is_available() else "cpu")

    print(f"Loading {args.model}...")
    if args.model == 'me_attention':
        model = load_me_attention(args.model_dir, args.checkpoint, device)
    elif args.model == 'dual':
        model = load_dual(args.model_dir, args.checkpoint, device)

    print(f"Input frames: {n_frames}, ROI: {2*args.roi_half}x{2*args.roi_half}")
    print(f"Tb: a={args.tb_a}, b={args.tb_b}\n")

    _READOUT["support_half"] = args.roi_half
    if args.sigma_c_deg is not None:
        if args.px_per_deg is None:
            parser.error("--sigma_c_deg requires --px_per_deg")
        _READOUT["sigma_c"] = args.sigma_c_deg * args.px_per_deg
    elif args.sigma_c is not None:
        _READOUT["sigma_c"] = args.sigma_c
    else:
        _READOUT["sigma_c"] = 2.0 * args.roi_half / math.sqrt(12.0)
    print(f"Readout: Gaussian sigma_c={_READOUT['sigma_c']:.3f} px, "
          f"support half={args.roi_half} px, skip_threshold={args.skip_threshold}")

    mat_files = sorted(glob.glob(
        os.path.join(args.stimuli_dir, "**", "frames.mat"),
        recursive=True))
    if not mat_files:
        print(f"No frames.mat found under {args.stimuli_dir}")
        return
    print(f"Found {len(mat_files)} stimulus files.\n")

    out_root = os.path.join(args.output_dir, args.model)
    os.makedirs(out_root, exist_ok=True)

    n_ok = n_fail = 0
    t_start = time.time()

    for i, mat_path in enumerate(mat_files):
        rel = os.path.relpath(os.path.dirname(mat_path), args.stimuli_dir)
        out_dir = os.path.join(out_root, rel)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "result.mat")

        try:
            data = load_mat(mat_path)
            if "tensor" not in data:
                print(f"  [{i+1}/{len(mat_files)}] SKIP {rel}: no tensor")
                n_fail += 1
                continue

            tensor = data["tensor"]

            t0 = time.time()
            result = process_stimulus(
                model, args.model, tensor, device,
                n_input_frames=n_frames,
                roi_half=args.roi_half,
                tb_a=args.tb_a, tb_b=args.tb_b,
                skip_threshold=args.skip_threshold)
            dt = time.time() - t0

            if result is None:
                print(f"  [{i+1}/{len(mat_files)}] SKIP {rel}: "
                      f"too few frames (need {n_frames})")
                n_fail += 1
                continue

            print(f"  [{i+1}/{len(mat_files)}]  {rel:40s}  "
                  f"Tb={result['tb_dir']:7.1f}  "
                  f"uni={result['uniform_dir']:7.1f}  "
                  f"cw={result['weighted_dir']:7.1f}  "
                  f"({result['n_windows_processed']}/{result['n_windows_total']} win)  "
                  f"[{dt:.1f}s]")

            save_dict = dict(result)
            save_dict["mode"] = "stINRF3_pipeline_gauss"
            save_dict["readout"] = "gaussian"
            save_dict["sigma_c_px"] = float(_READOUT["sigma_c"]) \
                if _READOUT["sigma_c"] else 2.0*args.roi_half/math.sqrt(12.0)
            save_dict["support_half_px"] = float(args.roi_half)
            save_dict["skip_threshold"] = float(args.skip_threshold)
            save_dict["flow_convention"] = "u right+, v down+"
            save_dict["model_name"] = args.model
            save_dict["roi_half"] = args.roi_half
            save_dict["tb_a"] = args.tb_a
            save_dict["tb_b"] = args.tb_b
            save_dict["n_input_frames"] = n_frames
            sio.savemat(out_path, save_dict)
            n_ok += 1

        except Exception as e:
            print(f"  [{i+1}/{len(mat_files)}] FAIL {rel}: {e}")
            import traceback
            traceback.print_exc()
            n_fail += 1

    elapsed = time.time() - t_start
    print(f"\nDone. {n_ok} OK, {n_fail} failed. "
          f"Total time: {elapsed/60:.1f} min.")


if __name__ == "__main__":
    main()
