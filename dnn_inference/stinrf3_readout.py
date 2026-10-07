#!/usr/bin/env python3
"""
stinrf3_readout.py

Numerical port of the stINRF3 canonical readout so that the DNN pipelines
report direction the same way stINRF3 does:

    1. Temporal pooling with Tb  (beta-shaped window, a=9, b=1)
    2. Spatial averaging with the Gaussian filter c_s (sigma_c, truncated
       to a square support)

Reference MATLAB:
    stINRF3_compute.m        -> makeTb()
    stINRF3_readout_gaussian.m
    run_gratings_2D_par.m    -> do_readout(), case 'gauss'

SIGN CONVENTION
    stINRF3 works in maths convention: 0 deg = rightward, +90 deg = upward.
    Optical-flow networks return v positive DOWNWARD (image convention).
    Feed the readout with  O_H = flow_u,  O_V = -flow_v  and the returned
    angle is directly comparable with stINRF3's.
    Helper flow_to_channels() does this for you.
"""
import math

import numpy as np

# Canonical values from stINRF3_config.m / run_gratings_2D_par.m
TB_A = 9.0
TB_B = 1.0
LEGACY_HALF_REGION = 100                            # px, 200x200 support
LEGACY_SIGMA_C = 2 * LEGACY_HALF_REGION / math.sqrt(12)   # = 57.735 px


# ----------------------------------------------------------------------
# Temporal pooling
# ----------------------------------------------------------------------
def make_tb(n, a=TB_A, b=TB_B):
    """Beta-shaped temporal pooling weights, matching makeTb() in MATLAB.

    MATLAB:
        t  = 0 : n_frames-1;
        Tb = (t/n_frames).^lpa .* (1 - t/n_frames).^lpb;
        Tb = Tb / (sum(Tb) + eps);

    NOTE: run_multiframe_pipeline.py currently uses a midpoint grid
    (t + 0.5/n), which is NOT what stINRF3 does. This function reproduces
    the MATLAB grid exactly.

    Parameters:
        n (int) number of samples (frames, pairs or windows)
        a, b (float) beta shape parameters

    Returns:
        tb (ndarray, n) weights summing to 1
    """
    if n < 1:
        raise ValueError("make_tb: n must be >= 1")
    t = np.arange(n, dtype=np.float64) / float(n)
    tb = (t ** a) * ((1.0 - t) ** b)
    s = tb.sum()
    return tb / s if s > 0 else np.full(n, 1.0 / n)


def make_tb_estimates(n_frames, n_est, grid="midpoints", a=TB_A, b=TB_B):
    """Tb weights for n_est flow estimates derived from an n_frames clip.

    stINRF3 has one response per frame, so Tb has n_frames samples. A flow
    network gives fewer estimates than that: n_frames-1 for a 2-frame model,
    n_frames-10 for an 11-frame window. Two ways to carry Tb across:

      'midpoints' (default) evaluate the Tb function on the ORIGINAL
          n_frames time base, at the instant each estimate refers to (the
          centre of the frames it consumed). This applies the same temporal
          weighting profile to the same physical clip, so it is the faithful
          mapping.
      'pairs'     apply makeTb() directly to the sequence of n_est estimates,
          as if they were the frames.

    The two agree closely (centre of mass within 0.006 of the clip for
    n_frames >= 31) but differ in the weight on the last estimate for short
    clips (0.42 vs 0.45 at n_frames = 12).

    Parameters:
        n_frames (int) frames in the clip
        n_est (int) number of flow estimates
        grid (str) 'midpoints' | 'pairs'

    Returns:
        tb (ndarray, n_est) weights summing to 1
    """
    if n_est < 1:
        raise ValueError("make_tb_estimates: n_est must be >= 1")
    if grid == "pairs":
        return make_tb(n_est, a, b)
    if grid != "midpoints":
        raise ValueError(f"unknown grid '{grid}'")
    span = n_frames - n_est            # frames consumed minus one
    t = (np.arange(n_est, dtype=np.float64) + span / 2.0) / float(n_frames)
    tb = (t ** a) * ((1.0 - t) ** b)
    s = tb.sum()
    return tb / s if s > 0 else np.full(n_est, 1.0 / n_est)


def renormalise(tb, keep):
    """Zero the weights of samples that were not processed and renormalise."""
    keep = np.asarray(keep, dtype=bool)
    w = np.asarray(tb, dtype=np.float64) * keep
    s = w.sum()
    if s > 0:
        return w / s
    if keep.any():
        return keep.astype(np.float64) / keep.sum()
    return np.zeros_like(w)


def make_tb_masked(n, keep, a=TB_A, b=TB_B):
    """Tb weights restricted to the samples actually processed.

    Samples that were skipped (e.g. frame pairs below the contrast
    threshold) get weight 0 and the remaining weights are renormalised.
    Without this the skipped samples act as zeros and bias the pooled
    vector towards the origin.

    Parameters:
        n (int) total number of samples
        keep (ndarray of bool, n) True where the sample was processed

    Returns:
        tb (ndarray, n) weights summing to 1 over the kept samples
    """
    keep = np.asarray(keep, dtype=bool)
    tb = make_tb(n, a, b) * keep
    s = tb.sum()
    if s > 0:
        return tb / s
    # Tb is ~0 wherever we kept samples (can happen if only early samples
    # survive): fall back to a uniform average over them.
    if keep.any():
        return keep.astype(np.float64) / keep.sum()
    return np.zeros(n)


# ----------------------------------------------------------------------
# Spatial readout
# ----------------------------------------------------------------------
def gaussian_weights(shape, sigma_c, support_half=LEGACY_HALF_REGION,
                     row_half=None, centre=None):
    """Gaussian readout weight map W, matching stINRF3_readout_gaussian.m.

    MATLAB uses 1-based indexing with the centre at round(cols/2)+0.5.
    This function reproduces the same sub-pixel centre in 0-based indexing,
    so the weight map is identical element by element.

    Parameters:
        shape (tuple) (rows, cols)
        sigma_c (float) Gaussian sigma [px]
        support_half (float or None) half-width of the square truncation
            along x [px]; None = no truncation
        row_half (float or None) half-width along y [px]; None = support_half
        centre (tuple or None) (cx, cy) in 0-based pixel units; None = the
            stINRF3 convention

    Returns:
        W (ndarray, rows x cols) weights summing to 1
    """
    rows, cols = shape
    if centre is None:
        # MATLAB: cx = round(cols/2) + 0.5  (1-based)  ->  0-based: -1
        # MATLAB round() is half-away-from-zero; floor(x+0.5) matches it
        # for the positive values involved here (np.round is half-to-even).
        cx = math.floor(cols / 2.0 + 0.5) - 0.5
        cy = math.floor(rows / 2.0 + 0.5) - 0.5
    else:
        cx, cy = centre
    if row_half is None:
        row_half = support_half

    yy, xx = np.mgrid[0:rows, 0:cols]
    dx = xx - cx
    dy = yy - cy
    r2 = dx * dx + dy * dy

    W = np.exp(-r2 / (2.0 * sigma_c ** 2)) / (2.0 * np.pi * sigma_c ** 2)
    if support_half is not None and np.isfinite(support_half):
        W[(np.abs(dx) > support_half) | (np.abs(dy) > row_half)] = 0.0

    s = W.sum()
    if s <= 0:
        raise ValueError("gaussian_weights: empty support (sigma_c or "
                         "support_half too small for this canvas)")
    return W / s


def gaussian_readout(O_H, O_V, sigma_c, support_half=LEGACY_HALF_REGION,
                     row_half=None, centre=None, W=None):
    """Spatial readout with the Gaussian filter c_s.

    Port of stINRF3_readout_gaussian.m (scalar sigma_c branch).

    Parameters:
        O_H, O_V (ndarray, rows x cols) Tb-integrated channels, maths
            convention (O_V positive upward)
        sigma_c (float) Gaussian sigma [px]
        W (ndarray or None) precomputed weight map, to avoid rebuilding it
            for every stimulus of the same size

    Returns:
        dict with direction_deg, magnitude, mean_OH, mean_OV, n_eff
    """
    if W is None:
        W = gaussian_weights(O_H.shape, sigma_c, support_half, row_half,
                             centre)
    mean_OH = float(np.sum(W * O_H))
    mean_OV = float(np.sum(W * O_V))
    return {
        "direction_deg": float(np.degrees(np.arctan2(mean_OV, mean_OH))),
        "magnitude": float(math.hypot(mean_OH, mean_OV)),
        "mean_OH": mean_OH,
        "mean_OV": mean_OV,
        "n_eff": float(1.0 / np.sum(W ** 2)),
    }


def sigma_c_from_deg(sigma_c_deg, px_per_deg):
    """Convert a readout width in degrees to pixels.

    The legacy sigma_c = 57.735 px was fixed on stINRF3's own canvas. When
    the DNN export uses a different px_per_deg (DH: 180, Opponency: 200),
    keeping 57.735 px would read out a different area of the visual field
    than stINRF3 did. Specify the width in degrees to keep the two
    comparable, e.g.

        sigma_c_deg = LEGACY_SIGMA_C / 100.0     # DH: stINRF3 ran at 100 px/deg
        sigma_c_px  = sigma_c_from_deg(sigma_c_deg, 180)
    """
    return sigma_c_deg * px_per_deg


# ----------------------------------------------------------------------
# Convenience: full temporal + spatial readout on a stack of flow fields
# ----------------------------------------------------------------------
def flow_to_channels(flow_u, flow_v):
    """Image-convention flow -> stINRF3 channels (O_V positive upward)."""
    return flow_u, -flow_v


def integrate_and_read(flow_u_stack, flow_v_stack, keep=None,
                       sigma_c=LEGACY_SIGMA_C,
                       support_half=LEGACY_HALF_REGION, row_half=None,
                       a=TB_A, b=TB_B, W=None):
    """Tb-integrate a stack of flow fields, then read out with the Gaussian.

    This is the order stINRF3 uses: pool over time FIRST, keeping the
    spatial map, then average in space. Both operations are linear, so the
    scalar result is the same either way -- but doing it in this order also
    gives you the Tb-integrated map, which is what you plot.

    Parameters:
        flow_u_stack, flow_v_stack (ndarray, T x rows x cols) per-pair (or
            per-window) flow, image convention (v positive downward)
        keep (ndarray of bool, T or None) which samples were processed
        sigma_c (float) Gaussian sigma [px]

    Returns:
        (result_dict, O_H, O_V) where O_H/O_V are the Tb-integrated maps
        in maths convention (O_V positive upward)
    """
    T = flow_u_stack.shape[0]
    if keep is None:
        keep = np.ones(T, dtype=bool)
    tb = make_tb_masked(T, keep, a, b)

    O_H = np.tensordot(tb, flow_u_stack, axes=(0, 0))
    O_V = -np.tensordot(tb, flow_v_stack, axes=(0, 0))

    res = gaussian_readout(O_H, O_V, sigma_c, support_half, row_half, W=W)
    res["tb_weights"] = tb
    return res, O_H.astype(np.float32), O_V.astype(np.float32)


# ----------------------------------------------------------------------
# Self-test
# ----------------------------------------------------------------------
def _selftest():
    print("Tb (n=10, a=9, b=1), MATLAB grid t = 0:n-1 over n:")
    print("  ", np.array2string(make_tb(10), precision=5))
    print("  peak at index", int(np.argmax(make_tb(10))), "of 9")

    # A uniform field must read back exactly, whatever the weights.
    O_H = np.full((436, 1024), 3.0)
    O_V = np.full((436, 1024), 4.0)
    r = gaussian_readout(O_H, O_V, LEGACY_SIGMA_C)
    print(f"\nuniform field (3,4): dir={r['direction_deg']:.6f} deg "
          f"(expect 53.130102), mag={r['magnitude']:.6f} (expect 5)")
    print(f"  n_eff={r['n_eff']:.1f} px  "
          f"(square 200x200 would be 40000)")

    # Weight map properties
    W = gaussian_weights((436, 1024), LEGACY_SIGMA_C)
    nz = W > 0
    ys, xs = np.nonzero(nz)
    print(f"\nsupport: rows {ys.min()}..{ys.max()} "
          f"({ys.max()-ys.min()+1} px), "
          f"cols {xs.min()}..{xs.max()} ({xs.max()-xs.min()+1} px)")
    print(f"  sum(W)={W.sum():.10f}")
    # Symmetry of the support about the canvas centre
    print(f"  centred: row mid={(ys.min()+ys.max())/2:.1f} "
          f"(canvas {(436-1)/2:.1f}), "
          f"col mid={(xs.min()+xs.max())/2:.1f} (canvas {(1024-1)/2:.1f})")

    # Masked Tb renormalises
    keep = np.zeros(20, dtype=bool)
    keep[5:15] = True
    tbm = make_tb_masked(20, keep)
    print(f"\nmake_tb_masked: sum={tbm.sum():.10f}, "
          f"zeros outside mask={np.all(tbm[~keep] == 0)}")


if __name__ == "__main__":
    _selftest()
