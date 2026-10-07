#!/usr/bin/env python3
"""Draw every figure in the paper from the export produced by export_all.m.

MATLAB computes and regenerates the stimuli; this script only draws, so the
two cannot disagree about any number. One script means one font and one
palette across all the figures, which is why the drawing lives here.

The figures are monochrome. Stimuli, panels, frames and bars are black, white
and grey, and the one colour left is a blue that means a single thing: the
direction the stimulus as a whole moves in, that is, the drift of a grating
and the pattern direction of a plaid. Anything that is not that is grey, so a
reader never has to ask what a colour stands for.

Each subfolder of the export holds a meta.json naming its kind, plus the CSVs
that kind needs. Supported kinds:

    panels       stimulus panels plus response bars, either as a filmstrip of
                 successive frames or as a space-time plot
    st_row       space-time panels side by side, one per swept value
    frames_pair  a forward and a reversed stimulus, drawn as a pair of
                 space-time plots (or as filmstrips where the meta asks)
    directions   predicted against stimulus direction, with stimulus icons
    curves       one curve per condition against a swept parameter
    masking      the masking function
    contrast     the contrast response function, with stimulus icons

Usage
-----
    python plot_figures.py                       every figure
    python plot_figures.py ReversePhi            just that one
    python plot_figures.py --export DIR          a different export folder
    python plot_figures.py --out DIR             write the PNGs somewhere else

With no folder given it reads the STINRF_EXPORT environment variable, and
failing that an export folder beside this script, so a clone of the
repository draws its own figures with no setting up.

Naming figures redraws only those, so restyling one after rerunning its
experiment does not mean rebuilding the other eleven.

Writes Fig_<Name>.png at 600 dpi into the export folder.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Polygon

# ------------------------------------------------------------------ palette
INK   = "#000000"     # panel frames, axes, bar outlines, curves
TEXT  = "#333333"     # annotation text
REF   = "#737373"     # reference lines and the unity diagonal
BFILL = "#BFBFBF"     # bar interiors
PATTN = "#0073B3"     # the direction the whole stimulus moves in
COMP1 = "#009E73"     # plaid component 1
COMP2 = "#E69F00"     # plaid component 2

DEG = "°"

# Degrees of visual angle and degrees of orientation are the same unit, so
# both take the degree sign; "deg" and the sign are two spellings of one
# thing and having both in one set of figures reads as an oversight.
SPACE_LAB = f"Space ({DEG})"
TIME_LAB  = "Time (ms)"

N_FRAMES   = 5
FRAME_LW   = 1.0      # panels are framed in a hairline, not in a colour
ARROW_FRAC = 0.75

# ------------------------------------------------------- panel geometry
# Inches. A journal gives a figure 7.2 in across for a double column, and
# the earlier layouts ran from 7.4 to 10.3 -- they were being scaled down to
# fit, which is what made the lettering small and the panels muddy. Every
# figure here is built to sit inside that width at the size it is drawn.
FIG_MAX_W = float(os.environ.get('FIGW', '7.2'))
PANEL_W, PANEL_H = 1.32, 1.76
PANEL_GAP, ROW_GAP = 0.384, 0.544
LAB_H     = 0.544      # condition label above each panel
YLAB_W    = 1.056      # "Time (ms)" and its tick labels, leftmost panel only
XTICK_H   = 0.336
XLAB_H    = 0.432
SB_GAP    = 1.216      # panels to bars, clearing the bars' y label
MAR_T, MAR_B, MAR_R = 0.160, 0.192, 0.448

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Econ Sans Condensed", "DejaVu Sans Condensed", "Liberation Sans Narrow",
    "Nimbus Sans Narrow", "DejaVu Sans",
]
plt.rcParams["axes.linewidth"] = 1.0


def load(folder: Path):
    meta = json.loads((folder / "meta.json").read_text())
    return meta


def as_list(v):
    """MATLAB collapses a one-element cell or array to a scalar in JSON."""
    if v is None:
        return []
    return v if isinstance(v, list) else [v]


def degrees(t: str) -> str:
    """Orientations and directions are in degrees of angle, so they take the
    degree sign. Distances on the retina keep 'deg' for degrees of visual
    angle, which is a different quantity wearing the same word."""
    return t.replace("(deg)", f"({DEG})")


def save(fig, folder: Path, name: str, out: Path = None):
    """PNG only, at 600 dpi. The figures are drawn at the size they are meant
    to be printed, up to a double column, so 600 dpi puts every one of them
    above 3000 pixels across -- past what a journal asks for."""
    stem = (out or folder.parent) / f"Fig_{name}"
    fig.savefig(f"{stem}.png", dpi=600, facecolor="w")
    plt.close(fig)
    w, h = fig.get_size_inches()
    flag = "" if w <= FIG_MAX_W + 1e-6 else "   WIDER THAN A DOUBLE COLUMN"
    print(f"wrote {stem}.png  ({w:.2f} x {h:.2f} in, "
          f"{round(w*600)} x {round(h*600)} px){flag}")


# ------------------------------------------------------------------- panels
def frame_panel(ax, colour=INK, lw=FRAME_LW):
    for side in ax.spines.values():
        side.set_color(colour)
        side.set_linewidth(lw)


def spacetime(ax, stim, x, t_ms, clim):
    """One space-time panel: space across, time down, in grey levels."""
    ax.imshow(stim, cmap="gray", vmin=clim[0], vmax=clim[1], aspect="auto",
              interpolation="nearest", extent=[x[0], x[-1], t_ms[-1], t_ms[0]])
    frame_panel(ax)
    ax.tick_params(colors=INK, direction="out", length=2.5, labelsize=14.40)


def space_axis(meta, ncols):
    """Space in degrees, and the columns the meta asks to keep.

    Falls back to the column count when the geometry is missing, which is what
    the older exports of the time-reversal experiments look like."""
    N = meta.get("N") or ncols
    dov = meta.get("dov")
    if dov is None:
        return np.arange(ncols), np.ones(ncols, bool), False
    x = np.arange(N) * dov / N - dov / 2
    sp = meta.get("space_win")
    keep = np.abs(x) <= sp if sp is not None else np.ones_like(x, bool)
    return x, keep, True


# ------------------------------------------------------------------- bars
def draw_bars(fig, rect, mu, sem, labels, ylabel="Mean response (a.u.)",
              title="Model responses", ylim=None, yticks=None,
              yticklabels=None):
    """Response panel shared by every stimulus layout.

    The bars carry no value printed above them: the axis already says what
    they are, and the printed numbers were the first thing every reader read
    instead of the shape. A zero SEM means the experiment randomised nothing
    -- one fixed phase, one stimulus, one answer -- so no whisker is drawn
    rather than one of length zero."""
    ax = fig.add_axes(rect)
    ax.bar(range(len(mu)), mu, width=0.55, color=BFILL, edgecolor=INK,
           linewidth=1.0)
    for i, s in enumerate(sem):
        if s > 0:
            ax.errorbar(i, mu[i], yerr=s, fmt="none", ecolor=INK,
                        elinewidth=1.1, capsize=4, capthick=1.1)
    ax.axhline(0, color=INK, lw=0.9)
    ax.set_xticks(range(len(mu)))
    ax.set_xticklabels([wrap(l) for l in labels], fontsize=13.60, rotation=(35 if len(labels) > 2 else 0), ha=("right" if len(labels) > 2 else "center"), rotation_mode="anchor")
    ax.set_ylabel(ylabel, fontsize=14.40)
    ax.set_title(title, fontsize=15.20)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", labelsize=13.60)

    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo - 0.10 * (hi - lo), hi + 0.10 * (hi - lo))
    if yticks is not None:
        ax.set_yticks(yticks)
        if yticklabels is not None:
            ax.set_yticklabels(yticklabels)
    return ax


def wrap(t: str) -> str:
    """Break a long tick label at its middle space so neighbours never touch."""
    if len(t) <= 12 or " " not in t:
        return t
    ks = [i for i, ch in enumerate(t) if ch == " "]
    j = min(ks, key=lambda k: abs(k - len(t) / 2))
    return t[:j] + "\n" + t[j + 1:]


# -------------------------------------------------------------- space-time
def bar_width_for(n, labels):
    """Room for n bars and their tick labels.

    Driven by the longest label as well as the count: three bars reading
    "10 ms" need far less width than four reading "Rightward, light bar",
    and sizing on the count alone gave the short-labelled figures an inch of
    white space they then had to be scaled down to fit."""
    longest = max((len(wrap(l).split("\n")[0]) for l in labels), default=6)
    return min(3.0, max(1.55, n * max(0.42, 0.055 * longest) + 0.95) * (1.2 if n <= 2 else 1.6))


def panel_row(fig, wfig, hfig, stims, x, t_ms, clim, labels, *,
              ncols, x0, ytop, pw=PANEL_W, ph=PANEL_H):
    """Lay stimulus panels out in a grid and return the strip's width.

    Only the leftmost panel of each row carries the time axis and its tick
    labels, and the whole strip carries one space label instead of one per
    panel. Repeating them cost about an inch and a half of width and half an
    inch of height per figure, and said the same thing every time."""
    nP = len(stims)
    nrows = int(np.ceil(nP / ncols))
    strip_w = ncols * pw + (ncols - 1) * PANEL_GAP
    for i in range(nP):
        col, row = i % ncols, i // ncols
        px = x0 + col * (pw + PANEL_GAP)
        py = ytop - (row + 1) * (LAB_H + ph) - row * (XTICK_H + ROW_GAP)
        ax = fig.add_axes([px / wfig, py / hfig, pw / wfig, ph / hfig])
        spacetime(ax, stims[i], x, t_ms, clim)
        ax.tick_params(labelsize=12.00)
        ax.xaxis.set_major_locator(plt.MaxNLocator(3))
        if col == 0:
            ax.set_ylabel(TIME_LAB, color=INK, fontsize=14.40, labelpad=2)
        else:
            ax.set_yticklabels([])
        last_row = row == nrows - 1
        if not last_row:
            ax.set_xticklabels([])
        fig.text((px + pw / 2) / wfig, (py + ph + 0.05) / hfig,
                 "\n".join(as_list(labels[i])), color=INK, fontsize=14.40,
                 fontweight="bold", ha="center", va="bottom", linespacing=1.15)
    return strip_w


def panels_figure(stims, x, t_ms, clim, labels, bars, *, title=None,
                  ncols=None, header=None, captions=None):
    """One row (or, for four conditions, two) of space-time panels with the
    response bars beside them. Every panelled figure in the paper is this
    shape, so they are all built here and come out the same size."""
    nP = len(stims)
    if ncols is None:
        ncols = nP if nP <= 3 else 2
    nrows = int(np.ceil(nP / ncols))

    bar_w = bar_width_for(len(bars["mu"]), bars["labels"])
    # The panels take what is left of a double column once the bars and the
    # axis furniture have had theirs, so the figure is never drawn wider than
    # it can be printed and then shrunk, which is what was making the
    # lettering small. They never grow past PANEL_W, only shrink.
    avail = FIG_MAX_W - YLAB_W - SB_GAP - bar_w - MAR_R
    pw = min(PANEL_W, (avail - (ncols - 1) * PANEL_GAP) / ncols)
    ph = min(PANEL_H, pw * 1.9)
    strip_w = ncols * pw + (ncols - 1) * PANEL_GAP
    x0 = YLAB_W
    wfig = x0 + strip_w + SB_GAP + bar_w + MAR_R

    cap_h = 0.384 if captions else 0.0
    hdr_h = 0.448 if header else 0.0
    ttl_h = 0.480 if title else 0.0
    block = nrows * (LAB_H + ph + XTICK_H) + (nrows - 1) * ROW_GAP
    hfig = MAR_T + ttl_h + hdr_h + block + XLAB_H + cap_h + MAR_B

    fig = plt.figure(figsize=(wfig, hfig), facecolor="w")
    if title:
        fig.text(0.5, 1 - (MAR_T + 0.10) / hfig, title, fontsize=17.60,
                 fontweight="bold", ha="center", va="top")
    ytop = hfig - MAR_T - ttl_h - hdr_h
    if header:
        fig.text((x0 + strip_w / 2) / wfig, (ytop + 0.06) / hfig, header,
                 color=INK, fontsize=15.20, fontweight="bold",
                 ha="center", va="bottom")
    panel_row(fig, wfig, hfig, stims, x, t_ms, clim, labels,
              ncols=ncols, x0=x0, ytop=ytop, pw=pw, ph=ph)
    ybot = ytop - block
    fig.text((x0 + strip_w / 2) / wfig, (ybot - 0.03) / hfig, SPACE_LAB,
             color=INK, fontsize=14.40, ha="center", va="top")

    if captions:
        ybot = ybot - XLAB_H
        for i, c in enumerate(captions):
            cx = x0 + (i % ncols) * (pw + PANEL_GAP) + pw / 2
            fig.text(cx / wfig, (ybot + 0.02) / hfig,
                     c.replace("sigma_t", "\u03c3$_t$"), color=TEXT,
                     fontsize=14.40, ha="center", va="bottom")

    # The bars sit inside the same band as the panels, so their title lines
    # up with the condition labels rather than running off the top edge.
    usable = block - LAB_H
    bar_h = min(usable - 0.02, 2.20 if nrows == 1 else 3.60)
    bar_b = ytop - block + (usable - bar_h) / 2
    draw_bars(fig, [(x0 + strip_w + SB_GAP) / wfig, bar_b / hfig,
                    bar_w / wfig, bar_h / hfig],
              bars["mu"], bars["sem"], bars["labels"], **bars.get("kw", {}))
    return fig


def fig_spacetime(folder: Path, meta):
    stims = [np.loadtxt(folder / f"stim_{i+1}.csv", delimiter=",")
             for i in range(len(as_list(meta["dirs"])))]
    x, keep, real = space_axis(meta, stims[0].shape[1])
    t_ms = np.arange(stims[0].shape[0]) / meta["framerate"] * 1000
    clim = (min(s.min() for s in stims), max(s.max() for s in stims))
    bars = bar_spec(meta)
    bars["labels"] = as_list(meta["labels"])
    return panels_figure([st[:, keep] for st in stims], x[keep], t_ms, clim,
                         meta["rowlab"], bars)


def fig_st_row(folder: Path, meta):
    """Panels side by side, one per swept value, where the point is how the
    picture changes across the sweep. The stimulus drifts the same way in
    every one, so the direction is stated once above the row."""
    caps = as_list(meta["captions"])
    stims = [np.loadtxt(folder / f"stim_{i+1}.csv", delimiter=",")
             for i in range(len(caps))]
    x, keep, real = space_axis(meta, stims[0].shape[1])
    t_ms = np.arange(stims[0].shape[0]) / meta["framerate"] * 1000
    clim = (min(s.min() for s in stims), max(s.max() for s in stims))
    word = str(meta["direction_word"]).lower()
    bars = bar_spec(meta)
    bars["labels"] = as_list(meta["labels"])
    return panels_figure([st[:, keep] for st in stims], x[keep], t_ms, clim,
                         [[c.replace("sigma_t", "\u03c3$_t$")] for c in caps],
                         bars, header=f"Motion direction: {word}")


def fig_spacetime_pair(folder: Path, meta):
    """A forward stimulus beside its time reverse.

    A filmstrip of five frames shows five instants and leaves the reader to
    infer what happened between them, and for a pattern that drifts by a
    fraction of a bar per frame there is nothing to see in any one of them.
    The space-time plot shows the whole trial at once, and the time reversal
    is then the same picture flipped top to bottom, which is the claim the
    figure is making."""
    fwd = np.loadtxt(folder / "stim_1.csv", delimiter=",")
    rev = np.loadtxt(folder / "stim_2.csv", delimiter=",")
    mism = np.abs(rev - fwd[::-1]).max()
    if mism > 1e-9 * max(1.0, np.abs(fwd).max()):
        print(f"  warning: stim_rev is not the exact time reverse of stim_fwd "
              f"(max difference {mism:.3g}); the two panels may not mirror")

    x, keep, real = space_axis(meta, fwd.shape[1])
    fr = meta.get("framerate")
    t_ms = (np.arange(fwd.shape[0]) / fr * 1000 if fr
            else np.arange(fwd.shape[0]))
    clim = (min(fwd.min(), rev.min()), max(fwd.max(), rev.max()))
    labels = as_list(meta["labels"])
    bars = {"mu": np.array(as_list(meta["mean"]), float),
            "sem": np.array(as_list(meta["sem"]), float),
            "labels": labels}
    return panels_figure([fwd[:, keep], rev[:, keep]], x[keep], t_ms, clim,
                         [[l] for l in labels], bars,
                         title=meta.get("title", meta["name"]))


# --------------------------------------------------------------- filmstrip
def fig_frames(folder: Path, meta):
    dirs = as_list(meta["dirs"])
    stims = [np.loadtxt(folder / f"stim_{i+1}.csv", delimiter=",")
             for i in range(len(dirs))]
    rowlab = meta["rowlab"]
    nR = len(stims)

    x, keep, real = space_axis(meta, stims[0].shape[1])
    nf = stims[0].shape[0]
    span = min(max(N_FRAMES - 1, round(meta["frame_span"] * nf)), nf - 1)
    idx = np.round(np.linspace(0, span, N_FRAMES)).astype(int)
    clim = (min(s.min() for s in stims), max(s.max() for s in stims))

    # Laid out in inches and converted at the end. Working in figure fractions
    # made the bottom row's arrow fall off the canvas as soon as a fourth row
    # was added, because the margins did not scale with the row count.
    side, fgap = 0.92, 0.10
    strip_w = N_FRAMES * side + (N_FRAMES - 1) * fgap
    mar_l, gap_sb, bar_w, mar_r = 1.05, 0.95, 2.90, 0.35
    wfig = mar_l + strip_w + gap_sb + bar_w + mar_r

    num_h, a_gap, t_lab = 0.26, 0.16, 0.28
    row_block = num_h + side + a_gap + t_lab
    mar_t, mar_b, r_gap = 0.24, 0.22, 0.24
    hfig = mar_t + nR * row_block + (nR - 1) * r_gap + mar_b

    fig = plt.figure(figsize=(wfig, hfig), facecolor="w")
    fx = [mar_l + k * (side + fgap) for k in range(N_FRAMES)]
    ytop = hfig - mar_t
    y_bot = []
    yc = ytop
    for r in range(nR):
        yc -= num_h
        row_y = yc - side
        y_bot.append(row_y)
        for k, fi in enumerate(idx):
            ax = fig.add_axes([fx[k] / wfig, row_y / hfig,
                               side / wfig, side / hfig])
            ax.imshow(np.tile(stims[r][fi, keep], (2, 1)), cmap="gray",
                      vmin=clim[0], vmax=clim[1], aspect="auto",
                      interpolation="nearest")
            ax.set_xticks([]); ax.set_yticks([])
            frame_panel(ax)
            ax.set_title(str(k + 1), color=INK, fontsize=16.80,
                         fontweight="bold", pad=3)
            if k == 0:
                ax.set_ylabel("\n".join(as_list(rowlab[r])), color=INK,
                              fontsize=16.00, fontweight="bold", labelpad=8,
                              linespacing=1.25)
        ya = (row_y - a_gap) / hfig
        aw = ARROW_FRAC * strip_w / wfig
        ax1 = (mar_l + (strip_w - ARROW_FRAC * strip_w) / 2) / wfig
        fig.add_artist(FancyArrow(ax1, ya, aw, 0, transform=fig.transFigure,
                                  width=0.0009, head_width=0.008,
                                  head_length=0.009,
                                  length_includes_head=True, color=TEXT))
        fig.text((mar_l + strip_w / 2) / wfig, ya - 0.012, "time", color=TEXT,
                 fontsize=16.00, style="italic", ha="center", va="top")
        yc = row_y - a_gap - t_lab
        if r < nR - 1:
            yc -= r_gap

    avail = ytop - y_bot[-1]
    bar_h = min(avail, 3.2)
    bar_b = y_bot[-1] + (avail - bar_h) / 2
    draw_bars(fig, [(mar_l + strip_w + gap_sb) / wfig, bar_b / hfig,
                    bar_w / wfig, bar_h / hfig],
              np.array(as_list(meta["mean"])), np.array(as_list(meta["sem"])),
              as_list(meta["labels"]))
    return fig


# ------------------------------------------------------------------- icons
def stim_icon(ax, comp_dirs, cycles, comp_cols=None, pattern_dir=None,
              pattern_col=None, amp=1.0, pad=2.35):
    """Circular stimulus icon in the style of Movshon plaid figures: a light
    grey disc carrying the square-wave components, one arrow per component,
    and a longer arrow for the pattern direction where there is one. The wave
    vector of each component is aligned with its drift, so its bars come out
    perpendicular to its own arrow.

    comp_cols of None draws no arrows at all, for the figures whose x axis is
    not direction and where an arrow on every icon said the same thing five
    times over. pad then comes down to fit the bare disc."""
    n, bg, duty = 420, 0.86, 0.42
    v = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(v, v)
    inside = X ** 2 + Y ** 2 <= 1
    comp_dirs = list(np.atleast_1d(comp_dirs))
    cycles = np.broadcast_to(np.atleast_1d(cycles), (len(comp_dirs),))

    L = np.zeros_like(X)
    for th_deg, cy in zip(comp_dirs, cycles):
        th = np.deg2rad(th_deg)
        u = X * np.cos(th) + Y * np.sin(th)
        L += (np.mod(u * cy, 1) >= duty).astype(float)
    img = 0.5 * bg + amp * bg * (L / len(comp_dirs) - 0.5)
    img = np.where(inside, img, np.nan)

    ax.imshow(img, cmap="gray", vmin=0, vmax=1, extent=[-1, 1, -1, 1],
              origin="lower", interpolation="nearest")
    a = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(a), np.sin(a), "-", color=INK, lw=1.0)
    if comp_cols is not None:
        comp_cols = list(np.atleast_1d(comp_cols))
        if len(comp_cols) < len(comp_dirs):
            comp_cols = comp_cols * len(comp_dirs)
        for th_deg, cc in zip(comp_dirs, comp_cols):
            th = np.deg2rad(th_deg)
            u = np.array([np.cos(th), np.sin(th)])
            arrow(ax, 1.08 * u, 1.62 * u, cc, 1.8, 0.34)
        if pattern_dir is not None:
            th = np.deg2rad(pattern_dir)
            u = np.array([np.cos(th), np.sin(th)])
            arrow(ax, 1.08 * u, 2.05 * u, pattern_col or PATTN, 2.6, 0.42)
    ax.set_xlim(-pad, pad); ax.set_ylim(-pad, pad)
    ax.set_aspect("equal"); ax.axis("off")


def arrow(ax, p0, p1, colour, lw, head):
    d = (p1 - p0) / max(np.linalg.norm(p1 - p0), 1e-12)
    stem = p1 - 0.85 * head * d
    ax.plot([p0[0], stem[0]], [p0[1], stem[1]], "-", color=colour, lw=lw)
    nvec = np.array([-d[1], d[0]])
    base = p1 - head * d
    ax.add_patch(Polygon([p1, base + 0.5 * head * nvec,
                          base - 0.5 * head * nvec],
                         closed=True, facecolor=colour, edgecolor=colour))


# ----------------------------------------------------------------- masking
def fig_masking(folder: Path, meta):
    d = np.loadtxt(folder / "curve.csv", delimiter=",")
    x, y = d[:, 0], d[:, 1]
    sig = meta["signal_cpd"]
    # The signal-alone level comes from the dedicated no-mask run, exported
    # beside this one. A constant copied into the masking script goes stale
    # the moment the model or the parameters change.
    ref = None
    sib = folder.parent / "MaskingNoMask" / "meta.json"
    if sib.exists():
        ref = json.loads(sib.read_text()).get("value")
    else:
        print("  no MaskingNoMask export found, drawing without the "
              "signal-alone line; run run_Masking_nomask_1D")

    fig = plt.figure(figsize=(6.6, 3.8), facecolor="w")
    ax = fig.add_axes([0.125, 0.185, 0.83, 0.72])
    lo = min(list(y) + ([ref] if ref else []))
    hi = max(list(y) + ([ref] if ref else []))
    rng = hi - lo
    ax.set_ylim(lo - 0.10 * rng, hi + 0.20 * rng)

    if ref:
        ax.axhline(ref, ls="--", color=REF, lw=1.1)
        ax.text(x.max() * 1.3, ref, "Signal alone, no mask ", color=REF,
                fontsize=12.80, ha="right", va="bottom")
    ax.axvline(sig, ls=":", color=INK, lw=0.9)
    ax.text(sig, ax.get_ylim()[1], f" Signal {sig:g} cpd", fontsize=12.80,
            va="top")

    # No marker singled out at the minimum and no label on it. Where the
    # function bottoms out is the reading to be taken off the axis, and
    # picking the point out in advance made it a claim about one frequency
    # of the thirteen that were run.
    ax.plot(x, y, "-", color=INK, lw=1.7)
    ax.plot(x, y, "o", color=INK, mfc="w", ms=5, mew=1.3)

    ax.set_xscale("log")
    ax.set_xlim(x.min() * 0.7, x.max() * 1.4)
    # Enough rungs in between that the eye can place a point on the axis. The
    # old four ticks left the whole decade between 0.5 and 2 unmarked, which
    # is where the function actually turns over.
    ticks = [t for t in [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 6, 12]
             if ax.get_xlim()[0] <= t <= ax.get_xlim()[1]]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel("Mask spatial frequency (cpd)", fontsize=16.00)
    ax.set_ylabel("Response (a.u.)", fontsize=16.00)
    ax.set_title("Motion masking", fontsize=17.60)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")
    return fig


# ---------------------------------------------------------------- contrast
def fig_contrast(folder: Path, meta):
    d = np.loadtxt(folder / "curve.csv", delimiter=",")
    x, y = d[:, 0], d[:, 1]

    fig = plt.figure(figsize=(5.6, 4.05), facecolor="w")
    ax = fig.add_axes([0.165, 0.265, 0.785, 0.635])
    ax.plot(x, y, "-", color=INK, lw=1.7)
    ax.plot(x, y, "o", color=INK, mfc="w", ms=4, mew=1.2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Contrast", fontsize=16.00)
    ax.set_ylabel("Response (a.u.)", fontsize=16.00)
    ax.set_title("Contrast response", fontsize=17.60)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")

    want = np.logspace(np.log10(x.min()), np.log10(x.max()), 5)
    idx = sorted({int(np.argmin(abs(x - w))) for w in want})
    ax.set_xticks(x[idx])
    ax.set_xticklabels([f"{v:.2g}" for v in x[idx]])
    ax.minorticks_off()

    p = ax.get_position()
    wfig, hfig = fig.get_size_inches()
    xl = np.log10(ax.get_xlim())
    fx = (np.log10(x[idx]) - xl[0]) / (xl[1] - xl[0])
    # The icons carry no arrow: this axis sweeps contrast, the drift is the
    # same in all five, and five copies of one arrow said nothing. Without
    # them the disc no longer has to share its box, so the box is drawn from
    # the diameter wanted rather than the other way round -- here a tenth
    # larger than the arrowed icons used to come out.
    iw = min(0.098, 0.92 * min(np.diff(fx)) * p.width) if len(fx) > 1 else 0.098
    disc_in = 1.10 * iw * wfig * (2 / 2.35 / 2)
    pad = 1.08
    box_in = disc_in * pad
    icon_top = p.y0 - 0.55 / hfig          # clear of the x axis label
    for f, i in zip(fx, idx):
        a = fig.add_axes([p.x0 + f * p.width - box_in / wfig / 2,
                          icon_top - box_in / hfig,
                          box_in / wfig, box_in / hfig])
        stim_icon(a, [0], [2.2], None, amp=x[i], pad=pad)
    return fig


def fig_directions(folder: Path, meta):
    """Predicted against stimulus direction, with a stimulus icon under each
    condition. Only the fixed-phase condition is exported, so there is no
    dispersion and no error bar: the claim is that the model reports this
    direction, not that it reports it to a given precision.

    Blue is the direction the stimulus as a whole moves in. For a plaid that
    is the pattern direction and the components are drawn in their own two
    colours; for a single grating there is nothing but the global direction,
    so the grating's one arrow is blue as well rather than a third colour
    meaning the same thing."""
    x = np.array(as_list(meta["x"]), float)
    y = np.array(as_list(meta["y"]), float)
    icons = meta["icons"] if isinstance(meta["icons"], list) else [meta["icons"]]
    single = meta["line_role"] == "comp1"
    line_col = PATTN
    icon_cols = [PATTN, COMP2] if single else [COMP1, COMP2]

    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    pad = 0.12 * (hi - lo)
    lim = (lo - pad, hi + pad)

    # Laid out in inches so the canvas ends where the icons do. Where the
    # conditions crowd, the icons alternate between two rows rather than
    # shrinking -- that is what used to make the Type-II icons half the size
    # of the others -- and the figure grows by exactly that second row.
    wfig, pw = 6.0, 4.78
    icon_in = 0.58
    mar_l, ph, mar_t = 0.92, 2.85, 0.24
    xlab_h, igap, mar_b = 0.42, 0.08, 0.14
    fx = (x - lim[0]) / (lim[1] - lim[0])
    gapmin = min(np.diff(np.sort(fx))) * pw if len(x) > 1 else pw
    stagger = gapmin < 1.02 * icon_in
    drop = 0.62 * icon_in if stagger else 0.0
    hfig = mar_b + icon_in + drop + igap + xlab_h + ph + mar_t

    fig = plt.figure(figsize=(wfig, hfig), facecolor="w")
    py = mar_b + icon_in + drop + igap + xlab_h
    ax = fig.add_axes([mar_l / wfig, py / hfig, pw / wfig, ph / hfig])
    ax.plot(lim, lim, "--", color=REF, lw=1.0)
    ax.plot(x, y, "-", color=line_col, lw=1.6)
    ax.plot(x, y, "o", color=line_col, mfc=line_col, ms=5)
    ax.set_xlim(lim); ax.set_ylim(lim)

    # keep only ticks far enough apart to stay legible
    tk, last = [], -np.inf
    for v in np.unique(np.round(x)):
        if v - last >= 0.07 * (lim[1] - lim[0]):
            tk.append(v); last = v
    ax.set_xticks(tk); ax.set_yticks(tk)
    ax.set_xlabel(degrees(meta["xlabel"]), fontsize=16.00)
    ax.set_ylabel(f"Predicted direction ({DEG})", fontsize=16.00)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")

    for i in range(len(x)):
        dy = drop if (stagger and i % 2 == 1) else 0.0
        a = fig.add_axes([(mar_l + fx[i] * pw - icon_in / 2) / wfig,
                          (mar_b + drop - dy) / hfig,
                          icon_in / wfig, icon_in / hfig])
        ic = icons[i]
        stim_icon(a, as_list(ic["comp"]), 2.2, icon_cols,
                  pattern_dir=ic.get("pattern") or None)
    return fig


def bar_spec(meta):
    """What the response panel plots.

    Where the experiment exported the share of trials that came out in the
    direction the stimulus actually moved, that is what the bars show, as a
    percentage on a fixed 0-100 axis: the question the experiment asks is
    whether the model gets the direction right, and a mean response in
    arbitrary units answers it only once the reader has worked out which sign
    means which way. Where that was not exported, the bars fall back to the
    mean response."""
    if meta.get("correct") is None:
        return {"mu": np.array(as_list(meta["mean"]), float),
                "sem": np.array(as_list(meta["sem"]), float),
                "kw": {}}
    pc = 100.0 * np.array(as_list(meta["correct"]), float)
    se = 100.0 * np.array(as_list(meta.get("correct_sem")), float) \
        if meta.get("correct_sem") is not None else np.zeros_like(pc)
    return {"mu": pc, "sem": se,
            "kw": {"ylabel": "Correct answers (%)",
                   "title": "Model responses",
                   "ylim": (0, 105), "yticks": [0, 25, 50, 75, 100]}}


def period_frames(stim):
    """How many frames the pattern takes to return to where it started.

    Correlates every frame with the first and takes the first lag at which the
    match comes back, after having dropped away. Taking the best match instead
    lands on a multiple of the period, which is what put five frames across
    several cycles and made them look unrelated.

    Returns a twelfth of the trial when nothing periodic is found, which is a
    sensible window for a stimulus that does not repeat.
    """
    n = stim.shape[0]
    a = stim - stim.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(a, axis=1)
    norm[norm == 0] = 1.0
    r = (a @ a[0]) / (norm * norm[0])
    dropped = False
    for t in range(2, n):
        if not dropped:
            if r[t] < 0.5:
                dropped = True
        elif r[t] > 0.9:
            return t
    return max(N_FRAMES - 1, int(0.12 * n))


def fig_frames_pair(folder: Path, meta):
    """A forward filmstrip over a reversed one, for the time-reversal stimuli.

    The reversed row is read at mirrored indices. Taking the same time indices
    in both rows makes the middle frames differ by one sample whenever the
    frame count is even, so the pair numbered alike would not be alike."""
    fwd = np.loadtxt(folder / "stim_1.csv", delimiter=",")
    rev = np.loadtxt(folder / "stim_2.csv", delimiter=",")

    n = fwd.shape[0]
    win = meta.get("frame_window")
    if win:
        # Set in the experiment's export block when the automatic estimate
        # does not suit, for example frame_window = [1 121].
        first, last = int(win[0]) - 1, int(win[1]) - 1
        print(f"  frame window {win[0]}-{win[1]} of {n}, set by hand")
    else:
        # These stimuli are periodic in space, so the pattern returns to where
        # it started every period. Frames spread over the whole trial land on
        # nearly the same phase each time and come out almost identical. One
        # period, minus a hair so the last frame does not repeat the first,
        # gives a clean progression instead.
        per = period_frames(fwd)
        first, last = 0, int(round(0.97 * per))
        print(f"  stimulus period {per} frames of {n}; showing frames "
              f"{first + 1}-{last + 1}")
    last = min(max(last, N_FRAMES - 1), n - 1)
    idx = np.clip(np.round(np.linspace(first, last, N_FRAMES)).astype(int),
                  0, n - 1)
    idx_rev = (n - 1) - idx[::-1]

    mism = np.abs(rev - fwd[::-1]).max()
    if mism > 1e-9 * max(1.0, np.abs(fwd).max()):
        print(f"  warning: stim_rev is not the exact time reverse of stim_fwd "
              f"(max difference {mism:.3g}); paired frames may not match")

    mu = np.array(as_list(meta["mean"]), float)
    sem = np.array(as_list(meta["sem"]), float)
    clim = (min(fwd.min(), rev.min()), max(fwd.max(), rev.max()))

    side, fgap = 0.92, 0.10
    strip_w = N_FRAMES * side + (N_FRAMES - 1) * fgap
    mar_l, gap_sb, bar_w, mar_r = 1.05, 0.95, 2.60, 0.35
    wfig = mar_l + strip_w + gap_sb + bar_w + mar_r
    num_h, a_gap, t_lab = 0.26, 0.16, 0.28
    row_block = num_h + side + a_gap + t_lab
    mar_t, mar_b, r_gap = 0.55, 0.22, 0.34
    hfig = mar_t + 2 * row_block + r_gap + mar_b

    fig = plt.figure(figsize=(wfig, hfig), facecolor="w")
    fig.suptitle(meta.get("title", meta["name"]), fontsize=22.40,
                 fontweight="bold", y=1 - 0.18 / hfig)

    ytop = hfig - mar_t
    yc = ytop
    y_bot = []
    labels = as_list(meta["labels"])
    rows = [(labels[0], fwd, idx, range(1, N_FRAMES + 1)),
            (labels[1], rev, idx_rev, range(N_FRAMES, 0, -1))]
    for label, stim, ids, nums in rows:
        yc -= num_h
        row_y = yc - side
        y_bot.append(row_y)
        for k, (fi, num) in enumerate(zip(ids, nums)):
            ax = fig.add_axes([(mar_l + k * (side + fgap)) / wfig,
                               row_y / hfig, side / wfig, side / hfig])
            ax.imshow(np.tile(stim[fi, :], (2, 1)), cmap="gray",
                      vmin=clim[0], vmax=clim[1], aspect="auto",
                      interpolation="nearest")
            ax.set_xticks([]); ax.set_yticks([])
            frame_panel(ax)
            ax.set_title(str(num), color=INK, fontsize=16.80,
                         fontweight="bold", pad=3)
            if k == 0:
                ax.set_ylabel(label, color=INK, fontsize=17.60,
                              fontweight="bold", labelpad=8)
        ya = (row_y - a_gap) / hfig
        aw = ARROW_FRAC * strip_w / wfig
        ax1 = (mar_l + (1 - ARROW_FRAC) * strip_w / 2) / wfig
        fig.add_artist(FancyArrow(ax1, ya, aw, 0, transform=fig.transFigure,
                                  width=0.0009, head_width=0.008,
                                  head_length=0.009,
                                  length_includes_head=True, color=TEXT))
        fig.text((mar_l + strip_w / 2) / wfig, ya - 0.012, "time", color=TEXT,
                 fontsize=16.00, style="italic", ha="center", va="top")
        yc = row_y - a_gap - t_lab - r_gap

    avail = ytop - y_bot[-1]
    bar_h = min(avail, 3.2)
    draw_bars(fig, [(mar_l + strip_w + gap_sb) / wfig,
                    (y_bot[-1] + (avail - bar_h) / 2) / hfig,
                    bar_w / wfig, bar_h / hfig], mu, sem, labels)
    return fig


SERIES_MARK = ["o", "s", "D", "^", "v"]
SERIES_FILL = ["w", INK, "w", INK, "w"]


def fig_curves(folder: Path, meta):
    """One curve per condition against the envelope width, with SEM whiskers.
    The zero line is drawn because the result of interest is where a curve
    crosses it: that is the envelope width at which the reported direction
    reverses. The curves are told apart by their marker rather than by their
    colour, so the figure survives being printed in black and white."""
    x = np.array(as_list(meta["x"]), float)
    series = meta["series"] if isinstance(meta["series"], list) else [meta["series"]]

    fig = plt.figure(figsize=(7.2, 5.2), facecolor="w")
    ax = fig.add_axes([0.135, 0.135, 0.815, 0.80])
    for i, sr in enumerate(series):
        mu = np.array(as_list(sr["mean"]), float)
        se = np.array(as_list(sr["sem"]), float)
        ax.errorbar(x, mu, yerr=se, fmt="-" + SERIES_MARK[i % len(SERIES_MARK)],
                    color=INK, ecolor=INK, lw=1.5, ms=6,
                    mfc=SERIES_FILL[i % len(SERIES_FILL)], mew=1.4,
                    capsize=4, capthick=1.2,
                    label=sr.get("label") or sr.get("name"))
    ax.axhline(0, color=INK, ls="--", lw=1.0)
    if meta.get("xlog"):
        ax.set_xscale("log")
        # Label only the widths that are far enough apart on the log axis.
        # Ticking every swept value put 85, 89.2, 95 and 100 on top of one
        # another, which is what made the axis unreadable.
        lx = np.log10(x)
        span = lx.max() - lx.min()
        keep, last = [], -np.inf
        for i, v in enumerate(lx):
            if v - last >= 0.055 * span:
                keep.append(i); last = v
        if keep[-1] != len(x) - 1 and lx[-1] - lx[keep[-1]] >= 0.03 * span:
            keep.append(len(x) - 1)
        ax.set_xticks(x[keep])
        ax.set_xticklabels([f"{v:g}" for v in x[keep]], fontsize=14.40)
        ax.minorticks_off()
    ax.set_xlabel(degrees(meta["xlabel"].replace("sigma_t", "σ$_t$")),
                  fontsize=16.00)
    ax.set_ylabel(degrees(meta["ylabel"]), fontsize=16.00)
    ax.set_title(meta.get("title", ""), fontsize=17.60)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")
    ax.legend(frameon=False, fontsize=14.40)
    return fig


def main(root: str, only=None, out=None):
    root = Path(root)
    dirs = sorted(p for p in root.iterdir() if p.is_dir())
    if only:
        wanted = {w.lower() for w in only}
        missing = wanted - {p.name.lower() for p in dirs}
        for m in sorted(missing):
            print(f"no figure named {m} in {root}")
        dirs = [p for p in dirs if p.name.lower() in wanted]
        if not dirs:
            sys.exit(f"nothing to draw; {root} holds: "
                     + ", ".join(p.name for p in sorted(root.iterdir())
                                 if p.is_dir()))
    for d in dirs:
        meta = load(d)
        kind = meta["kind"]
        if kind == "reference":
            continue        # read while drawing another figure, not drawn
        if kind == "masking":
            fig = fig_masking(d, meta)
        elif kind == "contrast":
            fig = fig_contrast(d, meta)
        elif kind == "directions":
            fig = fig_directions(d, meta)
        elif kind == "curves":
            fig = fig_curves(d, meta)
        elif kind == "st_row":
            fig = fig_st_row(d, meta)
        elif kind == "frames_pair":
            # Space-time unless the experiment asks for a filmstrip.
            fig = (fig_frames_pair(d, meta) if meta.get("panel") == "frames"
                   else fig_spacetime_pair(d, meta))
        elif meta["panel"] == "spacetime":
            fig = fig_spacetime(d, meta)
        else:
            fig = fig_frames(d, meta)
        save(fig, d, meta["name"], out)


def default_export_dir() -> Path:
    """Where the export lives unless told otherwise."""
    env = os.environ.get("STINRF_EXPORT")
    return Path(env) if env else Path(__file__).resolve().parent / "export"


if __name__ == "__main__":
    args = sys.argv[1:]
    root = default_export_dir()
    if "--export" in args:
        i = args.index("--export")
        try:
            root = Path(args[i + 1])
        except IndexError:
            sys.exit("--export needs a folder")
        del args[i:i + 2]
    out = None
    if "--out" in args:
        i = args.index("--out")
        try:
            out = Path(args[i + 1])
        except IndexError:
            sys.exit("--out needs a folder")
        del args[i:i + 2]
        out.mkdir(parents=True, exist_ok=True)
    if not root.is_dir():
        sys.exit(f"no export folder at {root}\n"
                 "Run the experiments first, or pass --export DIR.")
    main(str(root), args or None, out)
