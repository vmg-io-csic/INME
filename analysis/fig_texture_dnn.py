"""Predicted vs stimulus direction for the five DNNs: original grating vs
textured grating, one figure in the style of Fig_Gratings."""
import sys, csv
import os
from pathlib import Path

# Repository root, resolved from this file so the script runs from anywhere.
REPO = Path(__file__).resolve().parent.parent
# Folder holding stimuli and results. Override with the INME_OUT variable.
OUT_ROOT = Path(os.environ.get("INME_OUT", REPO / "output"))

import numpy as np
import scipy.io as sio
sys.path.insert(0, str(REPO / 'model_stinrf'))
import plot_figures as pf
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

D = str(OUT_ROOT) + os.sep
DIRS = np.arange(0, 91, 15)
MODELS = [('raft_sintel', 'RAFT', '#E69F00', 'o'),
          ('pwcnet_sintel', 'PWC-Net', '#56B4E9', 's'),
          ('liteflownet2_sintel', 'LiteFlowNet2', '#009E73', '^'),
          ('me_attention', 'ME-Attention', '#CC79A7', 'D'),
          ('dual', 'Dual', '#D55E00', 'v')]
wrap = lambda a: (np.asarray(a) + 180) % 360 - 180

# original grating: c = 1.00, 34.13 px/frame (results_v25_long.csv)
orig = {}
for r in csv.DictReader(open(D + 'results_v25_long.csv')):
    if r['stimulus'].startswith('gratings_dir'):
        orig.setdefault(r['model'], {})[int(r['stimulus'][-3:])] = float(r['tb_dir'])
orig = {m: np.array([orig[m][d] for d in DIRS]) for m, *_ in MODELS}
# textured grating: c = 0.40 + texture 0.12 RMS, 16 px/frame (results_v16)
tex = {m: np.array([float(sio.loadmat(f'{D}results_v16/{m}/gratingTex_dir{d:03d}_ct120/result.mat')['tb_dir'].squeeze())
                    for d in DIRS]) for m, *_ in MODELS}

def unwrap_to(y, x):
    """Express each prediction as the angle closest to its ground truth, so a
    direction of -177 deg for a 75 deg target is drawn at 183 deg rather than
    across the whole axis. Values stay within [x-180, x+180]."""
    return x + wrap(y - x)

wfig = 7.2
mar_l, mar_r, gap = 0.62, 0.12, 0.42
pw = (wfig - mar_l - mar_r - gap) / 2
ph, mar_t, title_h = 2.75, 0.10, 0.36
xlab_h, icon_in, igap, mar_b = 0.40, 0.33, 0.04, 0.06
hfig = mar_b + icon_in + igap + xlab_h + ph + title_h + mar_t
fig = plt.figure(figsize=(wfig, hfig), facecolor='w')
xlim = (-9, 99); ylim = (-195, 195)
py = mar_b + icon_in + igap + xlab_h

panels = [('Original grating', 'contrast 1.00, 34.1 px/frame', orig),
          ('Textured grating', 'contrast 0.40 + texture, 16 px/frame', tex)]
for k, (title, sub, data) in enumerate(panels):
    px = mar_l + k * (pw + gap)
    ax = fig.add_axes([px / wfig, py / hfig, pw / wfig, ph / hfig])
    ax.plot([-9, 99], [-9, 99], '--', color=pf.REF, lw=1.0, zorder=1)
    handles = []
    for m, lab, col, mk in MODELS:
        y = data[m]
        mae = np.abs(wrap(y - DIRS)).mean()
        ax.plot(DIRS, y, '-', color=col, lw=1.4, zorder=2)
        ax.plot(DIRS, y, mk, color=col, ms=4.2, mec='w', mew=0.5, zorder=3)
        handles.append(Line2D([], [], color=col, marker=mk, ms=4.2, mec='w', mew=0.5, lw=1.4,
                              label=f'{lab}  ({mae:.1f}{pf.DEG})'))
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_xticks(DIRS); ax.set_yticks(np.arange(-180, 181, 45))
    ax.tick_params(direction='out', labelsize=8.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel(f'Stimulus direction ({pf.DEG})', fontsize=9.5)
    if k == 0:
        ax.set_ylabel(f'Predicted direction ({pf.DEG})', fontsize=9.5)
    else:
        ax.set_yticklabels([])
    leg = ax.legend(handles=handles, loc='lower left', fontsize=7.4, frameon=False,
                    title=f'Model (MAE)', title_fontsize=7.6, handlelength=2.0,
                    borderaxespad=0.2, labelspacing=0.35)
    leg._legend_box.align = 'left'
    fig.text((px + pw / 2) / wfig, (py + ph + title_h * 0.52) / hfig, title,
             fontsize=10.5, fontweight='bold', ha='center', va='bottom')
    fig.text((px + pw / 2) / wfig, (py + ph + title_h * 0.12) / hfig, sub,
             fontsize=8.5, color=pf.TEXT, ha='center', va='bottom')
    fx = (DIRS - xlim[0]) / (xlim[1] - xlim[0])
    for i, d in enumerate(DIRS):
        a = fig.add_axes([(px + fx[i] * pw - icon_in / 2) / wfig, mar_b / hfig,
                          icon_in / wfig, icon_in / hfig])
        pf.stim_icon(a, [d], 2.2, [pf.INK])

out = D + 'Fig_TextureGratings_DNN'
fig.savefig(out + '.png', dpi=600, facecolor='w')
fig.savefig(out + '.pdf', facecolor='w')
for m, lab, *_ in MODELS:
    print(lab, 'orig MAE %.2f' % np.abs(wrap(orig[m] - DIRS)).mean(), ' tex MAE %.2f' % np.abs(wrap(tex[m] - DIRS)).mean())
