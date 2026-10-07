"""Renderiza los estímulos del benchmark a mp4 (H.264) desde los frames.mat."""
import subprocess, csv
import os
from pathlib import Path

# Repository root, resolved from this file so the script runs from anywhere.
REPO = Path(__file__).resolve().parent.parent
# Folder holding stimuli and results. Override with the INME_OUT variable.
OUT_ROOT = Path(os.environ.get("INME_OUT", REPO / "output"))

import numpy as np, h5py

D = str(OUT_ROOT) + os.sep
OUT = str(OUT_ROOT / 'supp_videos') + os.sep
os.makedirs(OUT, exist_ok=True)
PLAY=25                      # fps de reproducción, igual para todos

def tensor(path):
    with h5py.File(path+'/frames.mat') as h:
        t=np.array(h['tensor']).squeeze()
        meta={k: float(np.array(h[k]).squeeze()) for k in h.keys()
              if k!='tensor' and np.array(h[k]).size==1 and np.array(h[k]).dtype.kind in 'fiu'}
    if t.ndim!=3: raise ValueError(path+' ndim %d'%t.ndim)
    # h5py entrega (T, W, H) -> (T, H, W)
    return np.transpose(t,(0,2,1)), meta

def write(name, frames, fps_real, title):
    v=frames.astype(np.float64)
    v=(v-v.min())/max(v.max()-v.min(),1e-9)*255.0
    v=v.astype(np.uint8)
    T,H,W=v.shape
    reps=max(1,int(np.ceil(2.5*PLAY/T)))     # clips cortos: repetir la secuencia
    v=np.tile(v,(reps,1,1))
    sw, sh = (W+1)//2*2, (H+1)//2*2
    p=subprocess.Popen(['ffmpeg','-y','-loglevel','error','-f','rawvideo','-pix_fmt','gray',
        '-s','%dx%d'%(W,H),'-r',str(PLAY),'-i','-','-vf','pad=%d:%d:0:0:color=gray'%(sw,sh),
        '-c:v','libx264','-preset','slow','-crf','18','-pix_fmt','yuv420p',OUT+name+'.mp4'],
        stdin=subprocess.PIPE)
    p.communicate(v.tobytes())
    slow=fps_real/PLAY
    return dict(file=name+'.mp4', title=title, frames=T, size='%dx%d'%(W,H),
                fps_stimulus='%.2f'%fps_real, fps_playback=PLAY,
                slowdown='%.1f'%slow, repeats=reps)

JOBS=[]
def add(name,path,title,fps=None):
    JOBS.append((name,path,title,fps))

V=D+'v25_meattn/stimuli_by_model/raft/v25/'
S=D+'stim_speedfix/raft/v25/'
R=D+'stim_reps/'
T16=D+'stim_v16/stimuli_by_model/raft/v25/'
for d in range(0,91,15):
    add('VideoS_gratings_dir%03d'%d, V+'gratings_dir%03d'%d, 'Grating drifting at %d deg'%d)
    add('VideoS_plaidI_dir%03d'%d,  V+'plaidI_dir%03d'%d,  'Type I plaid, pattern direction %d deg'%d)
    add('VideoS_texture_dir%03d'%d, T16+'gratingTex_dir%03d_ct120'%d, 'Textured grating, %d deg'%d)
for ioc in range(-30,31,10):
    import glob
    g=sorted(glob.glob(V+'plaidII_ioc%+04d_c*'%ioc))
    if g: add('VideoS_plaidII_ioc%+03d'%ioc, g[0], 'Type II plaid, IOC %+d deg'%ioc)
for ci,c in [(1,0.01),(10,0.0886),(20,1.0)]:
    add('VideoS_contrast_c%05.3f'%c, S+'Contrast_sf1.000_c%02d_%.4f'%(ci,c), 'Contrast saturation, contrast %.3f'%c)
for ty,lab in [(1,'square_wave'),(-1,'missing_fundamental')]:
    add('VideoS_missing_%s'%lab, S+'Missing_sf1.000_type%+d'%ty, 'Missing fundamental: %s'%lab.replace('_',' '))
for dr in [1,-1]:
    for inv in [0,1]:
        add('VideoS_opponency_dir%+d_inv%d'%(dr,inv), S+'Opponency_bw0.250_dir%+d_inv%d'%(dr,inv),
            'Opponency, %s bar moving %s'%('white' if inv else 'black','right' if dr>0 else 'left'))
for r in [0,0.02,0.05,0.2,0.25,0.5,1,1.5,2,2.5,3,4,6,12]:
    nm='Masking_nomask_rep01' if r==0 else 'Masking_maskr%05.2f_rep01'%r
    add('VideoS_masking_%s'%('nomask' if r==0 else 'mask%05.2fcpd'%r), R+'Masking/raft/v25/'+nm,
        'Motion masking, %s'%('no mask' if r==0 else 'mask at %g cpd'%r))
for ty,lab in [(1,'rightward'),(-1,'leftward')]:
    add('VideoS_secondorder_%s'%lab, R+'SecondOrder/raft/v25/SO_grain016_type%+d_rep001'%ty,
        'Second-order motion, envelope moving %s'%lab)
    add('VideoS_%s'%('phi' if ty>0 else 'reversephi'), R+'ReversePhi/raft/v25/RP_grain016_type%+d_rep01'%ty,
        'Phi motion' if ty>0 else 'Reverse-phi motion')
for s in [10,25,37.5,100]:
    add('VideoS_interaction_sigma%05.1fms'%s, R+'DH/raft/v25/DH_CSr_sig%05.1f_rep01'%s,
        'Interaction among scales, CSr, sigma_t = %g ms'%s)
for stim,base in [('LightEdges','LightEdges_v2'),('UnevenSquareWave','Uneven_v2')]:
    for rev in ['normal','reversed']:
        add('VideoS_timereversal_%s_%s'%(stim.lower(),rev), R+'%s/raft/v25/TR_%s_%s_ph01'%(base,stim,rev),
            '%s, %s playback'%(stim,rev))

rows=[]
for name,path,title,fps in JOBS:
    if not os.path.isdir(path): print('FALTA', path); continue
    t,meta=tensor(path)
    f=fps or meta.get('framerate') or 25
    rows.append(write(name,t,f,title))
    print(name, rows[-1]['size'], rows[-1]['frames'],'frames', rows[-1]['fps_stimulus'],'fps')
with open(OUT+'video_list.csv','w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('TOTAL', len(rows))
