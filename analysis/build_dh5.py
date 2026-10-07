import numpy as np, glob, re
from pathlib import Path
# Reuse the loading helpers defined at the top of build_tables.py.
exec(open(Path(__file__).resolve().parent / 'build_tables.py').read().split('out=[]')[0])
sig=[10,25,37.5,100]; CN=['S1','S3','CS','CSr','CM']
LAB={'S1':'S1','S3':'S3','CS':'CS','CSr':'CSr','CM':'CM'}
def gt(c,i): return 180 if (c=='CSr' and i==0) else 0
rows=[]; allerr={m:[] for m,_ in M+[('stinrf','')]}; conderr={c:{m:[] for m,_ in M+[('stinrf','')]} for c in CN}
for c in CN:
    base='results_reps/DH' if c=='CSr' else 'results_reps/DH_conds'
    St=np.loadtxt('stinrf_pertrial/'+('interaction.csv' if c=='CSr' else f'interaction_{c}.csv'),delimiter=',')
    for i,s in enumerate(sig):
        pre=f"\\multirow{{4}}{{*}}{{{LAB[c]}}}" if i==0 else ""
        cells=[pre, f"${s:g}$ (${gt(c,i)}^\\circ$)"]
        for m,_ in M:
            v=[rd(f) for f in sorted(glob.glob(f'{base}/{m}/DH_{c}_sig{s:05.1f}_rep*/result.mat'))]
            assert len(v)==25,(c,s,m,len(v))
            d=np.array([x[0] for x in v]); mu,sd=circ(d); e=list(np.abs(wrap(d-gt(c,i))))
            allerr[m]+=e; conderr[c][m]+=e
            cells+=[f"${f1(mu)} \\pm {sd:.1f}$", f"${fm(np.mean([x[1] for x in v]))}$"]
        x=St[i]; dd=np.where(x<0,0.0,180.0); mu,sd=circ(dd); mu=abs(mu) if abs(abs(mu)-180)<1e-6 else mu
        e=list(np.abs(wrap(dd-gt(c,i)))); allerr['stinrf']+=e; conderr[c]['stinrf']+=e
        cells+=[f"${f1(mu)} \\pm {sd:.1f}$", f"${fm(np.abs(x).mean())}$"]
        rows.append(" & ".join(cells)+r" \\")
    if c!=CN[-1]: rows.append("\\midrule")
keys=[m for m,_ in M]+['stinrf']
mae=lambda L:(np.mean(L),np.std(L))
hdr=" & ".join([r"\multicolumn{2}{c|}{%s}"%l for _,l in M])+r" & \multicolumn{2}{c}{stINRF}"
sub=" & ".join([r"$\hat\theta \pm$ SD & mag"]*6)
maerow="\\multicolumn{2}{l|}{MAE $\\pm$ SD} & "+" & ".join([r"\multicolumn{2}{c|}{$%.1f \pm %.1f$}"%mae(allerr[k]) for k in keys[:-1]]+[r"\multicolumn{2}{c}{$%.1f \pm %.1f$}"%mae(allerr['stinrf'])])+r" \\"
t1=("\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{ll|rr|rr|rr|rr|rr|rr}\n\\toprule\n"
 f"Condición & $\\sigma_t$ (ms) & {hdr} \\\\\n & & {sub} \\\\\n\\midrule\n"+"\n".join(rows)+f"\n\\midrule\n{maerow}\n\\bottomrule\n\\end{{tabular}}\n}}\n"
 "\\caption{Condiciones: S1, solo la componente de $1$~cpd en movimiento; S3, solo la de $3$~cpd; CS, $1$~cpd en movimiento $+$ $3$~cpd estática; CSr, $3$~cpd en movimiento $+$ $1$~cpd estática; CM, ambas en movimiento. Las componentes se mueven a $12$~Hz. Entre paréntesis, la dirección esperada: $0^\\circ$ es la dirección del movimiento. DNNs y stINRF: $25$ runs por condición y duración, con la fase de cada componente aleatoria. En stINRF el estímulo se desplaza hacia la izquierda y su respuesta se expresa respecto a la dirección del movimiento. El MAE se calcula sobre los $500$ runs de las $5$ condiciones.}\n"
 "\\label{tab:res_dh}\n\\end{table}\n")
crow=[]
for c in CN:
    crow.append(f"{c} & "+" & ".join(["$%.1f \\pm %.1f$"%mae(conderr[c][k]) for k in keys])+r" \\")
t2=("\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{l|rrrrrr}\n\\toprule\n"
 "Condición & RAFT & PWC-Net & LFN2 & ME-Att & Dual & stINRF \\\\\n\\midrule\n"+"\n".join(crow)+"\n\\bottomrule\n\\end{tabular}\n"
 "\\caption{MAE $\\pm$ SD (°) por condición de interacción entre escalas, sobre los $100$ runs de cada condición ($4$ duraciones $\\times$ $25$ runs).}\n\\label{tab:res_dh_mae}\n\\end{table}\n")
t=(t1+"\n"+t2).replace('\\pm -0.0','\\pm 0.0')
open('dh5_tables.tex','w').write(t)
for c in CN: print(c,[ "%.1f"%mae(conderr[c][k])[0] for k in keys])
print('ALL',["%.1f"%mae(allerr[k])[0] for k in keys])
