import scipy.io as sio, numpy as np, glob, csv, re
def wrap(a): return (np.asarray(a)+180)%360-180
def circ(deg):
    z=np.mean(np.exp(1j*np.radians(deg))); R=abs(z)
    return float(np.degrees(np.angle(z))), (float(np.degrees(np.sqrt(max(-2*np.log(R),0)))) if R>1e-12 else float('inf'))
def rd(f): d=sio.loadmat(f); return float(d['tb_dir'].squeeze()), float(d['tb_mag'].squeeze())
M=[('raft_sintel','RAFT'),('pwcnet_sintel','PWC-Net'),('liteflownet2_sintel','LFN2'),('me_attention','ME-Att'),('dual','Dual')]
f1=lambda x: (f"{x:.1f}").replace('-0.0','0.0')
def fm(x): return f"{x:.4f}" if abs(x)<0.1 else f"{x:.2f}"
def hdr(first, reps):
    h=" & ".join([r"\multicolumn{2}{c|}{%s}"%l for _,l in M])+r" & \multicolumn{2}{c}{stINRF}"
    sub=" & ".join([(r"$\hat\theta \pm$ SD & mag" if reps else r"$\hat\theta$ & mag")]*6)
    return f"{first} & {h} \\\\\n & {sub} \\\\"
def maerow(vals):
    cs=[r"\multicolumn{2}{c|}{$%.1f \pm %.1f$}"%v for v in vals[:-1]]+[r"\multicolumn{2}{c}{$%.1f \pm %.1f$}"%vals[-1]]
    return "MAE $\\pm$ SD & "+" & ".join(cs)+r" \\"
def table(sec, first, reps, rows, maes, cap, label, colspec="r|rr|rr|rr|rr|rr|rr"):
    return (f"% ----------------------------------------------------------------\n\\section{{{sec}}}\n% ----------------------------------------------------------------\n\n"
      "\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n\\resizebox{\\textwidth}{!}{%\n"
      f"\\begin{{tabular}}{{{colspec}}}\n\\toprule\n{hdr(first,reps)}\n\\midrule\n"+"\n".join(rows)+
      f"\n\\midrule\n{maerow(maes)}\n\\bottomrule\n\\end{{tabular}}\n}}\n\\caption{{{cap}}}\n\\label{{{label}}}\n\\end{{table}}\n\n")
out=[]

# ---------------- INTERACCION ----------------
sig=[10,25,37.5,100]; gt=[180,0,0,0]
I=np.loadtxt('stinrf_pertrial/interaction.csv',delimiter=',')      # 4 x 25, estimulo a la izquierda
rows=[]; errs={m:[] for m,_ in M}
for i,s in enumerate(sig):
    c=[f"${s:g}$ (${gt[i]}^\\circ$)"]
    for m,_ in M:
        v=[rd(f) for f in sorted(glob.glob(f'results_reps/DH/{m}/DH_CSr_sig{s:05.1f}_rep*/result.mat'))]
        d=np.array([x[0] for x in v]); mu,sd=circ(d); errs[m]+=list(np.abs(wrap(d-gt[i])))
        c+=[f"${f1(mu)} \\pm {sd:.1f}$", f"${fm(np.mean([x[1] for x in v]))}$"]
    x=I[i]; dd=np.where(x>0,0.0,180.0); mu,sd=circ(dd)
    c+=[f"${f1(mu%360 if mu>90 else mu)} \\pm {sd:.1f}$", f"${fm(np.abs(x).mean())}$"]
    rows.append(" & ".join(c)+r" \\")
st_ref=[0,0,0,180]; st_err=np.concatenate([np.abs(wrap(np.where(I[i]>0,0.0,180.0)-st_ref[i])) for i in range(4)])
maes=[(np.mean(errs[m]),np.std(errs[m])) for m,_ in M]+[(st_err.mean(),st_err.std())]
out.append(table("Interacción entre escalas","$\\sigma_t$ (ms)",True,rows,maes,
 "Condición CSr. DNNs y stINRF: $25$ runs por anchura, con la fase de cada componente aleatoria. En la columna stINRF el estímulo se desplaza hacia la izquierda, de modo que $0^\\circ$ es respuesta invertida y $180^\\circ$ respuesta correcta.","tab:res_dh"))

# ---------------- CONTRASTE ----------------
Cs=np.logspace(-2,0,20); stc=[float(r[1]) for r in csv.reader(open('stinrf_figures/export/Contrast/curve.csv'))]
rows=[]; errs={m:[] for m,_ in M}
for i,cc in enumerate(Cs):
    c=[f"${cc:.4f}$"]
    for m,_ in M:
        sf='0.500' if m=='raft_sintel' else '1.000'
        d,mg=rd(glob.glob(f'results_speedfix/{m}/Contrast_sf{sf}_c{i+1:02d}_*/result.mat')[0]); errs[m].append(abs(wrap(d)))
        c+=[f"${f1(d)}$", f"${fm(mg)}$"]
    c+=["$0.0$", f"${fm(stc[i])}$"]; rows.append(" & ".join(c)+r" \\")
maes=[(np.mean(errs[m]),np.std(errs[m])) for m,_ in M]+[(0.0,0.0)]
raft1=[abs(wrap(rd(glob.glob(f'results_speedfix/raft_sintel/Contrast_sf1.000_c{i+1:02d}_*/result.mat')[0])[0])) for i in range(20)]
t=table("Saturación por contraste","Contraste",False,rows,maes,
 "Los $20$ contrastes del código original. DNNs y stINRF: $1$ estímulo por contraste, con fase $0$; señal a $8^\\circ$/s.","tab:res_contrast")
t=t.replace(r"\multicolumn{2}{c|}{RAFT}",r"\multicolumn{2}{c|}{RAFT$^{*}$}",1)
out.append(t+"\\paragraph{*} RAFT aparece con frecuencia espacial de $0.5$~cpd, la condición en la que\nobtuvo su mejor resultado, mientras que los otros cuatro modelos lo hacen con $1$~cpd.\n\n")

# ---------------- MASKING ----------------
ratios=[0,0.02,0.05,0.2,0.25,0.5,1,3,4,6,12]
Gs=np.loadtxt('stinrf_pertrial/masking_sweep_gauss.csv',delimiter=','); M1=np.loadtxt('stinrf_pertrial/masking_1cpd.csv',delimiter=',')
sw={0:Gs[0],0.02:Gs[1],0.05:Gs[2],0.2:Gs[3],0.5:Gs[4],1:Gs[5],3:Gs[6],6:Gs[7],12:Gs[8],0.25:M1[1],4:M1[3]}
rows=[]; errs={m:[] for m,_ in M}; st_err=[]
for r in ratios:
    c=[f"${r:.2f}$ ($0^\\circ$)"]
    for m,_ in M:
        pat=f'Masking_nomask_rep*' if r==0 else f'Masking_maskr{r:05.2f}_rep*'
        v=[rd(f) for f in sorted(glob.glob(f'results_reps/Masking/{m}/{pat}/result.mat'))]
        d=np.array([x[0] for x in v]); mu,sd=circ(d); errs[m]+=list(np.abs(wrap(d)))
        c+=[f"${f1(mu)} \\pm {sd:.1f}$", f"${fm(np.mean([x[1] for x in v]))}$"]
    x=sw[r]; dd=np.where(x>0,0.0,180.0); mu,sd=circ(dd); st_err+=list(np.abs(wrap(dd)))
    c+=[f"${f1(mu)} \\pm {sd:.1f}$", f"${fm(np.abs(x).mean())}$"]; rows.append(" & ".join(c)+r" \\")
maes=[(np.mean(errs[m]),np.std(errs[m])) for m,_ in M]+[(np.mean(st_err),np.std(st_err))]
out.append(table("Masking","Ratio máscara",True,rows,maes,
 "Señal a $1$~cpd y $10^\\circ$/s; el ratio es la frecuencia espacial de la máscara dividida entre la de la señal ($0$ = sin máscara). DNNs y stINRF: $50$ runs por ratio, con fase aleatoria; la dirección de la máscara se sortea al azar en cada fotograma.","tab:res_masking"))

# ---------------- MISSING ----------------
ph0=np.loadtxt('stinrf_pertrial/contrast_missing_phase0.csv',delimiter=',')
rows=[]; errs={m:[] for m,_ in M}
for ty,gtv,lab,sti in [(1,0,'Onda cuadrada ($0^\\circ$)',5),(-1,180,'Missing fund.\\ ($180^\\circ$)',6)]:
    c=[lab]
    for m,_ in M:
        d,mg=rd(glob.glob(f'results_speedfix/{m}/Missing_sf1.000_type{ty:+d}/result.mat')[0]); errs[m].append(abs(wrap(d-gtv)))
        c+=[f"${f1(d)}$", f"${fm(mg)}$"]
    v=ph0[sti]; c+=[f"${0.0 if v>0 else 180.0:.1f}$", f"${fm(abs(v))}$"]; rows.append(" & ".join(c)+r" \\")
maes=[(np.mean(errs[m]),np.std(errs[m])) for m,_ in M]+[(0.0,0.0)]
out.append(table("Missing fundamental","Tipo",False,rows,maes,
 "DNNs y stINRF: $1$ estímulo por tipo, con fase $0$; $1$~cpd y $6^\\circ$/s ($24$~px por salto de $40$~ms en las DNNs).","tab:res_missing"))

# ---------------- SEGUNDO ORDEN ----------------
SO=np.loadtxt('stinrf_secondorder_trials.csv',delimiter=',')      # col0 derecha, col1 izquierda
rows=[]; errs={m:[] for m,_ in M}; st_err=[]
for ty,gtv,lab,col in [(1,0,'Derecha ($0^\\circ$)',0),(-1,180,'Izquierda ($180^\\circ$)',1)]:
    c=[lab]
    for m,_ in M:
        v=[rd(f) for f in sorted(glob.glob(f'results_reps/SecondOrder/{m}/SO_grain016_type{ty:+d}_rep*/result.mat'))]
        d=np.array([x[0] for x in v]); mu,sd=circ(d); errs[m]+=list(np.abs(wrap(d-gtv)))
        c+=[f"${f1(mu)} \\pm {sd:.1f}$", f"${fm(np.mean([x[1] for x in v]))}$"]
    x=SO[:,col]; dd=np.where(x>0,0.0,180.0); mu,sd=circ(dd); st_err+=list(np.abs(wrap(dd-gtv)))
    c+=[f"${f1(abs(mu))} \\pm {sd:.1f}$", f"${fm(np.abs(x).mean())}$"]; rows.append(" & ".join(c)+r" \\")
maes=[(np.mean(errs[m]),np.std(errs[m])) for m,_ in M]+[(np.mean(st_err),np.std(st_err))]
out.append(table("Movimiento de segundo orden","Dirección",True,rows,maes,
 "DNNs y stINRF: $100$ portadoras de ruido binario aleatorio por dirección, con la fase de la envolvente también aleatoria ($8^\\circ$/s). En las DNNs la portadora va en bloques de $16$~px. El MAE se calcula sobre las $200$ medidas individuales.","tab:res_secondorder"))

# ---------------- OPONENCIA ----------------
rows=[]; errs={m:[] for m,_ in M}
for dr,iv,lab in [(1,1,'$0^\\circ$, negra'),(1,0,'$0^\\circ$, blanca'),(-1,1,'$180^\\circ$, negra'),(-1,0,'$180^\\circ$, blanca')]:
    gtv=0 if dr>0 else 180; c=[lab]
    for m,_ in M:
        d,mg=rd(glob.glob(f'results_speedfix/{m}/Opponency_bw0.250_dir{dr:+d}_inv{iv}/result.mat')[0]); errs[m].append(abs(wrap(d-gtv)))
        c+=[f"${f1(d)}$", f"${fm(mg)}$"]
    c+=[f"${gtv:.1f}$","$87.18$"]; rows.append(" & ".join(c)+r" \\")
maes=[(np.mean(errs[m]),np.std(errs[m])) for m,_ in M]+[(0.0,0.0)]
out.append(table("Oponencia","Condición",False,rows,maes,
 "Barra de $0.25^\\circ$ a $8^\\circ$/s. DNNs y stINRF: $1$ estímulo por condición; el estímulo no tiene componentes aleatorias.","tab:res_opponency",colspec="l|rr|rr|rr|rr|rr|rr"))

# ---------------- SIMETRIA TEMPORAL ----------------
trs={}
for nm in ['LightEdges','Uneven']:
    R=[[float(x) for x in r] for r in csv.reader(open(f'stinrf_figures/export/{nm}/trials.csv'))]; trs[nm]=np.array(R)
lines=[]
for key,lab,stn in [('LightEdges','Light edges','LightEdges'),('UnevenSquareWave','Uneven sq.\\ wave','Uneven')]:
    first=True
    for m,ml in M+[('stinrf','stINRF')]:
        if m=='stinrf':
            fw=trs[stn][:,1]; rv=trs[stn][:,2]
            dF=np.where(fw>0,0.0,180.0); dR=np.where(rv>0,0.0,180.0); mF=np.abs(fw).mean(); mR=np.abs(rv).mean()
        else:
            vF=[rd(f) for f in sorted(glob.glob(f'results_reps/TimeReversal/{m}/TR_{key}_normal_ph*/result.mat'))]
            vR=[rd(f) for f in sorted(glob.glob(f'results_reps/TimeReversal/{m}/TR_{key}_reversed_ph*/result.mat'))]
            dF=np.array([x[0] for x in vF]); dR=np.array([x[0] for x in vR]); mF=np.mean([x[1] for x in vF]); mR=np.mean([x[1] for x in vR])
        muF,sdF=circ(dF); muR,sdR=circ(dR)
        e=np.concatenate([np.abs(wrap(dF)),np.abs(wrap(dR-180))])
        muR_s=f1(abs(muR)) if m=='stinrf' else f1(muR)
        pre=f"\\multirow{{6}}{{*}}{{{lab}}}" if first else ""; first=False
        lines.append(f"{pre} & {ml} & ${f1(muF)} \\pm {sdF:.1f}$ & ${fm(mF)}$ & ${muR_s} \\pm {sdR:.1f}$ & ${fm(mR)}$ & ${mF-mR:.2f}$ & ${e.mean():.1f} \\pm {e.std():.1f}$ \\\\")
    if key=='LightEdges': lines.append("\\midrule")
out.append("% ----------------------------------------------------------------\n\\section{Ruptura de la simetría temporal}\n% ----------------------------------------------------------------\n\n"
 "\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{ll|rr|rr|r|r}\n\\toprule\n"
 "Estímulo & Modelo & \\multicolumn{2}{c|}{Forward ($0^\\circ$)} & \\multicolumn{2}{c|}{Reversed ($180^\\circ$)} & $\\Delta$mag & MAE $\\pm$ SD \\\\\n"
 " & & $\\hat\\theta \\pm$ SD & mag & $\\hat\\theta \\pm$ SD & mag & fwd $-$ rev & (°) \\\\\n\\midrule\n"+"\n".join(lines)+
 "\n\\bottomrule\n\\end{tabular}\n}\n\\caption{DNNs y stINRF: $40$ fases aleatorias por estímulo y versión. $\\Delta$mag es la magnitud forward menos la reversed; la de stINRF está en unidades del modelo.}\n\\label{tab:res_simetria}\n\\end{table}\n\n")

doc=("\\documentclass[11pt,a4paper]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage[T1]{fontenc}\n\\usepackage[spanish]{babel}\n"
 "\\usepackage{amsmath}\n\\usepackage{booktabs}\n\\usepackage{float}\n\\usepackage{graphicx}\n\\usepackage{multirow}\n\\usepackage{geometry}\n\\geometry{margin=2.5cm}\n\n\\begin{document}\n\n"
 +"".join(out)+"\\end{document}\n")
open('resultados_repeticiones.tex','w').write(doc)
print('RAFT contraste a 1 cpd: MAE %.1f ± %.1f'%(np.mean(raft1),np.std(raft1)))
