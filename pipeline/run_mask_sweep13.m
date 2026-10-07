% Barrido de masking en stINRF3: 13 frecuencias de mascara (las 8 de Programs_ISP
% plus 0.25, 1.5, 2, 2.5 and 4 cpd, and the no-mask condition; 1 cpd signal,
% 50 trials.
% Both readouts are saved: the Gaussian pooling and the plain sensor mean.
cfg=stINRF3_config('Cont'); model=cfg.model;
N=400; dov=4; fr=120; nf=round(500/1000*fr); meanL=100; sigma_c=57.74;
j0=round(N/2); js=(j0-99):(j0+100); rng(3);
sfmask=[0.02 0.05 0.2 0.25 0.5 1 1.5 2 2.5 3 4 6 12]; nt=50;
G=zeros(numel(sfmask)+1,nt); S=G;
for si=0:numel(sfmask)
  ph=rand(1,nt)*2*pi; g=zeros(1,nt); s=zeros(1,nt);
  if si==0, msk=0; sfm=0.02; else, msk=1; sfm=sfmask(si); end
  for tr=1:nt
    mI=make_Masking_stimulus_1D(N,nf,fr,dov,1000,ph(tr),10,1,sfm,meanL,msk,0.8,0.8,10);
    O=stINRF3_compute_1D(mI,dov,model,fr,N,js);
    g(tr)=stINRF3_readout_gaussian_1D(O,js,sigma_c); s(tr)=mean(O);
  end
  G(si+1,:)=g; S(si+1,:)=s;
end
writematrix(G,'stinrf_pertrial/masking_sweep13_gauss.csv');
writematrix(S,'stinrf_pertrial/masking_sweep13_meansensors.csv');
disp('SWEEP DONE');
