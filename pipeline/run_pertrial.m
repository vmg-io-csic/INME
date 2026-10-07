% Model response of every individual run, for the experiments whose stimulus
% has a random spatial phase or a random noise sample.
out = inme_setup();
outd = out; if ~exist(outd,'dir'), mkdir(outd); end
N=400; meanL=100; sigma_c=57.74; j0=round(N/2); js=(j0-99):(j0+100);
rng(1);

%% CONTRASTE (50 runs x 5 contrastes)
cfg=stINRF3_config('Cont'); model=cfg.model; dov=4; fr=120; nf=round(500/1000*fr);
Cs=[0.010 0.032 0.100 0.316 1.000]; nt=50; R=zeros(numel(Cs),nt);
for si=1:numel(Cs)
  ph=rand(1,nt)*2*pi; r=zeros(1,nt); c=Cs(si);
  parfor tr=1:nt
    mI=make_Contrast_stimulus_1D(N,nf,fr,dov,100,ph(tr),8,1,meanL);
    mI=meanL*(1+c*2*(mI./abs(max(max(mI)))-0.5));
    O=stINRF3_compute_1D(mI,dov,model,fr,N,js); r(tr)=stINRF3_readout_gaussian_1D(O,js,sigma_c);
  end
  R(si,:)=r;
end
writematrix(R,[outd 'contrast.csv']); disp('contrast ok');

%% MASKING (50 runs x 4 ratios)
cfg=stINRF3_config('Cont'); cfg.model.sf0=2.5; model=cfg.model; dov=4; fr=120; nf=round(500/1000*fr);
ratios=[0 0.25 1 4]; nt=50; R=zeros(numel(ratios),nt);
for si=1:numel(ratios)
  ph=rand(1,nt)*2*pi; r=zeros(1,nt); msk=double(ratios(si)>0); sfm=max(ratios(si)*2.5,0.02);
  parfor tr=1:nt
    mI=make_Masking_stimulus_1D(N,nf,fr,dov,1000,ph(tr),10,2.5,sfm,meanL,msk,0.8,0.8,10);
    O=stINRF3_compute_1D(mI,dov,model,fr,N,js); r(tr)=stINRF3_readout_gaussian_1D(O,js,sigma_c);
  end
  R(si,:)=r;
end
writematrix(R,[outd 'masking.csv']); disp('masking ok');

%% MISSING FUNDAMENTAL (50 runs x 2 types: 1 square wave, -1 missing fundamental)
cfg=stINRF3_config('Cont'); model=cfg.model; dov=4; fr=120; nf=round(500/1000*fr);
Ty=[1 -1]; nt=50; R=zeros(2,nt);
for si=1:2
  ph=rand(1,nt)*2*pi; r=zeros(1,nt); ty=Ty(si);
  parfor tr=1:nt
    mI=make_Missing_stimulus_1D(N,nf,fr,dov,10000,ph(tr),6,1,meanL,ty);
    O=stINRF3_compute_1D(mI,dov,model,fr,N,js); r(tr)=stINRF3_readout_gaussian_1D(O,js,sigma_c);
  end
  R(si,:)=r;
end
writematrix(R,[outd 'missing.csv']); disp('missing ok');

%% REVERSE-PHI (20 runs x 2 tipos: 1 phi, -1 reverse)
cfg=stINRF3_config('Cont'); model=cfg.model; dov=16; fr=120; nf=round(500/1000*fr);
spd=4*(dov/N)*fr; nt=20; R=zeros(2,nt);
for si=1:2
  ph=rand(1,nt)*2*pi; r=zeros(1,nt); ty=Ty(si);
  parfor tr=1:nt
    mI=make_ReversePhi_stimulus_1D(N,nf,fr,dov,10000,ph(tr),spd,meanL,ty,0.8);
    O=stINRF3_compute_1D(mI,dov,model,fr,N,js); r(tr)=stINRF3_readout_gaussian_1D(O,js,sigma_c);
  end
  R(si,:)=r;
end
writematrix(R,[outd 'reversephi.csv']); disp('reversephi ok');

%% INTERACTION AMONG SCALES, CSr (25 runs x 4 envelope widths), leftward stimulus
cfg=stINRF3_config('DH'); model=cfg.model; dov=4; fr=1000; nf=500;
sig=[10 25 37.5 100]; nt=25; R=zeros(numel(sig),nt);
for si=1:numel(sig)
  ph1=rand(1,nt)*2*pi; ph2=rand(1,nt)*2*pi; r=zeros(1,nt); s=sig(si);
  parfor tr=1:nt
    mI=make_Interaction_stimulus_1D(N,nf,fr,dov,4,s,ph1(tr),ph2(tr),12/3,12/1,3,1,meanL,-1);
    O=stINRF3_compute_1D(mI,dov,model,fr,N,js); r(tr)=stINRF3_readout_gaussian_1D(O,js,sigma_c);
  end
  R(si,:)=r;
end
writematrix(R,[outd 'interaction.csv']); disp('interaction ok');
disp('ALL DONE');
