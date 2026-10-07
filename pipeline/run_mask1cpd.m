% Masking with the signal at 1 cpd, the default detector frequency and 50
% trials, read out with the Gaussian pooling.
out = inme_setup();
cfg=stINRF3_config('Cont'); model=cfg.model;
N=400; dov=4; fr=120; nf=round(500/1000*fr); meanL=100; sigma_c=57.74;
j0=round(N/2); js=(j0-99):(j0+100); rng(2);
ratios=[0 0.25 1 4]; sf_sig=1; nt=50; R=zeros(numel(ratios),nt);
for si=1:numel(ratios)
  ph=rand(1,nt)*2*pi; r=zeros(1,nt); msk=double(ratios(si)>0); sfm=max(ratios(si)*sf_sig,0.02);
  parfor tr=1:nt
    mI=make_Masking_stimulus_1D(N,nf,fr,dov,1000,ph(tr),10,sf_sig,sfm,meanL,msk,0.8,0.8,10);
    O=stINRF3_compute_1D(mI,dov,model,fr,N,js); r(tr)=stINRF3_readout_gaussian_1D(O,js,sigma_c);
  end
  R(si,:)=r;
end
writematrix(R,fullfile(out,'masking_1cpd.csv')); disp('written');
