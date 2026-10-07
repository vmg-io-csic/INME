% Model response with the spatial phase fixed at zero, for contrast saturation
% and the missing fundamental, with the same parameters as the random-phase runs.
out = inme_setup();
N=400; dov=4; fr=120; nf=round(500/1000*fr); meanL=100; sigma_c=57.74;
j0=round(N/2); js=(j0-99):(j0+100);
cfg=stINRF3_config('Cont'); model=cfg.model;
Cs=[0.010 0.032 0.100 0.316 1.000]; RC=zeros(1,numel(Cs));
for si=1:numel(Cs)
    mI=make_Contrast_stimulus_1D(N,nf,fr,dov,100,0,8,1,meanL);
    mI=meanL*(1+Cs(si)*2*(mI./abs(max(max(mI)))-0.5));
    O=stINRF3_compute_1D(mI,dov,model,fr,N,js); RC(si)=stINRF3_readout_gaussian_1D(O,js,sigma_c);
end
Ty=[1 -1]; RM=zeros(1,2);
for si=1:2
    mI=make_Missing_stimulus_1D(N,nf,fr,dov,10000,0,6,1,meanL,Ty(si));
    O=stINRF3_compute_1D(mI,dov,model,fr,N,js); RM(si)=stINRF3_readout_gaussian_1D(O,js,sigma_c);
end
fprintf('CONTRASTE fase0: %s\n', sprintf('%+.5f ', RC));
fprintf('MISSING fase0: cuadrada %+.4f  sin fundamental %+.4f\n', RM(1), RM(2));
writematrix([RC RM],fullfile(out,'contrast_missing_phase0.csv'));
disp('ZERO-PHASE RUN DONE');
