% The model's grating experiment, with the texture used for the networks.
out = inme_setup();
N=300; dov=300*(1.5/128); framerate=200; duration_ms=2000; meanL=100;
dirs=0:15:90; sf=0.25; hz=10; phase=0; half_region=100;
sigma_c=2*half_region/sqrt(12);
c_grating=0.4; c_tex=[0 0.04 0.08 0.12]; plo=40; phi=60; seed_tex=2;
n_frames=round(duration_ms/1000*framerate);
cfg=stINRF3_config('gratings'); model=cfg.model;
f_cpx=sf/(N/dov); sp=hz/framerate;
nd=numel(dirs); nc=numel(c_tex);
[ic_g,id_g]=ndgrid(1:nc,1:nd); icv=ic_g(:); idv=id_g(:); nj=numel(icv);
D=zeros(nj,1); M=zeros(nj,1); CL=false(nj,1); PK=zeros(nj,1); EX=zeros(nj,1);
parfor jj=1:nj
    [TCB,info]=make_textured_grating(N,n_frames,f_cpx,sp,dirs(idv(jj)),phase, ...
        c_grating,c_tex(icv(jj)),meanL,plo,phi,seed_tex);
    CL(jj)=info.clipped; PK(jj)=info.peak_texture; EX(jj)=info.max_excursion;
    [H,V]=stINRF3_compute(TCB,dov,model,framerate,N);
    [d,m]=stINRF3_readout_gaussian(H,V,sigma_c,'SupportHalf',half_region);
    D(jj)=d; M(jj)=m;
end
T=table(c_tex(icv)',dirs(idv)',D,M,CL,PK,EX,'VariableNames', ...
    {'c_texture','gt','dir','mag','clipped','peak_texture','max_excursion'});
writetable(T,fullfile(out,'stinrf_gratings_textura.csv'));
fprintf('recortados: %d de %d   pico textura max %.2f   excursion max %.3f\n', sum(CL), nj, max(PK), max(EX));
disp('DONE');
