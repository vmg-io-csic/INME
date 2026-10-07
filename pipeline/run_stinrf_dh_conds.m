% Model response for the interaction among scales, conditions S1, S3, CS and CM:
% 25 runs per envelope width, following the same protocol as the CSr condition.
ps=parallel.Settings; ps.Pool.AutoCreate=true;
out = inme_setup();
outd=out;
N=400; meanL=100; sigma_c=57.74; j0=round(N/2); js=(j0-99):(j0+100);
cfg=stINRF3_config('DH'); model=cfg.model; dov=4; fr=1000; nf=500;
sig=[10 25 37.5 100]; nt=25; CN={'S1','S3','CS','CSr','CM'};
rng(11);
for ci=[1 2 3 5]
  R=zeros(numel(sig),nt);
  for si=1:numel(sig)
    ph1=rand(1,nt)*2*pi; ph2=rand(1,nt)*2*pi; r=zeros(1,nt); s=sig(si);
    parfor tr=1:nt
      mI=make_Interaction_stimulus_1D(N,nf,fr,dov,ci,s,ph1(tr),ph2(tr),12/3,12/1,3,1,meanL,-1);
      O=stINRF3_compute_1D(mI,dov,model,fr,N,js); r(tr)=stINRF3_readout_gaussian_1D(O,js,sigma_c);
    end
    R(si,:)=r;
  end
  writematrix(R,[outd 'interaction_' CN{ci} '.csv']); fprintf('interaction %s ok\n',CN{ci});
end
disp('INTERACTION CONDITIONS DONE');
