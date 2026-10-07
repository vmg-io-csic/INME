% Interaction among scales: conditions S1, S3, CS and CM, with 25 pairs of
% random spatial phases each (the same seeds as CSr).
out = inme_setup();
G = struct('tag',{'raft','pwc_lfn2','meattn_dual'},'N_cols',{960,1024,768},'N_rows',{520,436,768});
ROOT=fullfile(out,'stimuli'); meanL=100;
CN={'S1','S3','CS','CSr','CM'};
for g=1:numel(G)
  N=G(g).N_cols; NR=G(g).N_rows; tag=G(g).tag;
  c=dnn_v25_config('DH',N); dov=N/c.px_per_deg; out=fullfile(ROOT,'DH_conds',tag,'v25');
  for ci=[1 2 3 5], for s=c.sigma_t_ms, for rep=1:25
    e=fullfile(out,sprintf('DH_%s_sig%05.1f_rep%02d',CN{ci},s,rep)); if exist(fullfile(e,'frames.mat'),'file'), continue; end
    rng(5000+rep); pm=rand*2*pi; pst=rand*2*pi;
    mI=make_DH_stimulus_1D_v25(N,c.n_frames,c.framerate,dov,ci,s,pm,pst,c.tf_hz,c.sfm,c.sfs,meanL);
    x=struct('experiment','DH_reps25','condition',CN{ci},'cond_idx',ci,'sigma_t_ms',s,'rep',rep,'phase_mov',pm,'phase_stat',pst, ...
      'px_per_deg',c.px_per_deg,'framerate',c.framerate,'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'tf_hz',c.tf_hz);
    dnn_v25_save(e,dnn_v25_broadcast(mI,NR),x);
  end, end, end
  fprintf('%s interaction conditions generated\n',tag);
end
disp('INTERACTION CONDITIONS GENERATED');
