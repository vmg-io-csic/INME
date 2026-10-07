out = inme_setup();
G = struct('tag',{'raft','pwc_lfn2','meattn_dual'},'N_cols',{960,1024,768},'N_rows',{520,436,768});
meanL=100;
for g=1:numel(G)
  N=G(g).N_cols; NR=G(g).N_rows; c=dnn_v25_config('Masking',N); fr=c.framerate; dov=N/c.px_per_deg;
  out=fullfile(fullfile(out,'stimuli','Masking_extra'),G(g).tag,'v25');
  for r=[1.5 2 2.5], for rep=1:50
    msk=1; sfm=r; nm=sprintf('Masking_maskr%05.2f_rep%02d',r,rep);
    rng(6000+rep); pm=rand*2*pi;
    mI=make_Masking_stimulus_1D(N,c.n_frames,fr,dov,c.sigma_t_ms,pm,c.speed_dps,1,sfm,meanL,msk,c.contrast_signal,c.contrast_mask,c.tf_mask_hz);
    x=struct('experiment','Masking_reps50','mask_ratio',r,'masking',msk,'sfm_cpd',1,'sfmask_cpd',sfm,'rep',rep,'phase_mov',pm, ...
      'px_per_deg',c.px_per_deg,'framerate',fr,'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'speed_dps',c.speed_dps, ...
      'speed_px_per_frame',c.speed_dps*c.px_per_deg/fr,'sigma_t_ms',c.sigma_t_ms,'gt_direction_deg',0);
    dnn_v25_save(fullfile(out,nm),dnn_v25_broadcast(mI,NR),x);
  end, end
  fprintf('%s masking extra generado\n',G(g).tag);
end
disp('GENERACION TERMINADA');
