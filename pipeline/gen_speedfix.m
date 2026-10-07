out = inme_setup();
G = struct('tag',{'raft','pwc_lfn2','meattn_dual'},'N_cols',{960,1024,768},'N_rows',{520,436,768});
meanL=100;
%% 1) MASKING, remaining repetitions: 10 deg/s, 50 phases x 11 mask frequencies
for g=1:numel(G)
  N=G(g).N_cols; NR=G(g).N_rows; c=dnn_v25_config('Masking',N); fr=c.framerate; dov=N/c.px_per_deg;
  out=fullfile(fullfile(out,'stimuli','Masking'),G(g).tag,'v25');
  for r=[0 c.mask_ratios], for rep=1:50
    if r==0, msk=0; sfm=1; nm=sprintf('Masking_nomask_rep%02d',rep); else, msk=1; sfm=r; nm=sprintf('Masking_maskr%05.2f_rep%02d',r,rep); end
    rng(6000+rep); pm=rand*2*pi;
    mI=make_Masking_stimulus_1D(N,c.n_frames,fr,dov,c.sigma_t_ms,pm,c.speed_dps,1,sfm,meanL,msk,c.contrast_signal,c.contrast_mask,c.tf_mask_hz);
    x=struct('experiment','Masking_reps50','mask_ratio',r,'masking',msk,'sfm_cpd',1,'sfmask_cpd',sfm,'rep',rep,'phase_mov',pm, ...
      'px_per_deg',c.px_per_deg,'framerate',fr,'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'speed_dps',c.speed_dps, ...
      'speed_px_per_frame',c.speed_dps*c.px_per_deg/fr,'sigma_t_ms',c.sigma_t_ms,'gt_direction_deg',0);
    dnn_v25_save(fullfile(out,nm),dnn_v25_broadcast(mI,NR),x);
  end, end
  fprintf('%s masking 10 deg/s\n',G(g).tag);
end
%% 2) CONTRASTE (20), OPONENCIA (barra 0.25), MISSING (1 cpd), fase 0
for g=1:numel(G)
  N=G(g).N_cols; NR=G(g).N_rows; tag=G(g).tag; out=fullfile(fullfile(out,'stimuli','speedfix'),tag,'v25');
  c=dnn_v25_config('Contrast',N); dov=N/c.px_per_deg; sfs=1; if strcmp(tag,'raft'), sfs=[0.5 1]; end
  for sf=sfs, for ci=1:numel(c.contrast_list)
    con=c.contrast_list(ci);
    mI=make_Contrast_stimulus_1D(N,c.n_frames,c.framerate,dov,c.sigma_t_ms,0,c.speed_dps,sf,meanL);
    mI=meanL*(1+con*2*(mI./abs(max(mI(:)))-0.5));
    x=struct('experiment','Contrast_v8','sfm_cpd',sf,'contrast_level',con,'phase_mov',0,'px_per_deg',c.px_per_deg,'framerate',c.framerate, ...
      'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'speed_dps',c.speed_dps,'speed_px_per_frame',c.speed_dps*c.px_per_deg/c.framerate, ...
      'sigma_t_ms',c.sigma_t_ms,'gt_direction_deg',0);
    dnn_v25_save(fullfile(out,sprintf('Contrast_sf%.3f_c%02d_%.4f',sf,ci,con)),dnn_v25_broadcast(mI,NR),x);
  end, end
  c=dnn_v25_config('Opponency',N); dov=N/c.px_per_deg; bw=0.25;
  for iv=c.inver, for dr=c.directions
    mI=make_Opponency_stimulus_1D_v25(N,c.n_frames,c.framerate,dov,c.sigma_t_ms,c.speed_dps,bw,meanL,dr,iv);
    mI=meanL*(1+2*(mI./abs(max(mI(:)))-0.5));
    x=struct('experiment','Opponency_v8','direction',dr,'bar_width',bw,'bar_width_px',bw*c.px_per_deg,'inver',iv, ...
      'px_per_deg',c.px_per_deg,'framerate',c.framerate,'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'speed_dps',c.speed_dps, ...
      'speed_px_per_frame',c.speed_dps*c.px_per_deg/c.framerate,'sigma_t_ms',c.sigma_t_ms,'gt_direction_deg',(dr<0)*180);
    dnn_v25_save(fullfile(out,sprintf('Opponency_bw%.3f_dir%+d_inv%d',bw,dr,iv)),dnn_v25_broadcast(mI,NR),x);
  end, end
  c=dnn_v25_config('Missing',N); dov=N/c.px_per_deg;
  for ty=c.types
    mI=make_Missing_stimulus_1D(N,c.n_frames,c.framerate,dov,c.sigma_t_ms,0,c.speed_dps,1,meanL,ty);
    x=struct('experiment','Missing_v6','type',ty,'sfm_cpd',1,'phase_mov',0,'px_per_deg',c.px_per_deg,'framerate',c.framerate, ...
      'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'speed_dps',c.speed_dps,'px_per_step',c.speed_dps*0.040*c.px_per_deg, ...
      'sigma_t_ms',c.sigma_t_ms,'gt_direction_deg',(ty<0)*180);
    dnn_v25_save(fullfile(out,sprintf('Missing_sf1.000_type%+d',ty)),dnn_v25_broadcast(mI,NR),x);
  end
  fprintf('%s contraste/oponencia/missing\n',tag);
end
disp('SPEEDFIX GENERADO');
