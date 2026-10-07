% Stimuli for the networks, with the same number of random repetitions used
% for the model.
out = inme_setup();
G = struct('tag',{'raft','pwc_lfn2','meattn_dual'},'N_cols',{960,1024,768},'N_rows',{520,436,768});
ROOT=fullfile(out,'stimuli'); meanL=100;
for g=1:numel(G)
  N=G(g).N_cols; NR=G(g).N_rows; tag=G(g).tag;

  %% PHI / REVERSE-PHI: 20 instances (seeds 3001..3020)
  ppd=25; fr=48; spd=19.2; nf=36; out=fullfile(ROOT,'ReversePhi',tag,'v25'); dov=N/ppd;
  for Ty=[1 -1], for rep=1:20
    e=fullfile(out,sprintf('RP_grain016_type%+d_rep%02d',Ty,rep)); if exist(fullfile(e,'frames.mat'),'file'), continue; end
    mI=make_ReversePhi_stimulus_1D_grain(N,nf,fr,dov,10000,0,spd,meanL,Ty,0.8,16,3000+rep);
    x=struct('experiment','RP_reps20','grain_px',16,'type',Ty,'rep',rep,'contrast',0.8,'px_per_deg',ppd,'framerate',fr, ...
      'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'speed_dps',spd,'speed_px_per_frame',spd*ppd/fr,'gt_direction_deg',(Ty<0)*180);
    dnn_v25_save(e,dnn_v25_broadcast(mI,NR),x);
  end, end

  %% INTERACTION AMONG SCALES, CSr: 25 pairs of random phases per envelope width
  c=dnn_v25_config('DH',N); dov=N/c.px_per_deg; out=fullfile(ROOT,'DH',tag,'v25');
  for s=c.sigma_t_ms, for rep=1:25
    e=fullfile(out,sprintf('DH_CSr_sig%05.1f_rep%02d',s,rep)); if exist(fullfile(e,'frames.mat'),'file'), continue; end
    rng(5000+rep); pm=rand*2*pi; pst=rand*2*pi;
    mI=make_DH_stimulus_1D_v25(N,c.n_frames,c.framerate,dov,4,s,pm,pst,c.tf_hz,c.sfm,c.sfs,meanL);
    x=struct('experiment','DH_reps25','condition','CSr','sigma_t_ms',s,'rep',rep,'phase_mov',pm,'phase_stat',pst, ...
      'px_per_deg',c.px_per_deg,'framerate',c.framerate,'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'tf_hz',c.tf_hz);
    dnn_v25_save(e,dnn_v25_broadcast(mI,NR),x);
  end, end

  %% BREAKING OF TIME-REVERSAL SYMMETRY: 40 phases (seeds 7001..7040)
  ppd=128/1.5; fr=40; nf=36; vd=10*fr/ppd; dov=N/ppd; out=fullfile(ROOT,'TimeReversal',tag,'v25');
  STIM=struct('name',{'LightEdges','UnevenSquareWave'},'period_deg',{1.17,2.34});
  for k=1:2, for ph=1:40
    rng(7000+ph); phase=rand*2*pi;
    S=make_TimeReversal_stimulus_1D(STIM(k).name,N,nf,fr,dov,phase,1,50,vd,STIM(k).period_deg);
    for rev=[0 1]
      if rev, M=flipud(S); suf='reversed'; else, M=S; suf='normal'; end
      e=fullfile(out,sprintf('TR_%s_%s_ph%02d',STIM(k).name,suf,ph)); if exist(fullfile(e,'frames.mat'),'file'), continue; end
      x=struct('experiment','TR_reps40','stim_name',STIM(k).name,'reversed',rev,'phase',phase,'phase_idx',ph,'v_dps',vd, ...
        'period_deg',STIM(k).period_deg,'px_per_deg',ppd,'framerate',fr,'dov',dov,'N_cols',N,'N_rows',NR,'meanL',50,'contrast',1, ...
        'speed_px_per_frame',vd*ppd/fr,'gt_direction_deg',rev*180);
      dnn_v25_save(e,dnn_v25_broadcast(M,NR),x);
    end
  end, end

  %% SECOND-ORDER MOTION: 100 noise carriers (seeds 4001..4100)
  c=dnn_v25_config('secondorder',N); dov=N/c.px_per_deg; out=fullfile(ROOT,'SecondOrder',tag,'v25');
  for Ty=[1 -1], for rep=1:100
    e=fullfile(out,sprintf('SO_grain016_type%+d_rep%03d',Ty,rep)); if exist(fullfile(e,'frames.mat'),'file'), continue; end
    rng(4000+rep); pm=rand*2*pi; pc=rand*2*pi;
    mI=make_SecondOrder_stimulus_1D_grain(N,c.n_frames,c.framerate,dov,c.sigma_t_ms,pm,pc,c.speed_dps,1.0,4.0,meanL,Ty,1,16);
    x=struct('experiment','SO_reps100','carrier_grain_px',16,'type',Ty,'rep',rep,'phase_mov',pm,'phase_carrier',pc, ...
      'px_per_deg',c.px_per_deg,'framerate',c.framerate,'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL, ...
      'speed_dps',c.speed_dps,'speed_px_per_frame',c.speed_dps*c.px_per_deg/c.framerate,'gt_direction_deg',(Ty<0)*180);
    dnn_v25_save(e,dnn_v25_broadcast(mI,NR),x);
  end, end

  %% MASKING: 50 random phases per mask frequency (1 cpd signal, 11 ratios)
  c=dnn_v25_config('Masking',N); fr=c.framerate; dov=N/c.px_per_deg; out=fullfile(ROOT,'Masking',tag,'v25');
  for r=[0 0.02 0.05 0.2 0.25 0.5 1 3 4 6 12], for rep=1:50
    if r==0, msk=0; sfm=1; nm=sprintf('Masking_nomask_rep%02d',rep); else, msk=1; sfm=r; nm=sprintf('Masking_maskr%05.2f_rep%02d',r,rep); end
    e=fullfile(out,nm); if exist(fullfile(e,'frames.mat'),'file'), continue; end
    rng(6000+rep); pm=rand*2*pi;
    mI=make_Masking_stimulus_1D(N,c.n_frames,fr,dov,c.sigma_t_ms,pm,c.speed_dps,1,sfm,meanL,msk,c.contrast_signal,c.contrast_mask,c.tf_mask_hz);
    x=struct('experiment','Masking_reps50','mask_ratio',r,'masking',msk,'sfm_cpd',1,'sfmask_cpd',sfm,'rep',rep,'phase_mov',pm, ...
      'px_per_deg',c.px_per_deg,'framerate',fr,'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'speed_dps',c.speed_dps, ...
      'speed_px_per_frame',c.speed_dps*c.px_per_deg/fr,'gt_direction_deg',0);
    dnn_v25_save(e,dnn_v25_broadcast(mI,NR),x);
  end, end
  fprintf('%s generado\n', tag);
end
disp('GENERACION TERMINADA');
