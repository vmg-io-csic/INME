% Uneven square wave at the physical speed used for the model: 1.17 deg/s,
% period 2.34 deg, 85.33 px/deg, 9.984 fps, i.e. 10 px per frame; 15 frames
% (1.5 s) and 40 phases (the same seeds).
out = inme_setup();
G = struct('tag',{'raft','pwc_lfn2','meattn_dual'},'N_cols',{960,1024,768},'N_rows',{520,436,768});
ppd=128/1.5; vd=1.17; per=2.34; fr=vd*ppd/10; nf=15; meanL=50;
for g=1:numel(G)
  N=G(g).N_cols; NR=G(g).N_rows; dov=N/ppd; out=fullfile(fullfile(out,'stimuli','Uneven_v2'),G(g).tag,'v25');
  for ph=1:40
    rng(7000+ph); phase=rand*2*pi;
    S=make_TimeReversal_stimulus_1D('UnevenSquareWave',N,nf,fr,dov,phase,1,meanL,vd,per);
    for rev=[0 1]
      if rev, M=flipud(S); suf='reversed'; else, M=S; suf='normal'; end
      x=struct('experiment','TimeReversal_Uneven_v2_reps40','stim_name','UnevenSquareWave','reversed',rev,'phase',phase,'phase_idx',ph, ...
        'v_dps',vd,'period_deg',per,'temporal_freq_hz',vd/per,'px_per_deg',ppd,'framerate',fr,'n_frames',nf,'duration_ms',1000*nf/fr, ...
        'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'contrast',1,'speed_px_per_frame',vd*ppd/fr,'gt_direction_deg',rev*180);
      dnn_v25_save(fullfile(out,sprintf('TR_UnevenSquareWave_%s_ph%02d',suf,ph)),dnn_v25_broadcast(M,NR),x);
    end
  end
  fprintf('%s uneven15\n',G(g).tag);
end
disp('UNEVEN15 GENERADO');
