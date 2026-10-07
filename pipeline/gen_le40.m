% Light edges, configuracion definitiva (LE_v2): 170.67 ppd, 160 fps, 9.375 deg/s -> 10 px/frame, 40 fases
out = inme_setup();
G = struct('tag',{'raft','pwc_lfn2','meattn_dual'},'N_cols',{960,1024,768},'N_rows',{520,436,768});
ppd=512/3; fr=160; nf=36; vd=9.375; per=1.17; meanL=50;
for g=1:numel(G)
  N=G(g).N_cols; NR=G(g).N_rows; dov=N/ppd; out=fullfile(fullfile(out,'stimuli','LightEdges_v2'),G(g).tag,'v25');
  for ph=1:40
    rng(7000+ph); phase=rand*2*pi;
    S=make_TimeReversal_stimulus_1D('LightEdges',N,nf,fr,dov,phase,1,meanL,vd,per);
    for rev=[0 1]
      if rev, M=flipud(S); suf='reversed'; else, M=S; suf='normal'; end
      x=struct('experiment','TimeReversal_LE_v2_reps40','stim_name','LightEdges','reversed',rev,'phase',phase,'phase_idx',ph, ...
        'v_dps',vd,'period_deg',per,'temporal_freq_hz',vd/per,'px_per_deg',ppd,'framerate',fr,'dov',dov,'N_cols',N,'N_rows',NR, ...
        'meanL',meanL,'contrast',1,'speed_px_per_frame',vd*ppd/fr,'gt_direction_deg',rev*180);
      dnn_v25_save(fullfile(out,sprintf('TR_LightEdges_%s_ph%02d',suf,ph)),dnn_v25_broadcast(M,NR),x);
    end
  end
  fprintf('%s LE40\n',G(g).tag);
end
disp('LE40 GENERADO');
