% Second-order motion, 100 noise carriers, at the sampling of the configuration
% table: 100 px/deg, 40 fps, 8 deg/s, i.e. 20 px per frame.
out = inme_setup();
G = struct('tag',{'raft','pwc_lfn2','meattn_dual'},'N_cols',{960,1024,768},'N_rows',{520,436,768});
ppd=100; fr=40; spd=8; nf=36; sig=10000; meanL=100;
for g=1:numel(G)
  N=G(g).N_cols; NR=G(g).N_rows; dov=N/ppd; out=fullfile(fullfile(out,'stimuli','SecondOrder'),G(g).tag,'v25');
  for Ty=[1 -1], for rep=1:100
    rng(4000+rep); pm=rand*2*pi; pc=rand*2*pi;
    mI=make_SecondOrder_stimulus_1D_grain(N,nf,fr,dov,sig,pm,pc,spd,1.0,4.0,meanL,Ty,1,16);
    x=struct('experiment','SO_reps100','carrier_grain_px',16,'type',Ty,'rep',rep,'phase_mov',pm,'phase_carrier',pc, ...
      'px_per_deg',ppd,'framerate',fr,'dov',dov,'N_cols',N,'N_rows',NR,'meanL',meanL,'speed_dps',spd, ...
      'speed_px_per_frame',spd*ppd/fr,'gt_direction_deg',(Ty<0)*180);
    dnn_v25_save(fullfile(out,sprintf('SO_grain016_type%+d_rep%03d',Ty,rep)),dnn_v25_broadcast(mI,NR),x);
  end, end
  fprintf('%s SO 40fps\n', G(g).tag);
end
disp('SO REGENERADO');
