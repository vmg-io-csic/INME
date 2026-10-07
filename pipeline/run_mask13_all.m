ps=parallel.Settings; ps.Pool.AutoCreate=false;

out = inme_setup();
run(fullfile(inme_root(), 'pipeline', 'run_mask_sweep13.m'));
run('run_Masking_1D.m'); disp('MASKING FIGURE DONE');

run('run_Masking_nomask_1D.m'); disp('NO-MASK FIGURE DONE');
