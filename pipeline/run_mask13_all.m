ps=parallel.Settings; ps.Pool.AutoCreate=false;

run('/tmp/claude-1000/-home-mberraquero/c2798504-67d8-481a-a009-5fc07ff6d143/scratchpad/run_mask_sweep13.m');
out = inme_setup();
run('run_Masking_1D.m'); disp('FIG MASKING TERMINADA');

run('run_Masking_nomask_1D.m'); disp('FIG NOMASK TERMINADA');
