out = inme_setup();
cfg = stINRF3_config('Cont'); model = cfg.model;
N=400; dov=4.0; framerate=120; duration_ms=500; speed_dps=8; meanL=100;
sfm_env=1; sfm_car=4; Carrier_type=1; Type=[1 -1]; sigma_t=10000;
n_trials=100; sigma_c=57.74;
n_frames=round(duration_ms/1000*framerate);
j0=round(N/2); j_sensors=(j0-99):(j0+100);
R=zeros(numel(Type),n_trials);
for si=1:numel(Type)
    r=zeros(1,n_trials); ty=Type(si);
    parfor tr=1:n_trials
        pm=rand*2*pi; pc=rand*2*pi;
        mI=make_SecondOrder_stimulus_1D(N,n_frames,framerate,dov,sigma_t,pm,pc, ...
              speed_dps,sfm_env,sfm_car,meanL,ty,Carrier_type);
        O=stINRF3_compute_1D(mI,dov,model,framerate,N,j_sensors);
        r(tr)=stINRF3_readout_gaussian_1D(O,j_sensors,sigma_c);
    end
    R(si,:)=r;
end
writematrix(R',fullfile(out,'stinrf_secondorder_trials.csv'));
fprintf('derecha: media %+.4f sd %.4f  positivos %d/%d\n', mean(R(1,:)), std(R(1,:)), sum(R(1,:)>0), n_trials);
fprintf('izquierda: media %+.4f sd %.4f  negativos %d/%d\n', mean(R(2,:)), std(R(2,:)), sum(R(2,:)<0), n_trials);
disp('written');
