function dnn_v25_save(edir, tensor, extra)
% DNN_V25_SAVE  Write a stimulus in the format the DNN pipelines expect.
%
%   Writes frames.mat with the variable `tensor` (N_rows x N_cols x
%   n_frames, single) plus every field of `extra` as its own variable, so
%   run_dnn_fullvideo_gauss.py / run_multiframe_pipeline.py can copy the
%   provenance into result.mat.
%
% Inputs
%   edir (char) output directory
%   tensor (single, rows x cols x frames)
%   extra (struct) metadata fields

    if ~exist(edir, 'dir'), mkdir(edir); end
    tensor = single(tensor); %#ok<NASGU>
    fn = fieldnames(extra);
    args = {'tensor'};
    for k = 1:numel(fn)
        eval([fn{k} ' = extra.' fn{k} ';']);
        args{end+1} = fn{k}; %#ok<AGROW>
    end
    save(fullfile(edir, 'frames.mat'), args{:}, '-v7.3');
end
