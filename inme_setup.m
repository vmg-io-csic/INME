function out = inme_setup()
%INME_SETUP  Put the repository on the MATLAB path and return the output folder.
%
%   Call this once at the start of any script in this repository. It locates the
%   repository from the position of this file, so it works from any working
%   directory and on any machine, and it adds the repository root and the folders
%   holding the model, the stimulus generators and the pipeline scripts to the
%   MATLAB path.
%
%   Outputs
%     out  char. Absolute path of the folder where results are written. It is
%          <repository>/output by default, or the value of the INME_OUT
%          environment variable if that variable is set. The folder is created
%          if it does not exist.
%
%   Example
%     out = inme_setup();
%     writematrix(R, fullfile(out, 'masking.csv'));

root = fileparts(mfilename('fullpath'));
addpath(root);                       % so inme_output is always reachable
addpath(fullfile(root, 'model_stinrf'));
addpath(fullfile(root, 'stimuli_dnn'));
addpath(fullfile(root, 'pipeline'));

out = inme_output();
end
