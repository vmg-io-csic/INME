function out = inme_output()
%INME_OUTPUT  Folder where every script in this repository writes its results.
%
%   It is <repository>/output by default, or the value of the INME_OUT
%   environment variable when that variable is set, which is the way to send
%   results to a different disk without editing any script. The folder is
%   created if it does not exist.
%
%   Outputs
%     out  char. Absolute path of the output folder.
%
%   Example
%     writematrix(R, fullfile(inme_output(), 'masking.csv'));

out = getenv('INME_OUT');
if isempty(out)
    out = fullfile(fileparts(mfilename('fullpath')), 'output');
end
if ~exist(out, 'dir')
    mkdir(out);
end
end
