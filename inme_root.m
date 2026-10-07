function root = inme_root()
%INME_ROOT  Absolute path of the repository root.
%
%   Resolved from the location of this file, so it is independent of the current
%   working directory and of where the repository has been cloned.
%
%   Outputs
%     root  char. Absolute path of the folder containing this file.

root = fileparts(mfilename('fullpath'));
end
