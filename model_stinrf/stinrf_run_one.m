function stinrf_run_one(script_name)
% STINRF_RUN_ONE  Run one experiment in a workspace of its own.
%
%   Every experiment begins with clear, so running them straight from a
%   driver wipes the driver's own variables along with theirs. Calling them
%   from inside a function keeps that clear confined to this function's
%   workspace, so the loop that called it survives.
%
% Inputs
%   script_name (char) name of the experiment script, without .m
%
% Outputs
%   nothing; the script saves and exports as it normally would

    eval(script_name);
end
