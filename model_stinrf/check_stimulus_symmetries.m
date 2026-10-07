function [d, names] = check_stimulus_symmetries(S)
% CHECK_STIMULUS_SYMMETRIES  Test the eight space/contrast/time reversal symmetries.
%
%   Applies the operators chi (space reversal), Gamma (contrast reversal)
%   and Theta (time reversal) and their combinations, and reports for each
%   the smallest mean absolute difference from the original allowing any
%   cyclic shift in space and time. A value of zero means the stimulus
%   possesses that symmetry.
%
%   IMPORTANT: S must cover an integer number of spatial periods and an
%   integer number of temporal cycles, otherwise the cyclic shifts cannot
%   express the phase offset and every symmetry will appear broken.
%
% Inputs
%   S (double, T x N) space-time profile, ideally in the +/-1 convention
%
% Outputs
%   d (1 x 7 double) minimum discrepancy for each operator combination
%   names (1 x 7 cell) operator labels, in the order of d

    S = S - mean(S(:));                 % contrast reversal acts on the profile
    [T, N] = size(S);
    combos = [0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
    names  = {'Theta', 'Gamma', 'GammaTheta', 'chi', 'chiTheta', ...
              'chiGamma', 'chiGammaTheta'};

    d = zeros(1, size(combos,1));
    for c = 1:size(combos,1)
        Tr = S;
        if combos(c,3), Tr = circshift(flipud(Tr), 1, 1); end   % Theta
        if combos(c,1), Tr = circshift(fliplr(Tr), 1, 2); end   % chi
        if combos(c,2), Tr = -Tr;                        end   % Gamma

        best = Inf;
        for dt = 0:T-1
            Tt = circshift(Tr, dt, 1);
            for dx = 0:N-1
                v = mean(abs(S - circshift(Tt, dx, 2)), 'all');
                if v < best, best = v; end
            end
        end
        d(c) = best;
    end
end
