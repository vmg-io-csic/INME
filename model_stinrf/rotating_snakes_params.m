function p = rotating_snakes_params(variant)
%ROTATING_SNAKES_PARAMS  Geometry and luminance profile of a snakes wheel.
%
%   The wheel is built from n_rings concentric annuli, each divided into
%   n_units identical angular units. Every unit carries the same luminance
%   profile, and successive rings are rotated with respect to one another by
%   angular_offset units, which is what gives the wheel its spiral appearance.
%
%   Inputs
%     variant  char. Which wheel to build:
%              'illusion'  the asymmetric wheel that elicits illusory rotation.
%                          Each unit runs through four luminance levels, mid
%                          grey, black, dark grey and white, joined by smooth
%                          transitions, and the rings are mutually offset.
%              'symmetric' the control wheel. Each unit is a square wave,
%                          half white and half black, and the rings are in
%                          phase, so the luminance gradient that drives the
%                          illusion is absent.
%
%   Outputs
%     p  struct with fields:
%        edges           1 x (n_rings+1) double. Ring boundaries, as a
%                        fraction of the wheel radius, from centre to rim.
%        n_rings         double. Number of concentric annuli.
%        n_units         double. Angular units per revolution.
%        angular_offset  double. Rotation between consecutive rings, in units.
%        unit_profile    1 x 30 double. Luminance of one angular unit, sampled
%                        at 30 equally spaced phases, in [0, 100].
%
%   Note on n_units: it must not be a multiple of four. The model's Gabor bank
%   has only horizontal and vertical orientation channels, so a wheel that is
%   invariant under a 90 degree rotation makes the tangential projection cancel
%   to numerical zero over the whole disc. With 23 units that degeneracy does
%   not arise.
%
%   Example
%     p = rotating_snakes_params('illusion');

if nargin < 1, variant = 'illusion'; end

% Ring boundaries, shared by both variants.
p.edges = [             0  0.08695652174   0.1471571906    0.220735786   0.3110367893 ...
             0.4180602007   0.5384615385   0.6789297659   0.8361204013              1];
p.n_rings = 9;
p.n_units = 23;

switch lower(variant)
    case 'illusion'
        p.angular_offset = -0.320864096305;
        % Four luminance levels per unit -- mid grey, black, dark grey, white --
        % with smooth transitions between them.
        p.unit_profile = [   56.89891434    57.13027718    57.13289349    57.11285105    57.36836195    53.09469536 ...
                             26.95269652    3.151241866   0.1006372631  0.08723703887   0.1365204776  0.09833822506 ...
                             0.890828723     12.5139196    32.60789255    37.90363001    37.82898631    37.88901713 ...
                             37.85238984    37.42327141    42.84181487    72.70612639    97.21675242    99.78295225 ...
                             99.72367774    99.75256457    99.70774939    99.06532648    86.12564064    63.44470085];

    case 'symmetric'
        p.angular_offset = 0;
        % Square wave: the first half of the unit white, the second half black.
        p.unit_profile = [100*ones(1,15) zeros(1,15)];

    otherwise
        error('rotating_snakes_params:variant', ...
              'Unknown variant "%s". Use ''illusion'' or ''symmetric''.', variant);
end
end
