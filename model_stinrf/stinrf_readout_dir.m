function d = stinrf_readout_dir(O_H, O_V, mode, half_region)
% STINRF_READOUT_DIR  Direction read out of the model's two output maps.
%
%   Shared by the gratings and plaid experiments, which used to carry a copy
%   of this switch each. One copy means the three cannot end up reading their
%   output in slightly different ways.
%
% Inputs
%   O_H, O_V (double) horizontal and vertical model outputs
%   mode (char) 'single', 'multi' or 'gauss'
%   half_region (double) half-width of the readout support [px]
%
% Outputs
%   d (double) direction [deg]

    switch mode
        case 'single'
            [d, ~] = stINRF3_readout_single(O_H, O_V);
        case 'multi'
            [d, ~, ~] = stINRF3_readout_multisensor(O_H, O_V, half_region);
        case 'gauss'
            sigma_c = 2*half_region / sqrt(12);
            [d, ~, ~] = stINRF3_readout_gaussian(O_H, O_V, sigma_c, ...
                            'SupportHalf', half_region);
        otherwise
            error('stinrf_readout_dir: unknown mode ''%s''.', mode);
    end
end
