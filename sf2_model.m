function shape = sf2_model(h, model)
%SF2_MODEL Isotropic SF2 shape; its arbitrary amplitude cancels in uncertainty.
%
%   shape = sf2_model(h, model)
%
%   exp:        shape(h) = 2 * [1 - exp(-h/ell)]
%   gaussian:   shape(h) = 2 * [1 - exp(-(h/ell)^2)]
%   power:      shape(h) = h^zeta, with 0 < zeta < 2
%
% The finite-scale shapes use unit field variance for normalization. The
% power-law normalization depends on the coordinate units. Common
% multiplicative factors cancel from fractional covariance.
% 'scale-free' and 'scalefree' remain aliases for 'power'.
% Model parameters are type, ell, and zeta.

validateattributes(h, {'numeric'}, {'real', 'nonnegative', 'finite'});
h = full(double(h));
if ~isstruct(model) || ~isscalar(model) || ~isfield(model, 'type')
    error('Model must be a scalar structure with a type field.');
end

switch lower(char(model.type))
    case {'exp', 'gaussian'}
        validateattributes(model.ell, {'numeric'}, ...
            {'scalar', 'real', 'positive', 'finite'});
        t = h / double(model.ell);
        if strcmpi(model.type, 'gaussian')
            t = t.^2;
        end
        % expm1 retains accuracy at separations much smaller than ell.
        shape = 2 * (-expm1(-t));

    case {'power', 'scale-free', 'scalefree'}
        validateattributes(model.zeta, {'numeric'}, ...
            {'scalar', 'real', '>', 0, '<', 2, 'finite'});
        shape = h.^double(model.zeta);

    otherwise
        error('Unknown model type. Use exp, gaussian, or power.');
end
end
