function out = sf2_estimate(xy, X, edges, model, varargin)
%SF2_ESTIMATE SF2 and Gaussian reference uncertainty for arbitrary 2D samples.
%
%   S = sf2_estimate(xy, X, edges, model)
%   S = sf2_estimate(..., 'CrossCovariance', true)
%   S = sf2_estimate(..., 'Pairs', sf2_pairs(xy, edges))
%   S = sf2_estimate(xy, X, [], model, 'Pairs', pp)
%
% Inputs
%   xy      Planar coordinates x + 1i*y, in any consistent distance units.
%   X       Real scalar observations. NaNs are omitted pairwise.
%   edges   Edge vector or m-by-2 [lower, upper] interval array. Intervals
%           may be unordered, separated by gaps, or overlapping. By default
%           all unordered pairs within each bin are used with equal weights. Bins include their lower edge.
%           Upper edges are excluded except the globally largest upper edge,
%           included in every bin ending there. Overlaps retain all memberships.
%           Use [] with supplied Pairs to infer bins from their labels.
%   model   Shape parameters only:
%             struct('type', 'exp',        'ell', ell)
%             struct('type', 'gaussian',   'ell', ell)
%             struct('type', 'power',      'zeta', zeta)
%           Model amplitude cancels from fractional covariance.
%           The power-law range is 0 < zeta < 2.
%
% Options
%   Pairs            Cell vector from sf2_pairs: pp{b} is an M_b-by-2
%                    [i,j] matrix using original observation indices, i < j.
%                    Skips enumeration; subsets and reordered rows allowed.
%                    Duplicates within a bin are rejected; pairs may appear
%                    in several bins. Membership is checked when edges are
%                    supplied. NaN pairs are omitted within each cell.
%                    With edges=[], numel(pp) defines the bin count, including
%                    trailing empty bins. Explicit {} or [] means no pairs;
%                    omit the option for automatic enumeration.
%   Method           'trace' (default): Appendix B trace formulation.
%                    'direct': direct pair-block sum for verification.
%   BlockSize        Pair-block size for Method = 'direct' (default 512).
%   CacheMB          Budget for cached graph products, in MiB (default 256).
%                    Additional working matrices are required; see README.
%   ConfidenceLevel  Pointwise chi-squared CI coverage (default 0.95).
%   CrossCovariance  Compute cross-bin covariance (default false).
%                    Off-diagonal entries are NaN when false.
%
% Principal outputs
%   D2                  Observed, equally weighted SF2 in each bin.
%   r                   Arithmetic mean of contributing pair separations,
%                       after omitting pairs with missing observations.
%                       NaN for empty bins. Prefer this plotting coordinate
%                       to nominal bin centers for irregular sampling.
%   edges               Input bin edges, describing the full bin extent;
%                       empty when edges were omitted. Without edges, bin
%                       count is numel(pp), preserving trailing empty cells.
%   fractionalCovariance Gaussian reference covariance divided by the
%                       product of model bin means; amplitude independent.
%   relativeVariance    Diagonal of fractionalCovariance.
%   covariance          Plug-in absolute covariance:
%                       fractionalCovariance .* (D2 * D2.').
%   variance, stderr    Diagonal of covariance and its square root.
%   effectiveDof        2 ./ relativeVariance (moment matching only).
%   D2shape             Bin-averaged model SF2 shape.
%   confidenceInterval  nb-by-2 array of pointwise [lower, upper] limits
%                       for population bin-averaged SF2. Moment-matched
%                       chi-squared approximation; not a simultaneous band.
%   confidenceLevel     Requested coverage probability.
%
% Finite-bin interpretation: D2 estimates the population SF2 averaged over
% the actual contributing pair separations, not necessarily SF2 at r.
% In general, mean(D2population(separation)) ~= D2population(mean(separation)).
% For wide bins, compare with a model averaged over the pair separations,
% as done for D2shape, rather than evaluating the model only at r.
%
% Assumptions: isotropic, homogeneous second-order increments with zero
% ensemble mean; jointly Gaussian increments for the fourth-order closure.
% Fractional covariance sums use the actual finite pair geometry. Absolute
% uncertainty substitutes observed SF2 for the unknown population values;
% this plug-in estimate is not an unbiased variance estimator in general.
% Parameter-fitting uncertainty and non-Gaussian contributions are excluded.

%% Validate inputs and options.
validateattributes(xy, {'numeric'}, {'vector', 'finite', 'nonempty'});
validateattributes(X, {'numeric'}, {'vector', 'real'});
xy = full(double(xy(:)));
X = full(double(X(:)));
if numel(X) ~= numel(xy) || any(isinf(X))
    error('X must match xy and contain finite values or NaNs.');
end
[~, binInfo] = sf2_pairs([], edges);
edges = binInfo.edges;
binLimits = binInfo.limits;
includeUpper = binInfo.includeUpper;
clear binInfo
hasEdges = ~isempty(binLimits);
nb = size(binLimits, 1);

ip = inputParser;
addParameter(ip, 'Pairs', {});
addParameter(ip, 'Method', 'trace');
addParameter(ip, 'BlockSize', 512);
addParameter(ip, 'CacheMB', 256);
addParameter(ip, 'ConfidenceLevel', 0.95);
addParameter(ip, 'CrossCovariance', false);
parse(ip, varargin{:});
options = ip.Results;
method = validatestring(options.Method, {'trace', 'direct'});
validateattributes(options.BlockSize, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
validateattributes(options.CrossCovariance, {'logical'}, {'scalar'});
validateattributes(options.CacheMB, {'numeric'}, ...
    {'scalar', 'real', 'nonnegative', 'finite'});
validateattributes(options.ConfidenceLevel, {'numeric'}, ...
    {'scalar', 'real', '>', 0, '<', 1, 'finite'});
options.BlockSize = double(options.BlockSize);
options.CacheMB = double(options.CacheMB);
options.ConfidenceLevel = double(options.ConfidenceLevel);
sf2_model(0, model); % Validate the model even when no pairs are retained.

%% Obtain pairs and calculate the observed and model-shape bin means.
validX = ~isnan(X);
originalIndex = find(validX);
usePairs = ~any(strcmp(ip.UsingDefaults, 'Pairs'));
if ~hasEdges && ~usePairs
    error('sf2_estimate:MissingEdges', ...
        'Nonempty edges are required unless the Pairs option is supplied.');
end
if usePairs
    pp = options.Pairs;
    % Cell position defines the bin label, including empty bins.
    if isnumeric(pp) && isempty(pp)
        pp = cell(0, 1);
    end
    if ~iscell(pp) || (~isempty(pp) && ~isvector(pp))
        error('sf2_estimate:InvalidPairs', 'Pairs must be a vector cell array of index matrices.');
    end
    pp = pp(:);
    if hasEdges
        if isempty(pp)
            pp = repmat({zeros(0, 2)}, nb, 1);
        elseif numel(pp) ~= nb
            error('sf2_estimate:InvalidPairs', 'Pairs must contain one cell per bin.');
        end
    else
        nb = numel(pp);
    end
else
    % Exclude missing observations before enumeration, then restore original
    % indices within each cell. This also handles fewer than two valid data.
    pp = sf2_pairs(xy(originalIndex), edges);
    for b = 1:nb
        pp{b} = reshape(originalIndex(pp{b}), size(pp{b}));
    end
end
% Release references before filtering cells and allocating dense products.
options = rmfield(options, 'Pairs');
clear ip varargin
count = zeros(nb, 1);
D2 = nan(nb, 1);
D2shape = nan(nb, 1);
r = nan(nb, 1);
for b = 1:nb
    p = pp{b};
    if usePairs
        validateattributes(p, {'numeric'}, {'2d', 'real', 'finite', 'integer'});
        if isempty(p) && isequal(size(p), [0, 0])
            p = zeros(0, 2);
        end
        if size(p, 2) ~= 2
            error('sf2_estimate:InvalidPairs', 'Each pair cell must have two columns [i,j].');
        end
        p = full(double(p));
        if any(p(:, 1) < 1 | p(:, 2) > numel(X) | p(:, 1) >= p(:, 2))
            error('sf2_estimate:InvalidPairs', 'Each pair requires 1 <= i < j <= numel(X).');
        end
        if size(unique(p, 'rows'), 1) ~= size(p, 1)
            error('sf2_estimate:DuplicatePairs', 'A pair may occur only once within each bin.');
        end
    end
    separation = abs(xy(p(:, 2)) - xy(p(:, 1)));
    if usePairs && hasEdges
        consistent = separation >= binLimits(b, 1) ...
            & (separation < binLimits(b, 2) ...
            | (includeUpper(b) & separation == binLimits(b, 2)));
        if any(~consistent)
            error('sf2_estimate:InconsistentPairs', 'Supplied pairs do not belong to bin %d.', b);
        end
    end
    % Filter before bin statistics. Validation precedes masking so invalid
    % indices or memberships cannot be hidden by missing observations.
    keep = validX(p(:, 1)) & validX(p(:, 2));
    if ~all(keep)
        p = p(keep, :);
        separation = separation(keep);
    end
    pp{b} = p;
    count(b) = size(p, 1);
    if count(b) > 0
        D2(b) = mean((X(p(:, 2)) - X(p(:, 1))).^2);
        D2shape(b) = mean(sf2_model(separation, model));
        r(b) = mean(separation);
    end
end
clear p separation validX keep consistent

% Normalize the covariance calculation by a representative shape amplitude.
% This changes no fractional result, but avoids squaring extremely small or
% large arbitrary shape amplitudes (especially after changing distance units).
positiveShape = D2shape(isfinite(D2shape) & D2shape > 0);
if isempty(positiveShape)
    shapeScale = 1;
else
    shapeScale = max(positiveShape);
end

%% Evaluate the covariance using either exact formulation.
% The trace method is Appendix B extended to finite bins and cross-bin
% covariance. The direct method independently evaluates the original sum.
if strcmp(method, 'trace')
    shapeCovariance = trace_covariance(xy, pp, count, model, ...
        options.CrossCovariance, options.CacheMB, shapeScale);
else
    shapeCovariance = direct_covariance(xy, pp, count, model, options, shapeScale);
end
clear pp

%% Remove model amplitude and convert to plug-in absolute uncertainty.
% Both numerator and denominator scale as the amplitude squared. Using
% bin-averaged shapes retains finite-bin effects, rather than evaluating
% the model at a representative separation.
normalizedShape = D2shape / shapeScale;
fractionalCovariance = (shapeCovariance ./ normalizedShape) ./ normalizedShape.';
undefined = count == 0 | D2shape == 0;
fractionalCovariance(undefined, :) = NaN;
fractionalCovariance(:, undefined) = NaN;
if ~options.CrossCovariance
    fractionalCovariance(~logical(eye(nb))) = NaN;
end
% A nonzero Gaussian quadratic form has 2/rank <= F <= 2, with rank
% at most min(number of pairs, number of nonmissing observations minus one).
% Flag material numerical violations rather than return false certainty.
for a = find(~undefined).'
    lowerBound = 2 / min(count(a), numel(originalIndex) - 1);
    value = fractionalCovariance(a, a);
    if ~isfinite(value) || value < lowerBound * (1 - 1e-7) || value > 2 * (1 + 1e-7)
        warning('sf2_estimate:NumericalCovariance', ...
            ['Unreliable covariance in bin %d (numerical cancellation or ' ...
             'extreme model scales). Returning NaN uncertainty for this bin.'], a);
        fractionalCovariance(a, :) = NaN;
        fractionalCovariance(:, a) = NaN;
    end
end
covariance = fractionalCovariance .* (D2 * D2.');
relativeVariance = diag(fractionalCovariance);
variance = diag(covariance);
effectiveDof = 2 ./ relativeVariance;

%% Invert the moment-matched chi-squared pivot for population SF2.
% D2hat / D2population is approximated by chi2(nu)/nu, nu = 2/F.
% Inverting this ratio reverses the quantiles: the upper chi2 quantile
% gives the lower confidence limit, and the lower gives the upper limit.
% The chi-squared approximation is generally not exact for dependent
% increments. These are pointwise intervals, conditional on model shape;
% model-parameter uncertainty is not included.
confidenceInterval = nan(nb, 2);
valid = isfinite(D2) & isfinite(effectiveDof) & effectiveDof > 0;
nu = effectiveDof(valid);
tail = (1 - options.ConfidenceLevel) / 2;
% chi2inv(p, nu) = 2*gammaincinv(p, nu/2).
if any(valid)
    qLower = 2 * gammaincinv(tail * ones(size(nu)), nu / 2);
    qUpper = 2 * gammaincinv(tail * ones(size(nu)), nu / 2, 'upper');
    confidenceInterval(valid, 1) = D2(valid) .* (nu ./ qUpper);
    confidenceInterval(valid, 2) = D2(valid) .* (nu ./ qLower);
end
% Nonzero model increments cannot have zero fractional variance for a finite
% sample. Numerical failures have already been replaced by NaNs above.

out = struct('edges', edges, 'r', r, 'D2', D2, 'D2shape', D2shape, ...
    'fractionalCovariance', fractionalCovariance, ...
    'relativeVariance', relativeVariance, 'effectiveDof', effectiveDof, ...
    'covariance', covariance, 'variance', variance, 'stderr', sqrt(variance), ...
    'confidenceInterval', confidenceInterval, ...
    'confidenceLevel', options.ConfidenceLevel, ...
    'count', count, 'model', model, 'method', method);
end

function V = trace_covariance(xy, pp, count, model, crossCovariance, cacheMB, shapeScale)
%TRACE_COVARIANCE Observation-level trace calculation from Appendix B.
% For each bin a, L_a = B_a' * B_a / N_a is its normalized Laplacian.
% With S(i,j) = f(distance(i,j))/shapeScale,
% V(a,b) = trace(S * L_a * S * L_b) / 2.
nb = numel(count);
V = zeros(nb, nb);
if ~any(count)
    return
end

% Mark participating observations without concatenating or sorting all pairs.
used = false(numel(xy), 1);
for a = 1:nb
    used(pp{a}(:)) = true;
end
map = zeros(numel(xy), 1);
map(used) = (1:sum(used))';
xy = xy(used);
n = numel(xy);
S = sf2_model(abs(xy - xy.'), model) / shapeScale;
L = cell(nb, 1);
for a = 1:nb
    if count(a) == 0
        continue
    end
    p = reshape(map(pp{a}), size(pp{a}));
    i = p(:, 1);
    j = p(:, 2);
    L{a} = sparse([i; j; i; j], [i; j; j; i], ...
        [ones(2 * count(a), 1); -ones(2 * count(a), 1)], n, n) / count(a);
end

active = find(count > 0);
if ~crossCovariance
    % Compute diagonal variances one bin at a time.
    for a = active.'
        Ta = S * L{a};
        V(a, a) = 0.5 * sum(sum(Ta .* Ta.'));
    end
    return
end

% Cache a group of T_a = S*L_a matrices. Compare every member of that
% group with each later T_b, computing the latter only once per group.
% If the entire set fits, each sparse product is calculated just once.
% The budget covers cached products only: S, one streamed product, and
% trace/product temporaries require additional O(n^2) working memory.
% Even CacheMB = 0 retains one product, the minimum needed by this method.
bytesPerProduct = 8 * double(n)^2;
cacheCount = min(numel(active), ...
    max(1, floor(cacheMB * 2^20 / bytesPerProduct)));
for first = 1:cacheCount:numel(active)
    last = min(first + cacheCount - 1, numel(active));
    group = active(first:last);
    cached = cell(numel(group), 1);
    for u = 1:numel(group)
        cached{u} = S * L{group(u)};
    end

    % Diagonal and cross-bin traces within the cached group.
    for u = 1:numel(group)
        a = group(u);
        Ta = cached{u};
        V(a, a) = 0.5 * sum(sum(Ta .* Ta.'));
        for v = u + 1:numel(group)
            b = group(v);
            V(a, b) = 0.5 * sum(sum(Ta .* cached{v}.'));
            V(b, a) = V(a, b);
        end
    end

    % Cross-bin traces against later groups. Reuse each streamed product
    % for all cached bins; never form the dense matrix product Ta*Tb.
    for next = last + 1:numel(active)
        b = active(next);
        Tb = S * L{b};
        for u = 1:numel(group)
            a = group(u);
            V(a, b) = 0.5 * sum(sum(cached{u} .* Tb.'));
            V(b, a) = V(a, b);
        end
    end
    % Release references before the next group allocates dense products.
    clear cached Ta Tb
end
end

function V = direct_covariance(xy, pp, count, model, options, shapeScale)
%DIRECT_COVARIANCE Direct Gaussian covariance sum, retained for verification.
% Cov(delta_u^2, delta_v^2) = 2*G_uv^2. All four endpoint distances enter
% G_uv, including distances outside the SF bins. Blocks limit memory only.
nb = numel(count);
% Only this backend needs a flat list for its existing pair-block algorithm.
ij = vertcat(pp{:});
if isempty(ij)
    ij = zeros(0, 2);
end
bin = zeros(sum(count), 1);
first = 1;
for b = 1:nb
    if count(b) > 0
        bin(first:first + count(b) - 1) = b;
    end
    first = first + count(b);
end
np = size(ij, 1);
pairWeight = 1 ./ count(bin);
V = zeros(nb, nb);
blockSize = options.BlockSize;

for ia = 1:blockSize:np
    u = ia:min(ia + blockSize - 1, np);
    zu1 = xy(ij(u, 1));
    zu2 = xy(ij(u, 2));
    A = sparse((1:numel(u))', bin(u), pairWeight(u), numel(u), nb);

    % Process only the upper triangle of blocks, then add its transpose.
    for ja = ia:blockSize:np
        v = ja:min(ja + blockSize - 1, np);
        zv1 = xy(ij(v, 1)).'; % Nonconjugating transpose for coordinates.
        zv2 = xy(ij(v, 2)).';
        G = 0.5 * (sf2_model(abs(zu2 - zv1), model) / shapeScale ...
            + sf2_model(abs(zu1 - zv2), model) / shapeScale ...
            - sf2_model(abs(zu2 - zv2), model) / shapeScale ...
            - sf2_model(abs(zu1 - zv1), model) / shapeScale);
        T = 2 * G.^2;
        if ~options.CrossCovariance
            T = T .* (bin(u) == bin(v).');
        end
        B = sparse((1:numel(v))', bin(v), pairWeight(v), numel(v), nb);
        Q = full(A' * T * B);
        if ia == ja
            V = V + Q;
        else
            V = V + Q + Q';
        end
    end
end

end
