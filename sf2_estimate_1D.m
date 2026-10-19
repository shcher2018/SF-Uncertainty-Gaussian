function out = sf2_estimate_1D(x, X, edges, model, varargin)
%SF2_ESTIMATE_1D SF2 and Gaussian uncertainty for a regular 1D record.
%
%   S = sf2_estimate_1D(x, X, edges, model)
%   S = sf2_estimate_1D(..., 'CrossCovariance', false, ...
%       'ConfidenceLevel', 0.95, 'ReturnPairs', false)
%
% Inputs match sf2_estimate, except x must be a real, strictly monotone,
% uniformly spaced coordinate vector. Units match edges and model.ell.
% X is a real scalar record; NaNs are omitted pairwise without closing gaps.
% Model shapes: exp/gaussian with ell, or power with 0 < zeta < 2.
% All pairs in each bin receive equal weight. The final upper edge is included.
%
% Outputs match sf2_estimate: r, D2, D2shape, count, fractionalCovariance,
% relativeVariance, effectiveDof, covariance, variance, stderr,
% confidenceInterval, confidenceLevel, edges, model, and method.
% r is the mean separation of actual retained pairs, not the bin center.
% D2 estimates population SF2 averaged over those separations; a nonlinear
% model must likewise be averaged over pairs rather than evaluated only at r.
%
% Options
%   CrossCovariance  Full cross-bin covariance (default true). Set false
%                    for diagonal variances and pointwise intervals only.
%   ConfidenceLevel Pointwise chi-squared CI coverage (default 0.95).
%   ReturnPairs      Materialize diagnostic pairs/pairBin (default false).
%                    Otherwise these fields are empty, avoiding O(N*J)
%                    pair storage. Pair ordering may differ from sf2_estimate.
%   The graph-specific Method, CacheMB, and BlockSize options are not used.
%
% Additional outputs: spacing, lags, lagBin, and lagCount describe the grid
% and contributing integer lags. method is 'regular'.
%
% Algorithm: for every pair of integer lags p,q, sum the squared increment
% covariance over starting-point offsets h, weighted by the exact number
% of valid start pairs. Complete records use analytic overlap counts; NaNs
% use zero-padded FFT cross-correlations of the valid-start masks.
% No periodic wraparound, covariance truncation, or large-N approximation.
% For J retained lags, cost is O(N*J^2) for complete records and up to
% O(N*J^2*log(N)) with missing values. Working memory is O(N), plus O(J)
% lag metadata and O(number of bins^2) output; ReturnPairs adds O(N*J).
% Wide bins can contain many lags: this is not universally faster than graph.
%
% Assumptions and uncertainty interpretation are the same as sf2_estimate:
% zero-mean homogeneous Gaussian increments for the reference covariance;
% model shape fixed; absolute covariance uses the plug-in observed SF2.
% Chi-squared intervals are moment-matched, pointwise approximations.

%% Validate the grid, observations, bins, and common options.
validateattributes(x, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(X, {'numeric'}, {'vector', 'real'});
x = full(double(x(:)));
X = full(double(X(:)));
n = numel(x);
if numel(X) ~= n || any(isinf(X))
    error('X must match x and contain finite values or NaNs.');
end
if n < 2
    error('A regular sampling grid requires at least two coordinates.');
end
delta = diff(x);
if ~(all(delta > 0) || all(delta < 0))
    error('x must be strictly monotone.');
end
% Use endpoints to avoid accumulating local roundoff in the grid spacing.
dx = abs((x(end) - x(1)) / (n - 1));
validateattributes(dx, {'numeric'}, {'scalar', 'positive', 'finite'});
tolerance = 1e-8 * dx + 8 * eps(max(abs(x)));
if any(abs(abs(delta) - dx) > tolerance)
    error('x must be uniformly spaced; use sf2_estimate for irregular sampling.');
end
validateattributes(edges, {'numeric'}, ...
    {'vector', 'real', 'nonnegative', 'finite'});
edges = full(double(edges(:)));
nb = numel(edges) - 1;
if nb < 1 || any(diff(edges) <= 0)
    error('Edges must increase.');
end
ip = inputParser;
addParameter(ip, 'CrossCovariance', true);
addParameter(ip, 'ConfidenceLevel', 0.95);
addParameter(ip, 'ReturnPairs', false);
parse(ip, varargin{:});
options = ip.Results;
validateattributes(options.CrossCovariance, {'logical'}, {'scalar'});
validateattributes(options.ReturnPairs, {'logical'}, {'scalar'});
validateattributes(options.ConfidenceLevel, {'numeric'}, ...
    {'scalar', 'real', '>', 0, '<', 1, 'finite'});
options.ConfidenceLevel = double(options.ConfidenceLevel);
sf2_model(0, model);

%% Aggregate observations by integer lag, then by separation bin.
lags = (1:n - 1)';
lagBin = discretize(lags * dx, edges);
keep = ~isnan(lagBin);
lags = lags(keep);
lagBin = lagBin(keep);
validX = ~isnan(X);
lagCount = zeros(size(lags));
lagSum = zeros(size(lags));
pairChunks = cell(0, 1);
if options.ReturnPairs
    pairChunks = cell(numel(lags), 1);
end
for a = 1:numel(lags)
    p = lags(a);
    starts = (1:n - p)';
    good = validX(starts) & validX(starts + p);
    starts = starts(good);
    lagCount(a) = numel(starts);
    lagSum(a) = sum((X(starts + p) - X(starts)).^2);
    if options.ReturnPairs
        pairChunks{a} = [starts, starts + p, repmat(lagBin(a), numel(starts), 1)];
    end
end
% Zero-count lags contribute nothing, including when all data are missing.
keep = lagCount > 0;
lags = lags(keep);
lagBin = lagBin(keep);
lagCount = lagCount(keep);
lagSum = lagSum(keep);
if options.ReturnPairs && ~isempty(pairChunks)
    pp = vertcat(pairChunks{:});
else
    pp = zeros(0, 3);
end
if isempty(pp)
    pp = zeros(0, 3);
end
pairs = pp(:, 1:2);
pairBin = pp(:, 3);
clear pp pairChunks
count = zeros(nb, 1);
D2 = nan(nb, 1);
D2shape = nan(nb, 1);
r = nan(nb, 1);
if ~isempty(lags)
    count = accumarray(lagBin, lagCount, [nb, 1], @sum, 0);
    D2 = accumarray(lagBin, lagSum, [nb, 1], @sum, 0) ./ count;
    D2shape = accumarray(lagBin, lagCount .* sf2_model(lags * dx, model), ...
        [nb, 1], @sum, 0) ./ count;
    r = accumarray(lagBin, lagCount .* (lags * dx), [nb, 1], @sum, 0) ./ count;
end

%% Exact stationary offset sum for every pair of contributing lags.
positiveShape = D2shape(isfinite(D2shape) & D2shape > 0);
shapeScale = 1;
if ~isempty(positiveShape)
    shapeScale = max(positiveShape);
end
V = zeros(nb, nb);
if ~isempty(lags)
    % Every endpoint distance in this calculation is an integer grid lag
    % between 0 and N-1, so evaluate the model once and use table lookups.
    f = sf2_model((0:n - 1)' * dx, model) / shapeScale;
    complete = all(validX);
    if ~complete
        nfft = 2^nextpow2(2 * n - 1);
    end
    for a = 1:numel(lags)
        p = lags(a);
        np = n - p;
        ba = lagBin(a);
        if ~complete
            maskP = double(validX(1:np) & validX(1 + p:n));
            fftP = fft(maskP, nfft);
        end
        for b = a:numel(lags)
            bb = lagBin(b);
            if ~options.CrossCovariance && ba ~= bb
                continue
            end
            q = lags(b);
            nq = n - q;
            % h = start_q - start_p; retain both positive and negative h.
            h = (-(np - 1):nq - 1)';
            if complete
                multiplicity = max(0, min(np, nq - h) - max(1, 1 - h) + 1);
            else
                maskQ = double(validX(1:nq) & validX(1 + q:n));
                c = real(ifft(conj(fftP) .* fft(maskQ, nfft)));
                % Counts are integers; rounding removes FFT roundoff only.
                multiplicity = max(0, round(c(mod(h, nfft) + 1)));
            end
            G = 0.5 * (f(abs(h - p) + 1) + f(abs(h + q) + 1) ...
                - f(abs(h + q - p) + 1) - f(abs(h) + 1));
            value = 2 * sum(multiplicity .* G.^2) / count(ba) / count(bb);
            if a == b
                V(ba, bb) = V(ba, bb) + value;
            else
                % Reverse lag ordering has the same total covariance.
                % If both lags occupy one bin, both terms enter its variance.
                V(ba, bb) = V(ba, bb) + value;
                V(bb, ba) = V(bb, ba) + value;
            end
        end
    end
end

%% Normalize and report uncertainty using the same conventions as 2D.
normalizedShape = D2shape / shapeScale;
F = (V ./ normalizedShape) ./ normalizedShape.';
undefined = count == 0 | D2shape == 0;
F(undefined, :) = NaN;
F(:, undefined) = NaN;
if ~options.CrossCovariance
    F(~logical(eye(nb))) = NaN;
end
for a = find(~undefined).'
    lowerBound = 2 / min(count(a), sum(validX) - 1);
    value = F(a, a);
    if ~isfinite(value) || value < lowerBound * (1 - 1e-7) || value > 2 * (1 + 1e-7)
        warning('sf2_estimate_1D:NumericalCovariance', ...
            'Unreliable covariance in bin %d. Returning NaN uncertainty.', a);
        F(a, :) = NaN;
        F(:, a) = NaN;
    end
end
relativeVariance = diag(F);
covariance = F .* (D2 * D2.');
variance = diag(covariance);
effectiveDof = 2 ./ relativeVariance;
confidenceInterval = nan(nb, 2);
valid = isfinite(D2) & isfinite(effectiveDof) & effectiveDof > 0;
nu = effectiveDof(valid);
tail = (1 - options.ConfidenceLevel) / 2;
if any(valid)
    qLower = 2 * gammaincinv(tail * ones(size(nu)), nu / 2);
    qUpper = 2 * gammaincinv(tail * ones(size(nu)), nu / 2, 'upper');
    confidenceInterval(valid, 1) = D2(valid) .* (nu ./ qUpper);
    confidenceInterval(valid, 2) = D2(valid) .* (nu ./ qLower);
end
out = struct('edges', edges, 'r', r, 'D2', D2, 'D2shape', D2shape, ...
    'fractionalCovariance', F, 'relativeVariance', relativeVariance, ...
    'effectiveDof', effectiveDof, 'covariance', covariance, ...
    'variance', variance, 'stderr', sqrt(variance), ...
    'confidenceInterval', confidenceInterval, ...
    'confidenceLevel', options.ConfidenceLevel, 'count', count, ...
    'pairs', pairs, 'pairBin', pairBin, 'model', model, 'method', 'regular', ...
    'spacing', dx, 'lags', lags, 'lagBin', lagBin, 'lagCount', lagCount);
end
