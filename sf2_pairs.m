function pairs = sf2_pairs(xy, edges)
%SF2_PAIRS Enumerate all unordered observation pairs within separation bins.
%
%   pairs = sf2_pairs(xy, edges) returns rows [i, j, bin], with i < j.
%   Coordinates are planar x + 1i*y. Bins are left-closed/right-open,
%   except that the final bin includes its upper edge.
%   Coincident observations are included when zero is within a bin.
%   This helper is used internally by sf2_estimate; pairs cannot currently
%   be supplied to the estimator as an option.

validateattributes(xy, {'numeric'}, {'vector', 'finite', 'nonempty'});
validateattributes(edges, {'numeric'}, ...
    {'vector', 'real', 'nonnegative', 'finite'});
xy = full(double(xy(:)));
edges = full(double(edges(:)));
if numel(edges) < 2 || any(diff(edges) <= 0)
    error('Edges must increase.');
end

% Enumerate one endpoint at a time to avoid an n-by-n distance matrix.
n = numel(xy);
chunks = cell(n, 1);
for i = 1:n - 1
    j = (i + 1:n)';
    h = abs(xy(j) - xy(i));
    bin = discretize(h, edges);
    keep = ~isnan(bin);
    chunks{i} = [repmat(i, sum(keep), 1), j(keep), bin(keep)];
end
pairs = vertcat(chunks{:});
if isempty(pairs)
    pairs = zeros(0, 3);
end
end
