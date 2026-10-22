function [pairs, binInfo] = sf2_pairs(xy, edges)
%SF2_PAIRS Enumerate unordered observation pairs and define separation bins.
%
%   pairs = sf2_pairs(xy, edges) returns one cell per bin. pairs{b} is
%   an M_b-by-2 [i,j] matrix, with i < j. Empty bins contain zeros(0,2).
%   [~, binInfo] = sf2_pairs([], edges) validates the bins and returns their
%   normalized edges, limits, endpoint flags, and vectorEdges flag without
%   enumerating pairs. Empty xy returns empty cells, one per bin; empty
%   xy and empty edges return a zero-bin cell array.
%   Edge vectors use discretize; independent intervals use membership checks.
%   Coordinates are planar x + 1i*y. Bins are left-closed/right-open,
%   except that the final bin includes its upper edge for vector input.
%   edges may instead be m-by-2 [lower, upper] intervals, with arbitrary
%   order, gaps, and overlaps. Every upper edge is excluded except the
%   globally largest upper edge, which is included in all bins ending there.
%   Overlaps repeat [i,j] across cells, never within a cell.
%   Coincident observations are included when zero is within a bin.
%   Reuse the output with sf2_estimate(..., 'Pairs', pairs).
%   Pair indices refer to the original coordinate vector, including any
%   observations whose data values may subsequently be omitted as NaNs.

validateattributes(xy, {'numeric'}, {'2d', 'finite'});
if ~isempty(xy)
    validateattributes(xy, {'numeric'}, {'vector'});
end
xy = full(double(xy(:)));
[edges, limits, includeUpper, vectorEdges] = sf2_bin_spec(edges);
binInfo = struct('edges', edges, 'limits', limits, ...
    'includeUpper', includeUpper, 'vectorEdges', vectorEdges);
nb = size(limits, 1);
pairs = repmat({zeros(0, 2)}, nb, 1);
if isempty(xy)
    return
end
if nb == 0
    error('sf2:InvalidBins', 'Nonempty bins are required for pair enumeration.');
end

% Append to per-bin buffers with geometric capacity growth. This avoids a
% global [i,j,bin] table and repeated concatenation of each bin's full data.
n = numel(xy);
count = zeros(nb, 1);
capacity = zeros(nb, 1);
for i = 1:n - 1
    j = (i + 1:n)';
    h = abs(xy(j) - xy(i));
    rowsByBin = bin_membership(h, edges, limits, includeUpper, vectorEdges);
    for b = find(~cellfun('isempty', rowsByBin)).'
        selected = j(rowsByBin{b});
        needed = count(b) + numel(selected);
        if needed > capacity(b)
            capacity(b) = max([64, needed, 2 * capacity(b)]);
            pairs{b}(capacity(b), 2) = 0;
        end
        rows = count(b) + (1:numel(selected));
        pairs{b}(rows, :) = [repmat(i, numel(selected), 1), selected];
        count(b) = needed;
    end
end
for b = 1:nb
    pairs{b} = pairs{b}(1:count(b), :);
end
end

function [edges, limits, includeUpper, vectorEdges] = sf2_bin_spec(edges)
%SF2_BIN_SPEC Validate edge vectors or independent [lower, upper] intervals.
% Vector bins retain discretize's endpoint convention. For independent
% intervals, only bins reaching the globally largest upper edge are closed
% on the right. This is independent of row order and reproduces vector bins.
validateattributes(edges, {'numeric'}, {'2d', 'real', 'nonnegative', 'finite'});
edges = full(double(edges));
vectorEdges = isvector(edges);
if isempty(edges)
    edges = zeros(0, 1);
    limits = zeros(0, 2);
    includeUpper = false(0, 1);
    return
end
if vectorEdges
    edges = edges(:);
    if numel(edges) < 2 || any(diff(edges) <= 0)
        error('sf2:InvalidBins', 'Vector edges must contain at least two increasing values.');
    end
    limits = [edges(1:end - 1), edges(2:end)];
elseif size(edges, 2) == 2
    limits = edges;
    if any(limits(:, 1) >= limits(:, 2))
        error('sf2:InvalidBins', 'Each interval must satisfy lower < upper.');
    end
else
    error('sf2:InvalidBins', 'Bins must be an edge vector or an m-by-2 interval array.');
end
includeUpper = limits(:, 2) == max(limits(:, 2));
end

function rowsByBin = bin_membership(h, edges, limits, includeUpper, vectorEdges)
%BIN_MEMBERSHIP Group separation indices by bin without a flat pair table.
nb = size(limits, 1);
rowsByBin = cell(nb, 1);
if vectorEdges
    labels = discretize(h, edges);
    rows = find(~isnan(labels));
    if isempty(rows)
        return
    end
    % Group discretize's labels once, rather than scan h once per bin.
    [sortedLabels, order] = sort(labels(rows));
    rows = rows(order);
    first = [1; find(diff(sortedLabels) ~= 0) + 1];
    last = [first(2:end) - 1; numel(rows)];
    for k = 1:numel(first)
        rowsByBin{sortedLabels(first(k))} = rows(first(k):last(k));
    end
else
    for b = 1:nb
        keep = h >= limits(b, 1) & (h < limits(b, 2) ...
            | (includeUpper(b) & h == limits(b, 2)));
        rowsByBin{b} = find(keep);
    end
end
end
