# SF2 Gaussian reference toolbox

MATLAB R2016b or newer; no additional toolboxes. Add this directory to your MATLAB path, then run `demo_sf2` or `test_sf2`.

```matlab
model = struct('type', 'exp', 'ell', 20);
S = sf2_estimate(xy, X, 0:5:100, model);
```

`xy = x + 1i*y` contains planar coordinates, not longitude/latitude. Distance units must match `ell` and bin edges. `X` is a real scalar field. All unordered pairs within the bins receive equal weights. Observations with NaN values are removed before pair enumeration, preserving original endpoint indices in the outputs; coordinates must be finite. Explicit weighting is not supported. Pre-binned pair lists can be supplied as described below.

## Reuse pre-binned pairs

```matlab
pp = sf2_pairs(xy, edges); % M-by-3 rows [i, j, bin].
S = sf2_estimate(xy, X, edges, model, 'Pairs', pp);
```

Edges may be omitted when supplying pairs:

```matlab
S = sf2_estimate(xy, X, [], model, 'Pairs', pp);
```

In this case labels in `pp(:,3)` define the bins, `S.edges` is empty, and the number of bins is `max(pp(:,3))`, determined before NaN-pair filtering. Missing label numbers yield empty bins; trailing empty bins cannot be inferred. Empty supplied pairs with empty edges produce zero bins. Mean separations are still calculated from actual retained pairs. Bin assignments cannot be checked against edges when no edges are provided; supplied labels need not correspond to contiguous separation intervals. If edges are supplied, the existing consistency check remains active.

Indices refer to the original coordinate/observation vectors. The function checks integer indices, `1 <= i < j <= numel(X)`, valid bin numbers, unique unordered pairs, and (when edges are supplied) consistency of bin labels with the current coordinates and edges. Pairs involving NaN observations are omitted after validation; original indices are preserved in the output. Do not renumber indices to skip missing observations.

Reuse `pp` across data vectors or model shapes on the same coordinate layout. Subsets and reordered rows are allowed, defining an equally weighted estimator over exactly the supplied pairs. Explicit `[]` or `zeros(0,3)` returns an empty estimator; only omission of `Pairs` triggers automatic enumeration. The option works with both trace and direct covariance methods and does not apply to `sf2_estimate_1D`.

This skips O(n^2) pair enumeration. Validation still computes distances for the M supplied pairs and checks uniqueness, currently O(M log M); covariance is recomputed for each call. No geometry or covariance cache persists between function calls.

Run `test_sf2_prebinned` to compare automatic and supplied pairs, including reordered rows, missing observations, subsets, empty arrays, and invalid-input rejection.

## Regular 1D records

```matlab
x = linspace(0, 1000, 2000)';
model = struct('type', 'exp', 'ell', 5);
S = sf2_estimate_1D(x, X, edges, model, ...
    'CrossCovariance', true, 'ConfidenceLevel', 0.95);
```

`sf2_estimate_1D` accepts the same four positional inputs as `sf2_estimate`, except the first input is a real, strictly increasing or decreasing regular grid with at least two points. Row and column vectors are accepted. NaNs in X are omitted pairwise at their original grid locations. Linear and logarithmic separation bins are supported, including multiple integer lags within a bin.

The functions return the same principal scientific fields, including `r`, `D2`, `D2shape`, `count`, `fractionalCovariance`, `covariance`, `relativeVariance`, `variance`, `stderr`, `effectiveDof`, `confidenceInterval`, and `confidenceLevel`. Model amplitudes are not required. The 1D output additionally contains `spacing`, `lags`, `lagBin`, and `lagCount`; `method` is `'regular'`.

**Diagnostic difference:** `pairs` and `pairBin` are empty by default in the 1D function, to avoid potentially quadratic pair storage. Request `'ReturnPairs', true` if needed. Pair row ordering can differ; indices still refer to the original input record. Trace-specific `Method`, `CacheMB`, and `BlockSize` options do not apply to the 1D function.

For integer lags p and q, and start offset h = start_q - start_p, define the shape increment covariance

G_pq(h) = [f(|h-p| dx) + f(|h+q| dx) - f(|h+q-p| dx) - f(|h| dx)] / 2.

The code sums 2*G_pq(h)^2 times the number of retained start pairs having that offset, divides by the bin pair counts, and accumulates all lag combinations in the bins. Both offset signs are included. Complete records use exact analytic overlap counts. Records with missing values use integer counts from zero-padded FFT correlations of valid-start masks. This is an exact finite-record sum under the same Gaussian closure, without a periodic-record assumption or covariance truncation.

For N observations and J retained integer lags, complete-record covariance evaluation takes O(N J^2) time and O(N) working storage, plus lag metadata and the bin covariance output. With missing values, correlations can increase runtime to O(N J^2 log N). Wide bins or a large maximum separation may contain many lags, so the method is not universally faster than the trace backend. `'CrossCovariance', false` skips lag combinations belonging to different bins. No dense N-by-N matrices are formed. Explicit diagnostic pairs add O(N J) storage when requested.

Regularity is checked within a tolerance of `1e-8*spacing + 8*eps(max(abs(x)))`; accepted coordinate roundoff is treated as an exact grid with spacing set by the endpoints. In contrast, the arbitrary-geometry function uses actual pair distances, which can place nominally edge-adjacent lags in different bins due to coordinate roundoff. Avoid placing bin edges exactly on integer lags when comparing floating-point coordinate grids; half-lag boundaries are convenient for single-lag bins.

Run `demo_sf2_1D` for a 2000-point exponential-covariance example and `test_sf2_1D` for comparisons with the full function. The demo uses a stationary AR(1) recursion to generate the sampled exponential covariance without a dense Cholesky factorization.

## Shape parameters only

| Model | SF2 shape f(h) | Required parameters |
|---|---|---|
| `exp` | 2[1 - exp(-h/ell)] | ell > 0 |
| `gaussian` | 2[1 - exp(-(h/ell)^2)] | ell > 0 |
| `power` | h^zeta | 0 < zeta < 2 |

`'scale-free'` and `'scalefree'` remain aliases for `'power'`. No `sigma2`, `D0`, or `r0` is needed. The finite-scale normalizations correspond to unit field variance, purely for convenience. The power-law shape is not a stationary field covariance; it describes stationary increments without an outer scale. The Gaussian convention has no factor of 1/2 in its exponent. The one-dimensional spectral relation q = zeta + 1 should not be interpreted as a two-dimensional spectral-density exponent.

## Fractional covariance and plug-in absolute uncertainty

Write the population SF2 as D(h) = a f(h), with unknown common amplitude a. For pair u = (i,j) and pair v = (p,q), calculate the shape increment covariance

G_uv = [f(|j-p|) + f(|i-q|) - f(|j-q|) - f(|i-p|)] / 2.

For bins k and l, the shape-based covariance is

Vshape_kl = 2 sum_(u in k, v in l) G_uv^2 / (N_k N_l).

Let fbar_k be the average of f over the actual pair separations in bin k. The fractional covariance is

F_kl = Vshape_kl / (fbar_k fbar_l).

The unknown amplitude cancels exactly. Shared endpoints and dependence between disjoint pairs are retained. All four endpoint distances enter the model, including distances outside the SF bins. Bin means use the actual separations, not the model at the plotted mean separation.

The absolute covariance still requires SF2 amplitudes. The toolbox uses the observed bin estimates:

Vhat_kl = F_kl Dhat_k Dhat_l.

This is a **plug-in Gaussian reference estimate**, not a known-model variance or an unbiased variance estimator. It does not constrain the observed SF2 to a common fitted amplitude. Fractional covariance depends only on geometry and model shape; plug-in absolute covariance also depends on the realization. Zero observed SF2 can yield zero plug-in uncertainty even when the population variance is positive. Under a fully specified synthetic model, the exact reference covariance is recovered externally as `S.fractionalCovariance .* (Dpopulation * Dpopulation.')`, where `Dpopulation` contains the model bin means.

## Outputs

| Field | Meaning |
|---|---|
| `D2` | Observed SF2 in each bin |
| `D2shape` | Bin-averaged model shape, with arbitrary normalization |
| `r` | Mean pair separation, for plotting |
| `fractionalCovariance` | Amplitude-independent Gaussian reference covariance F |
| `relativeVariance` | Diagonal of F |
| `effectiveDof` | 2 / relativeVariance; moment-matched chi-square degrees of freedom |
| `covariance` | Plug-in absolute covariance F .* (D2 * D2.') |
| `variance`, `stderr` | Diagonal of covariance and its square root |
| `confidenceInterval` | Two columns: lower and upper pointwise chi-squared confidence limits |
| `confidenceLevel` | Requested coverage, default 0.95 |
| `count` | Number of contributing pairs |
| `pairs`, `pairBin` | Retained pair indices and bin assignments, for diagnostics |

`D2model` was replaced by `D2shape` because an absolute model prediction cannot be supplied without amplitude. Weight-related outputs were removed. All per-bin vectors are columns. Bins are left-closed/right-open, except the final upper edge is included. Empty bins have zero counts and NaN estimates and uncertainty. A bin containing only zero-separation pairs has undefined fractional covariance and therefore NaN uncertainty.

## Approximate chi-squared confidence intervals

```matlab
S = sf2_estimate(xy, X, edges, model, 'ConfidenceLevel', 0.95);
lower = S.confidenceInterval(:, 1);
upper = S.confidenceInterval(:, 2);
```

The approximation is Dhat / Dpopulation ~ chi2(nu) / nu, where nu = 2 / relativeVariance. For coverage 1-alpha, inversion gives

[nu * Dhat / chi2_quantile(1-alpha/2, nu),
 nu * Dhat / chi2_quantile(alpha/2, nu)].

These are intervals for the population SF2 averaged over the contributing pair separations, not prediction intervals for another realization. They are pointwise, not simultaneous across separation bins. Cross-bin correlation is not used to construct these marginal intervals.

The degrees of freedom depend on model shape and sampling geometry. The interval is generally a moment-matched approximation for dependent increments, even under Gaussian closure; it is exact for an equally weighted mean of independent zero-mean Gaussian squares with equal variances. Model-parameter uncertainty and non-Gaussian fourth-order contributions are excluded. The asymmetric intervals need not coincide with estimate +/- a multiple of the plug-in standard error.

Quantiles use MATLAB's base `gammaincinv`, so no Statistics Toolbox is required. Empty bins and bins with undefined fractional variance return NaN limits. An observed zero SF2 with valid degrees of freedom returns [0,0] from the inverted pivot; this boundary result should not be interpreted as evidence of zero population uncertainty for quantized or otherwise imperfect data.

## Assumptions and limitations

The reference calculation assumes homogeneous, isotropic second-order increments with zero ensemble mean, and jointly Gaussian increments for the fourth-order closure. Deterministic spatial trends violate the zero-mean assumption: squared-increment covariance then also contains 4 mu_u mu_v G_uv. Subtracting an estimated trend changes the sampling covariance and is not automatically accounted for.

Model shape parameters are treated as specified. Parameter-fitting uncertainty, non-Gaussian fourth-order contributions, and independent observation errors are excluded. Duplicate coordinates are treated as observations of the same modeled field value. For complex velocity, analyze individual real components; adding component variances omits their cross covariance. This is a spatial model: do not flatten observations from different times into a single spatial sample without justifying that model.

Integer and single-precision inputs are converted to full double precision before arithmetic. Coordinates must be finite even where observations are missing.

Covariance sums are internally rescaled to avoid overflow or underflow from squaring an arbitrary model amplitude. Material violations of the Gaussian fractional-variance bounds trigger warning `sf2_estimate:NumericalCovariance` and NaN uncertainties for affected bins, rather than a falsely collapsed confidence interval. The bound check does not guarantee numerical accuracy: very small increments between widely separated clusters can still suffer cancellation in four-point differences. Compare the graph and direct methods for demanding geometries; both can lose precision at extreme scale ratios.

The finite-geometry sums use no large-domain approximation. The Gaussian closure is an assumption; conversion to absolute uncertainty is a plug-in approximation.

## Computation: Appendix B graph formulation

The default algorithm uses the observation graph, with one normalized Laplacian per separation bin:

L_a = B_a' B_a / N_a,

where B_a is the signed pair-incidence matrix. With S_ij = f(distance(i,j)), the shape covariance is

Vshape_ab = trace(S L_a S L_b) / 2.

This is Appendix B's trace identity, extended to cross-bin covariance and normalized by actual bin-averaged shapes. It is exactly equivalent to the pair sum above. The code builds sparse Laplacians directly from endpoint indices; it never forms B_a or a pair-by-pair covariance matrix.

```matlab
S = sf2_estimate(xy, X, edges, model); % Default: Method = 'trace'.
S = sf2_estimate(xy, X, edges, model, 'CrossCovariance', false);
```

The graph calculation removes observations that contribute no retained pair from its internal matrices, while diagnostic pair indices still refer to the original data. It builds a dense n-by-n model-shape matrix and sparse Laplacians, where n is the number of participating observations.

Graph products T_a = S L_a are cached in groups, using a default cache budget of 256 MiB. Within a group, every product is computed once. Each product for a later bin is streamed once and compared with every member of the group. If all products fit in the cache, each is computed only once in total, instead of repeatedly for preceding bins.

```matlab
S = sf2_estimate(xy, X, edges, model, 'CacheMB', 512);
```

`CacheMB` controls cached products only, not total process memory. The model-shape matrix, one streamed product, pair arrays, sparse Laplacians, and arithmetic temporaries require additional memory. Each cached double-precision product occupies approximately 8*n^2 bytes. At n = 2000, one product occupies 30.5 MiB: the default cache holds eight products. `CacheMB = 0` retains one product, the minimum needed, and reproduces the earlier low-memory traversal. The option is unused by `Method = 'direct'`.

Computing T_a costs O(n nnz(L_a)); each trace is evaluated as `sum(sum(T_a .* T_b.'))`, at O(n^2) cost, without dense matrix multiplication. With K occupied bins, the old traversal computed K(K+1)/2 products. A cache that holds all K requires only K products; a smaller cache incurs some recomputation. All cross-bin traces must still be evaluated, so overall runtime need not improve by the same factor as the sparse-product count. The trace method generally benefits when retained pairs greatly outnumber observations, but still requires dense observation-level storage. Pair enumeration itself remains O(n_valid^2), where n_valid is the number of nonmissing observations.

The original pair-block implementation is retained as an independent verification option:

```matlab
Scheck = sf2_estimate(xy, X, edges, model, ...
    'Method', 'direct', 'BlockSize', 256);
```

For `Method = 'direct'`, temporary covariance memory scales as `BlockSize^2`, and runtime remains O(M^2). `BlockSize` is unused by the default trace method. No thinning or covariance truncation occurs in either method. `CrossCovariance = false` returns NaN off-diagonal entries; the trace method then skips cross-bin products and traces.

## Validation

`test_sf2` checks both exact algorithms for all three models against an independent incidence-matrix calculation, Monte Carlo population means and covariance, changes of cache capacity, block size, and data amplitude, power-law coordinate-unit invariance, missing observations, empty bins, the single-pair identity, chi-squared quantiles and confidence-level ordering, interval scaling, numeric input classes, missing-data index mapping, unequal bins and edge inclusion, tiny shape amplitudes, and Monte Carlo coverage in the exact single-pair case. Neither MATLAB nor Octave was available in the authoring environment, so these MATLAB tests remain unrun. Independent Python numerical checks validate the covariance algebra and normalization, but do not establish MATLAB runtime correctness.

Reference: Isserlis, L. (1918), On a formula for the product-moment coefficient of any order of a normal frequency distribution in any number of variables, Biometrika, 12, 134–139. https://doi.org/10.1093/biomet/12.1-2.134
