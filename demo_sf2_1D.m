% DEMO_SF2_1D Regular 1D sampling with exponential covariance.
rng(7);
N = 2000;
x = linspace(0, 1000, N)';
model = struct('type', 'exp', 'ell', 20);
sigma2True = 1;

% An AR(1) recursion generates the exact sampled exponential covariance.
% Initialize from its stationary distribution; no spinup or dense Cholesky.
phi = exp(-(x(2) - x(1)) / model.ell);
innovationSD = sqrt(sigma2True * (-expm1(2 * log(phi))));
X = zeros(n, 1);
X(1) = sqrt(sigma2True) * randn;
for j = 2:n
    X(j) = phi * X(j - 1) + innovationSD * randn;
end

edges = (0:5:55)';
tic;
S = sf2_estimate_1D(x, X, edges, model);
fprintf('Regular 1D calculation: %.3f s\n', toc);

figure;
errorbar(S.r, S.D2, S.D2 - S.confidenceInterval(:, 1), ...
    S.confidenceInterval(:, 2) - S.D2, 'o');
hold on;
plot(S.r, sigma2True * S.D2shape, 'k-');
xlabel('Separation');
ylabel('SF2');
legend('Estimate with approximate 95% pointwise CI', ...
    'Known bin-averaged SF2', 'Location', 'best');

% Logarithmic bins are also allowed. Narrow bins with no integer lag return NaN.
% edges = logspace(log10(0.5), log10(55), 16);
% S = sf2_estimate_1D(x, X, edges, model, 'CrossCovariance', false);
% X(100:110) = NaN; % Gaps remain on the original grid; no compression occurs.
