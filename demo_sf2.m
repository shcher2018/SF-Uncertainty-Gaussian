% DEMO_SF2 Irregular 2D sampling with amplitude-free uncertainty models.
rng(7);

%% Generate a synthetic exponential-covariance field.
N = 1e3;
% xy = 100 * (rand(N, 1) + 1i * rand(N, 1)); % Planar coordinates in km.
xy = 100*rand(N,1);
model = struct('type', 'exp', 'ell', 20);
sigma2True = 1; % Needed to simulate the field, not to estimate uncertainty.
C = sigma2True * exp(-abs(xy - xy.') / model.ell);
X = chol(C + 1e-12 * eye(numel(xy)), 'lower') * randn(numel(xy), 1);

%% Estimate SF2 and plug-in Gaussian reference uncertainty.
edges = 0:5:55;
S = sf2_estimate(xy, X, edges, model);

figure
loglog(S.r,S.D2,'k-','linewidth',2)
hold on
grid on
errorbar(S.r, S.D2, S.stderr, 'bo');
hold on;
errorbar(S.r*1.02, S.D2, S.D2 - S.confidenceInterval(:, 1), ...
    S.confidenceInterval(:, 2) - S.D2, 'ro');

% Only a simulation supplies a known population amplitude for comparison.
plot(S.r, sigma2True * S.D2shape, 'k-');
xlabel('Separation (km)');
ylabel('SF2');
legend('Estimate +/- 1 SE', ...
    '95% pointwise CI',...
    'Known bin-averaged SF2', 'Location', 'best');
