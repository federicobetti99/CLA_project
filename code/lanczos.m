%% clean
clear
close all
clc

%% load matrix and compute exact quantity, define number of Lanczos iterations
datastruct = load("../matrices/nos3.mat");
M = datastruct.Problem.A;
G = ichol(M);
diaginvM = diag(inv(M));
k = 500;
errors = zeros(k, 1);

%% compute Lanczos estimator for every value of j = 1, ..., k
[Ts, Vs] = lanczos_iterations(M, G, k);
for j = 1:k
    [est, W] = compute_lanczos_estimator(G, Ts(1:j, 1:j), Vs(:, 1:j));
    errors(j, :) = vecnorm(est-diaginvM) / norm(diaginvM);
end

%% plot and save figure
fig = figure();
x = (1:k);
loglog(x, errors, 'LineWidth', 3);
hold on
xlabel("$k$", 'interpreter', 'latex', 'FontSize', 15);
ylabel("$\frac{\| \mathbf{d}_{\mathrm{Lanczos}}^k- \mathrm{diag}(A^{-1}) \|_2}{\| \mathrm{diag}(A^{-1}) \|_2}$", ...
    'interpreter', 'latex', 'FontSize', 18);
title("Lanczos estimator", 'FontSize', 15);
saveas(fig, "../figures/lanczos_estimator", "epsc");