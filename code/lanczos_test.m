%% clean
clear
close all
clc

%% import utilities
addpath("utils")

%% define number of Lanczos iterations
k = 10;

%% load matrix and compute exact quantity
datastruct = load("../matrices/nos3.mat");
M = datastruct.Problem.A;
G = ichol(M);
diaginvM = diag(inv(M));
errors = zeros(k, 1);

%% compute Lanczos estimator for every value of k
for j=1:k
    [V, T] = lanczos(M, G, rand(size(M, 1), 1), j);
    L = chol(T);
    W = inv(G)' * V * inv(L)';
    est = vecnorm(W, 2, 2);
    error = vecnorm(est-diaginvM) / norm(diaginvM);
    errors(j, :) = error;
end


%% plot and save figure
fig = figure();
x = (1:k);
loglog(x, errors, 'LineWidth', 4);
hold on
xlabel("$k$", 'interpreter', 'latex', 'FontSize', 15);
ylabel("$\frac{\vert \vert \mathbf{d}_{\mathrm{Lanczos}}^k- diag(A^{-1}) \vert \vert_2}{\vert \vert diag(A^{-1}) \vert \vert_2}$", ...
    'interpreter', 'latex', 'FontSize', 18);
title("Lanczos estimator", 'FontSize', 15);
legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
saveas(fig, "../figures/lanczos_estimator", "epsc");