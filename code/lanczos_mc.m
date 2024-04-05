%% clean
clear
close all
clc

%% import utilities
addpath("utils")

%% define number of Lanczos iterations, number of samples and averages, width of CI
ks = [10; 50];
N = 1000;
avgs = 5;
alpha = 0.05;

%% load matrix and compute exact quantity
datastruct = load("../matrices/nos3.mat");
M = datastruct.Problem.A;
G = ichol(M);
diaginvM = diag(inv(M));
avg_errors = zeros(size(ks, 1), avgs, N);

%% compute Lanczos estimator for a fixed value of k
for j = 1:size(ks, 1)
    [Ts, Vs] = lanczos_iterations(M, G, ks(j));
    [est, W] = compute_lanczos_estimator(G, Ts, Vs);
    for l = 1:avgs
        ests = compute_mc_estimator(M, N, W);
        ests = ests + est; % add Lanczos estimate
        errors = vecnorm(ests-repmat(diaginvM, 1, N)) / norm(diaginvM);
        avg_errors(j, l, :) = errors;
    end
end

mean_errors = mean(avg_errors, 2);
std_dev = std(avg_errors, [], 2) / sqrt(avgs);
cdi = norminv(1-alpha/2);

%% plot and save figure
fig = figure();
fig_legend_string = [];
for j = 1:size(ks, 1)
    curve1 = mean_errors(j, :) + cdi * std_dev(j, :) / sqrt(avgs);
    curve2 = mean_errors(j, :) - cdi * std_dev(j, :) / sqrt(avgs);
    x = (1:N);
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    pl = loglog(x, mean_errors(j, :), 'LineWidth', 4);
    hold on
    kstr = sprintf("$k = %.0f$", string(ks(j)));
    fig_legend_string = [fig_legend_string, kstr];
    c = get(pl,'Color');
    h = fill(x2, inBetween, c);
    set(h, 'facealpha', .2);
    fig_legend_string = [fig_legend_string, ""];
    hold on
end

xlabel("$N$", 'interpreter', 'latex', 'FontSize', 15);
ylabel("$\frac{\vert \vert \mathbf{d}_{\mathrm{Lanczos-MC}}^{k, N} - diag(A^{-1}) \vert \vert_2}{\vert \vert diag(A^{-1}) \vert \vert_2}$", ...
    'interpreter', 'latex', 'FontSize', 18);
title("Lanczos Monte Carlo estimator", 'FontSize', 15);
legend(fig_legend_string, 'interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
saveas(fig, "../figures/lanczos_mc_estimator", "epsc");