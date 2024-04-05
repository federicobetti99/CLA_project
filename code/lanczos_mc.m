%% clean
clear
close all
clc

%% define number of Lanczos iterations, number of samples and averages, width of CI
ks = [10; 50];
N = 1000;
avgs = 10;
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
        avg_errors(j, l, :) = vecnorm(ests-repmat(diaginvM, 1, N)) / norm(diaginvM);
    end
end

mean_errors = mean(avg_errors, 2);
std_dev = std(avg_errors, [], 2) / sqrt(avgs);
cdi = norminv(1-alpha/2);
x = (1:N);
x2 = [x, fliplr(x)];

%% plot and save figure
fig = figure();
fig_legend_string = [];
for j = 1:size(ks, 1)
    curve1 = mean_errors(j, :) + cdi * std_dev(j, :) / sqrt(avgs);
    curve2 = mean_errors(j, :) - cdi * std_dev(j, :) / sqrt(avgs);
    inBetween = [curve1, fliplr(curve2)];
    pl = loglog(x, mean_errors(j, :), 'LineWidth', 3);
    hold on
    kstr = sprintf("$k = %.0f$", string(ks(j)));
    fig_legend_string = [fig_legend_string, kstr];
    h = fill(x2, inBetween, get(pl,'Color'));
    set(h, 'facealpha', .2);
    fig_legend_string = [fig_legend_string, ""];
    hold on
end

xlabel("$N$", 'interpreter', 'latex', 'FontSize', 15);
ylabel("$\frac{\| \mathbf{d}_{\mathrm{Lanczos-MC}}^{k, N} - \mathrm{diag}(A^{-1}) \|_2}{\| \mathrm{diag}(A^{-1}) \|_2}$", ...
    'interpreter', 'latex', 'FontSize', 18);
title("Lanczos Monte Carlo estimator", 'FontSize', 15);
legend(fig_legend_string, 'interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
saveas(fig, "../figures/lanczos_mc_estimator", "epsc");