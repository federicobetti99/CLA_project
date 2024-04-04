%% clean
clear
close all
clc

%% import utilities and fix seed
addpath("utils")
rng("default")
rng(0)

%% define number of samples, number of averages
N = 10000;
avgs = 10;

%% load matrix and compute exact quantity
datastruct = load("../matrices/nos3.mat");
M = datastruct.Problem.A;
diaginvM = diag(inv(M));
avg_errors = zeros(avgs, N);

%% compute Monte Carlo estimator, averaging for every value of N
for j=1:avgs
    fprintf("Carrying out average %d out of %d \n", j, avgs);
    ests = mc_estimator(M, N);
    errors = vecnorm(ests-repmat(diaginvM, 1, N)) / norm(diaginvM);
    avg_errors(j, :) = errors;
end

mean_errors = mean(avg_errors, 1);
std_dev = std(avg_errors, [], 1) / sqrt(avgs);
curve1 = mean_errors + std_dev;
curve2 = mean_errors - std_dev;
x = (1:N);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

%% plot and save figure
fig = figure();
fig_legend_string = ["$\mathcal{O}(1/\sqrt{N})$", ""];
loglog(x, 10 ./ sqrt(x), 'LineWidth', 2);
hold on
loglog(x, mean_errors, 'LineWidth', 4);
hold on
h = fill(x2, inBetween, 'b');
set(h, 'facealpha', .2);
hold on
xlabel("$N$", 'interpreter', 'latex', 'FontSize', 15);
ylabel("$\frac{\vert \vert \mathbf{d}_{\mathrm{MC}}^N - diag(A^{-1}) \vert \vert_2}{\vert \vert diag(A^{-1}) \vert \vert_2}$", ...
    'interpreter', 'latex', 'FontSize', 18);
title("Monte Carlo estimator", 'FontSize', 15);
legend(fig_legend_string, 'interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
saveas(fig, "../figures/mc_estimator", "epsc");