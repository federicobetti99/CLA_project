%% clean
clear
close all
clc

%% import utilities and fix seed
addpath("utils")
rng("default")
rng(0)

%% define number of samples, number of averages and width of CI
N = 10000;
avgs = 10;
alpha = 0.05;

%% load matrix and compute exact quantity
matname = "mhdb416";
matfile = sprintf("../matrices/%s.mat", matname);
datastruct = load(matfile);
M = datastruct.Problem.A;
n = size(M, 1);
diaginvM = diag(M \ eye(n, n));

%% compute Monte Carlo estimator, averaging for every value of N
avg_errors = zeros(avgs, N);
for j = 1:avgs
    fprintf("Computing %d average out of %d \n", j, avgs);
    ests = compute_mc_estimator(M, N);
    avg_errors(j, :) = vecnorm(ests-repmat(diaginvM, 1, N)) / norm(diaginvM);
end

mean_errors = mean(avg_errors, 1);
std_dev = std(avg_errors, [], 1) / sqrt(avgs);
cdi = norminv(1-alpha/2);
curve1 = mean_errors + cdi * std_dev / sqrt(avgs);
curve2 = mean_errors - cdi * std_dev / sqrt(avgs);
x = (1:N);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

%% plot and save figure
fig = figure();
fig_legend_string = ["$\mathcal{O}(1/\sqrt{N})$", ""];
loglog(x, 1 ./ sqrt(x), 'LineWidth', 1);
hold on
pl = loglog(x, mean_errors, 'LineWidth', 3);
hold on
h = fill(x2, inBetween, get(pl, 'Color'));
set(h, 'facealpha', .2);
hold on
xlabel("$N$", 'interpreter', 'latex', 'fontsize', 15);
ylabel("$\frac{\| \mathbf{d}_{\mathrm{MC}}^N - \mathrm{diag}(A^{-1}) \|_2}{\| \mathrm{diag}(A^{-1}) \|_2}$", ...
    'interpreter', 'latex', 'fontsize', 18);
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 15);
a = get(gca, 'YTickLabel');
set(gca, 'YTickLabel', a, 'fontsize', 15);
legend(fig_legend_string, 'interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
namefile = sprintf("../figures/%s/mc_estimator", matname);
saveas(fig, namefile, "epsc");