%% clean
clear
close all
clc

%% import utilities
addpath("utils")

%% load matrix and compute exact quantity, define number of Lanczos iterations
matname = "mesh3em5";
matfile = sprintf("../matrices/%s.mat", matname);

datastruct = load(matfile);
M = datastruct.Problem.A;
G = ichol(M);
n = size(M, 1);
diaginvM = diag(M \ eye(n, n));

switch matname
    case "nos3"
        k = 500;
    case "mesh3em5"
        k = 200;
    case "mhdb416"
        k = 400;
    otherwise
        disp("Invalid matrix name passed");
end

errors = zeros(k, 1);
precond_errors = zeros(k, 1);

%% compute k Lanczos iterations without preconditioner
[Ts, Vs] = lanczos_iterations(M, k);

%% compute Lanczos estimator for every value of j = 1, ..., k
for j = 1:k
    est = compute_lanczos_estimator(Ts(1:j, 1:j), Vs(:, 1:j));
    errors(j, :) = vecnorm(est-diaginvM) / norm(diaginvM);
end

%% compute k Lanczos iterations with preconditioner
[Ts, Vs] = lanczos_iterations(M, k, G);

%% compute Lanczos estimator for every value of j = 1, ..., k
for j = 1:k
    est = compute_lanczos_estimator(Ts(1:j, 1:j), Vs(:, 1:j), G);
    precond_errors(j, :) = vecnorm(est-diaginvM) / norm(diaginvM);
end

%% plot and save figure
fig = figure();
x = (1:k)';
semilogy(x, errors, 'LineWidth', 3);
hold on
semilogy(x, precond_errors, 'LineWidth', 3);
hold on
fig_legend_string = ["Lanczos", "preconditioned Lanczos"];
hold on
xlabel("$k$", 'interpreter', 'latex', 'FontSize', 15);
ylabel("$\frac{\| \mathbf{d}_{\mathrm{Lanczos}}^k- \mathrm{diag}(A^{-1}) \|_2}{\| \mathrm{diag}(A^{-1}) \|_2}$", ...
    'interpreter', 'latex', 'FontSize', 18);
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 13);
a = get(gca, 'YTickLabel');
set(gca, 'YTickLabel', a, 'fontsize', 13);
legend(fig_legend_string, 'interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
namefile = sprintf("../figures/%s/lanczos_estimator", matname);
saveas(fig, namefile, "epsc");