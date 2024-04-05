%% clean
clear
close all
clc

%% load matrix and compute exact quantity, define number of Lanczos iterations
matname = "mhdb416";
matfile = sprintf("../matrices/%s.mat", matname);
datastruct = load(matfile);
M = datastruct.Problem.A;
G = ichol(M, struct('type','ict','droptol', 1e-3));
diaginvM = diag(inv(M));
k = 200;
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
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 15);
a = get(gca, 'YTickLabel');
set(gca, 'YTickLabel', a, 'fontsize', 15);
namefile = sprintf("../figures/%s/lanczos_estimator", matname);
saveas(fig, namefile, "epsc");