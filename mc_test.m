n = 100;
N = 100;
M = load("nos3.mat");
diaginvM = diag(inv(M));
ests = mc_estimator(M, N);

fig_legend_string = ["", "$\mathcal{O}(1/\sqrt{N})$"];
loglog((1:N), vecnorm(ests-repmat(diaginvM, 1, N)) / norm(diaginvM), 'LineWidth', 4)
hold on
loglog((1:N), 10 ./ sqrt((1:N)), 'LineWidth', 2);
xlabel("N", 'FontSize', 12);
ylabel("$\frac{\vert \vert \mathbf{d}_{\mathrm{MC}}^N - diag(A^{-1}) \vert \vert_2}{\vert \vert diag(A^{-1}) \vert \vert_2}$", ...
    'interpreter', 'latex', 'FontSize', 15);
title("Monte Carlo estimator", 'FontSize', 15);
legend(fig_legend_string, 'interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
saveas(fig, "../figures/comparison", "epsc");