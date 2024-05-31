function [] = lanczos(M, k, namefile, eig_namefile)
% This function computes a Lanczos based estimator 
% for every value j from 1 to k
% 
% Inputs:
%   M: input matrix of size n x n
%   k: number of Lanczos iterations
%   namefile: where to save the plot
%   eig_namefile: where to save the eigenvalues plot
% 
% Outputs:

    %% compute incomplete Cholesky factor and compute exact quantity
    G = ichol(M);
    n = size(M, 1);
    diaginvM = diag(M \ eye(n, n));
    
    errors = zeros(k, 1);
    precond_errors = zeros(k, 1);
    
    %% compute k Lanczos iterations without preconditioner
    [Ts, Vs] = lanczos_iterations(M, k);
    
    %% compute Lanczos estimator for every value of j = 1, ..., k
    for j = 1:k
        [est, ~] = compute_lanczos_estimator(Ts(1:j, 1:j), Vs(:, 1:j));
        errors(j, :) = vecnorm(est-diaginvM) / norm(diaginvM);
    end
    
    %% compute k Lanczos iterations with preconditioner
    [Ts, Vs] = lanczos_iterations(M, k, G);
    
    %% compute Lanczos estimator for every value of j = 1, ..., k
    for j = 1:k
        [est, ~] = compute_lanczos_estimator(Ts(1:j, 1:j), Vs(:, 1:j), G);
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
    xlabel("$k$", 'interpreter', 'latex', 'FontSize', 18);
    ylabel("$\frac{\| \mathbf{d}_{\mathrm{Lanczos}}^k- \mathrm{diag}(A^{-1}) \|_2}{\| \mathrm{diag}(A^{-1}) \|_2}$", ...
        'interpreter', 'latex', 'FontSize', 25);
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 13);
    a = get(gca, 'YTickLabel');
    set(gca, 'YTickLabel', a, 'fontsize', 13);
    legend(fig_legend_string, 'interpreter', 'latex');
    legend('Location', 'northeast', 'FontSize', 18, 'NumColumns', 1);
    saveas(fig, namefile, "epsc");

    fig = figure();
    dM = eig(M);
    semilogy(1 ./ dM, 'LineWidth', 3);
    hold on
    ylabel('Inverse ordered eigenvalues', 'fontsize', 20);
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 13);
    a = get(gca, 'YTickLabel');
    set(gca, 'YTickLabel', a, 'fontsize', 13);
    saveas(fig, eig_namefile, "epsc");
    
end