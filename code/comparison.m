function [] = comparison(M, k, N, avgs, alpha, namefile)
% This function carries out a comparison between the three estimators
% investigated in the project for the estimation of diag(inv(M))
% 
% Inputs:
%   M: input matrix of size n x n
%   k: number of Lanczos iterations
%   N: number of Monte Carlo samples
%   avgs: number of averages for MC and Lanczos-MC estimates
%   alpha: width of the confidence interval
%   namefile: where to save the plot
% 
% Outputs:

    %% compute incomplete Cholesky factor
    G = ichol(M);
    n = size(M, 1);
    diaginvM = diag(M \ eye(n, n));
    
    %% compute Monte Carlo estimator, averaging for every value of N
    avg_mc_errors = zeros(avgs, N);
    for j = 1:avgs
        ests = compute_mc_estimator(M, N);
        avg_mc_errors(j, :) = vecnorm(ests-repmat(diaginvM, 1, N)) / norm(diaginvM);
    end
    
    mean_mc_errors = mean(avg_mc_errors, 1);
    std_mc_dev = std(avg_mc_errors, [], 1) / sqrt(avgs);
    cdi = norminv(1-alpha/2);
    curve1 = mean_mc_errors + cdi * std_mc_dev / sqrt(avgs);
    curve2 = mean_mc_errors - cdi * std_mc_dev / sqrt(avgs);
    inBetween_mc = [curve1, fliplr(curve2)];
    
    %% compute Lanczos estimator for k iterations
    [Ts, Vs] = lanczos_iterations(M, k, G);
    [lanczos_est, W] = compute_lanczos_estimator(Ts(1:k, 1:k), Vs(:, 1:k), G);
    lanczos_error = vecnorm(lanczos_est-diaginvM) / norm(diaginvM);
    
    %% compute Lanczos MC estimator
    avg_errors = zeros(avgs, N);
    for l = 1:avgs
        ests = compute_mc_estimator(M, N, W);
        ests = ests + lanczos_est; % add Lanczos estimate
        avg_errors(l, :) = vecnorm(ests-repmat(diaginvM, 1, N)) / norm(diaginvM);
    end
    
    mean_lanczos_mc_errors = mean(avg_errors, 1);
    std_lanczos_mc_dev = std(avg_errors, [], 1) / sqrt(avgs);
    cdi = norminv(1-alpha/2);
    curve1 = mean_lanczos_mc_errors + cdi * std_lanczos_mc_dev / sqrt(avgs);
    curve2 = mean_lanczos_mc_errors - cdi * std_lanczos_mc_dev / sqrt(avgs);
    inBetween_lanczos_mc = [curve1, fliplr(curve2)];
    
    %% plot and save figure
    fig = figure();
    fig_legend_string = ["MC", "", "LanczosMC", "", "Lanczos"];
    x = (1:N);
    x2 = [x, fliplr(x)];
    pl = semilogy(x, mean_mc_errors, 'LineWidth', 3);
    hold on
    h = fill(x2, inBetween_mc, get(pl, 'Color'));
    set(h, 'facealpha', .2);
    hold on
    pl = semilogy(x, mean_lanczos_mc_errors, 'LineWidth', 3);
    hold on
    h = fill(x2, inBetween_lanczos_mc, get(pl, 'Color'));
    set(h, 'facealpha', .2);
    hold on
    pl = semilogy(x, lanczos_error * ones(size(x)), 'LineWidth', 3);
    hold on
    xlabel("$N$", 'interpreter', 'latex', 'fontsize', 15);
    ylabel("Relative $\ell_2$ error", ...
        'interpreter', 'latex', 'fontsize', 18);
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 13);
    a = get(gca, 'YTickLabel');
    set(gca, 'YTickLabel', a, 'fontsize', 13);
    legend(fig_legend_string, 'interpreter', 'latex');
    legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
    saveas(fig, namefile, "epsc");
    
end    


