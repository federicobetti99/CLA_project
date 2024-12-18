function [] = lanczos_mc(M, N, ks, avgs, alpha, namefile)
% This function computes a Lanczos-MC estimator for diag(inv(M)) for every
% value of Lanczos iterations stored in ks
% 
% Inputs:
%   M: input matrix of size n x n
%   N: number of Monte Carlo samples
%   ks: numbers of Lanczos iterations
%   avgs: number of averages for MC and Lanczos-MC estimates
%   alpha: width of the confidence interval
%   namefile: where to save the plot
% 
% Outputs:

    %% load matrix and compute exact quantity
    G = ichol(M);
    diaginvM = diag(inv(M));
    avg_errors = zeros(size(ks, 1), avgs, N);
    
    %% compute Lanczos estimator for a fixed value of k
    for j = 1:size(ks, 1)
        [Ts, Vs] = lanczos_iterations(M, ks(j), G);
        [est, W] = compute_lanczos_estimator(Ts(1:ks(j), 1:ks(j)), Vs(:, 1:ks(j)), G);
        for l = 1:avgs
            fprintf("Lanczos-MC estimator for k = %d: computing average %d out of %d\n", ks(j), l, avgs);
            ests = compute_mc_estimator(M, N, W);
            ests = ests + est; % add Lanczos estimate
            avg_errors(j, l, :) = vecnorm(ests-repmat(diaginvM, 1, N)) / norm(diaginvM);
        end
    end
    
    mean_errors = mean(avg_errors, 2);
    std_dev = std(avg_errors, [], 2);
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
        pl = semilogy(x, mean_errors(j, :), 'LineWidth', 3);
        hold on
        kstr = sprintf("$k = %.0f$", string(ks(j)));
        fig_legend_string = [fig_legend_string, kstr];
        h = fill(x2, inBetween, get(pl,'Color'));
        set(h, 'facealpha', .2);
        fig_legend_string = [fig_legend_string, ""];
        hold on
    end
    
    xlabel("$N$", 'interpreter', 'latex', 'FontSize', 18);
    ylabel("$\frac{\| \mathbf{d}_{\mathrm{Lanczos-MC}}^{k, N} - \mathrm{diag}(A^{-1}) \|_2}{\| \mathrm{diag}(A^{-1}) \|_2}$", ...
        'interpreter', 'latex', 'FontSize', 25);
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 13);
    a = get(gca, 'YTickLabel');
    set(gca, 'YTickLabel', a, 'fontsize', 13);
    legend(fig_legend_string, 'interpreter', 'latex');
    legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
    saveas(fig, namefile, "epsc");

end