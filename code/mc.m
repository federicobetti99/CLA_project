function [] = mc(M, N, avgs, alpha, namefile)
% This function computes a MC estimator for diag(inv(M))
% 
% Inputs:
%   M: input matrix of size n x n
%   N: number of Monte Carlo samples
%   avgs: number of averages for MC and Lanczos-MC estimates
%   alpha: width of the confidence interval
%   namefile: where to save the plot
% 
% Outputs:  

    %% compute exact quantity
    n = size(M, 1);
    diaginvM = diag(M \ eye(n, n));
    
    %% compute Monte Carlo estimator, averaging for every value of N
    avg_errors = zeros(avgs, N);
    for j = 1:avgs
        fprintf("MC estimator: computing average %d out of %d\n", j, avgs);
        ests = compute_mc_estimator(M, N);
        avg_errors(j, :) = vecnorm(ests-repmat(diaginvM, 1, N)) / norm(diaginvM);
    end
    
    mean_errors = mean(avg_errors, 1);
    std_dev = std(avg_errors, [], 1);
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
    xlabel("$N$", 'interpreter', 'latex', 'fontsize', 18);
    ylabel("$\frac{\| \mathbf{d}_{\mathrm{MC}}^N - \mathrm{diag}(A^{-1}) \|_2}{\| \mathrm{diag}(A^{-1}) \|_2}$", ...
        'interpreter', 'latex', 'fontsize', 25);
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 18);
    a = get(gca, 'YTickLabel');
    set(gca, 'YTickLabel', a, 'fontsize', 18);
    legend(fig_legend_string, 'interpreter', 'latex');
    legend('Location', 'northeast', 'FontSize', 18, 'NumColumns', 1);
    saveas(fig, namefile, "epsc");

end