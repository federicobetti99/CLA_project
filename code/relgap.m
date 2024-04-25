function [] = relgap(M, namefile)

    % Eigenvalues of M
    d_M = eig(M);
    
    % Relative gaps for M
    rg_M = (diff(d_M))/(d_M(end) - d_M(1));
    
    % Preconditionner
    G = ichol(M);
    B = G' \ M;
    X = G \ B;
    
    % Eigenvalues of preconditioned M
    d_X = eig(full(X));
    
    % Relative gaps for preconditioned M
    rg_X = (diff(d_X)) / (d_X(end) - d_X(1));
    
    % Plots
    fig_legend_string = ["$A$", "$G^{-1}AG^{-T}$"];
    
    % Spectrums
    semilogy(d_M, 'LineWidth', 3);
    hold on
    semilogy(d_X, 'LineWidth', 3);
    hold on
    ylabel('$\mathrm{sp}(A)$', 'Interpreter', 'latex', 'fontsize', 18);
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 13);
    a = get(gca, 'YTickLabel');
    set(gca, 'YTickLabel', a, 'fontsize', 13);
    legend(fig_legend_string, 'interpreter', 'latex');
   
    % Reglap
    fig = figure();
    plot(rg_M);
    hold on
    plot(rg_X)
    hold on
    ylabel('Relative gap', 'fontsize', 18);
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 13);
    a = get(gca, 'YTickLabel');
    set(gca, 'YTickLabel', a, 'fontsize', 13);
    legend(fig_legend_string, 'interpreter', 'latex');
    legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
    saveas(fig, namefile, "epsc");
end