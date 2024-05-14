function [] = relgap(M, namefile)

    % Eigenvalues of M
    d_M = eig(M);
    
    % Preconditioner
    G = ichol(M);
    B = G' \ M;
    X = G \ B;
    
    % Eigenvalues of preconditioned M
    d_X = sort(real(eig(full(X))));
    
    % Plots
    fig_legend_string = ["$A$", "$G^{-1}AG^{-T}$"];
   
    % Reglap
    fig = figure();
    semilogy(d_M, 'LineWidth', 3);
    hold on
    semilogy(d_X, 'LineWidth', 3);
    hold on
    ylabel('Ordered eigenvalues', 'fontsize', 18);
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 13);
    a = get(gca, 'YTickLabel');
    set(gca, 'YTickLabel', a, 'fontsize', 13);
    legend(fig_legend_string, 'interpreter', 'latex');
    legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 1);
    saveas(fig, namefile, "epsc");
end