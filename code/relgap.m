function [] = relgap(M, namefile)

    % Eigenvalues of M
    d_M = eig(M);
   
    % Reglap
    fig = figure();
    semilogy(1 ./ d_M, 'LineWidth', 3);
    hold on
    ylabel('Inverse ordered eigenvalues', 'fontsize', 20);
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 13);
    a = get(gca, 'YTickLabel');
    set(gca, 'YTickLabel', a, 'fontsize', 13);
    saveas(fig, namefile, "epsc");
end