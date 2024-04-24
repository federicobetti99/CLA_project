matname = "nos3";
matfile = sprintf("../matrices/%s.mat", matname);
datastruct = load(matfile);
M = datastruct.Problem.A;

% Eigenvalues of M
d_M = eig(M);

% Relative gaps for M
rg_M = (diff(d_M))/(d_M(end) - d_M(1));

% Preconditionner
G = ichol(M);
Ginv = inv(G);
B = G'\M;
X = G\B;

% Eigenvalues of preconditionned M
d_X = eig(full(X));

% Relative gaps for preconditionned M
rg_X = (diff(d_X))/(d_X(end) - d_X(1));

% Plots

fig_legend_string = ["$A$", "$G^{-1}AG^{-T}$"];

% Spectrums
fig= figure();

semilogy(d_M, 'b-');
hold on,
semilogy(d_X, 'r-')
hold on
ylabel('relgap');
title(matname);
legend(fig_legend_string, 'interpreter', 'latex');

ylabel('$\mathrm{sp}(A)$', 'Interpreter', 'latex');

% Reglap
fig = figure();
plot(rg_M, 'b-');
hold on,
plot(rg_X, 'r-')
hold on
ylabel('relgap');
title(matname);
legend(fig_legend_string, 'interpreter', 'latex');