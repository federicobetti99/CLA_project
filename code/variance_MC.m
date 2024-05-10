% A code to understand the variance for the MC estimator.

n = 101;
N = 100;
avgs = 10;
alpha = 0.05;
Np = 5;
e = zeros(N, Np);
F = @(A) norm(inv(A)-diag(inv(A)),'fro');

fig = figure();

for k=1:Np
    p = 10^k;
    A = make_A(p, n);
    disp('Frobenius norm of inv(A) - diag(inv(A))')
    disp(F(A));
    disp('cond(A)')
    conditionA = cond(full(A));
    disp(conditionA);
    d = diag(inv(A));
    avg_errors = zeros(avgs, N);
    for j = 1:avgs
        fprintf("Computing %d average out of %d \n", j, avgs);
        ests = compute_mc_estimator(A, N);
        avg_errors(j, :) = vecnorm(ests-repmat(d, 1, N)) / norm(d);
    end
    
    mean_errors = mean(avg_errors, 1);
    std_dev = std(avg_errors, [], 1) / sqrt(avgs);
    cdi = norminv(1-alpha/2);
    curve1 = mean_errors + cdi * std_dev / sqrt(avgs);
    curve2 = mean_errors - cdi * std_dev / sqrt(avgs);
    x = (1:N);
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];

    pl = loglog(x, mean_errors, 'LineWidth', 3, 'DisplayName', ['$\|A^{-1}-\mathrm{diag}(A^{-1})\|_F=' num2str(F(A), '%.0e') '$']);
    hold on
    h = fill(x2, inBetween, get(pl, 'Color'));
    set(h, 'facealpha', .2, 'HandleVisibility', 'off');
    hold on
end

legend('Interpreter', 'latex');
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 15);
a = get(gca, 'YTickLabel');
set(gca, 'YTickLabel', a, 'fontsize', 15);  
xlabel("$N$", 'interpreter', 'latex', 'fontsize', 15);
ylabel("$\frac{\| \mathbf{d}_{\mathrm{MC}}^N - \mathrm{diag}(A^{-1}) \|_2}{\| \mathrm{diag}(A^{-1}) \|_2}$", ...
    'interpreter', 'latex', 'fontsize', 18);
title(['$\kappa(A)=' num2str(conditionA, '%.2e') '$'], 'interpreter', 'latex', 'fontsize', 15)
hold off;

function A = make_A(p, N)
    e = ones(N,1);
    A1 = spdiags([e 2*e e], -1:1, N, N);
    A = 6.55*p*A1/(3);
end