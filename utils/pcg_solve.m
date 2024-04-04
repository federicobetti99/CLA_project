function y = pcg_solve(A, b, tol, maxit)
    L = ichol(A);
    [y, ~, ~, ~, ~] = pcg(A, b, tol, maxit, L, L');
end