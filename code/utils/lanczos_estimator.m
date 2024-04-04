function [V, T] = lanczos_estimator(A, G, m)
%   This function performs k Lanczos iterations on the preconditioned matrix inv(G) A inv(G)'
%   
%   A: original matrix
%   G: preconditioner
%   x: starting vector
%   k: number of Lanczos iterations to be carried out
    n = size(A, 1);
    factor = inv(G) * A * inv(G)';
    V = zeros(n, m+1);
    alpha = zeros(m+2);
    beta = zeros(m+2);
    V(:, 2) = rand(n, 1);
    V(:, 2) = V(:, 2) / norm(V(:, 2), 2);
    beta(2) = 0;
    
    for j=2:m+2
        w = factor * V(:, j) - beta(j) * V(:, j-1);
        alpha(j) = w' * V(:, j);
        w = w - alpha(j) * V(:,j);
        beta(j+1) = norm(w, 2);
        V(:, j+1) = w / beta(j+1);
    end
    
    T = zeros(m+1, m+1);
    for i=2:m+1
        T(i-1, i-1) = alpha(i);
        T(i-1, i) = beta(i+1);
        T(i, i-1) = beta(i+1);
    end 

    T(m+1, m+1) = alpha(m+2);
    V = V(:, 2:end-1);
    
end