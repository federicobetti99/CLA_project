function [T, Q] = lanczos_iterations(mat, G, k)
% This function uses the Lanczos algorithm with full
% reorthogonalization to compute k x k symmetric tridiagonal
% matrix T that approximates mat up to rank k with respect to
% transformation Q. That is, mat = Q * T * Q'.
%               
% Inputs:       mat is a symmetric or Hermitian N x N matrix
%               G is the preconditioner
%               k is the number of Lanczos iterations to perform
%
% Outputs:      T is a k x k symmetric tridiagonal matrix. When k < N, T is
%               approximately similar to mat up to rank k. When k = N, T is
%               similar (in the linear algebra sense) to mat. That is, when
%               k = N, T has identical eigenvalues to mat
%
%               Q is the N x k similarity transformation such that
%               mat = Q * T * Q'. Note that equality holds only when k = N.
%               For k < N, mat ~ Q * T * Q'
%

    n = size(mat, 1);
    v = randn(n, 1);  % random initial vector
    
    % Initialize variables
    Q = nan(n, k);
    q = v / norm(v);
    Q(:, 1) = q;
    d = nan(k, 1);
    od = nan(k-1, 1);
    
    % Perform Lanczos iterations
    for i = 1:k
        z = G' \ q;
        z = mat * z;
        z = G \ z;
        d(i) = q' * z;
        
        z = z - Q(:, 1:i) * (Q(:, 1:i)' * z);
        z = z - Q(:, 1:i) * (Q(:, 1:i)' * z);
        
        if (i ~= k)
            od(i) = norm(z);
            q = z / od(i);
            Q(:, i+1) = q;
        end
    end

    % Construct T
    T = diag(d) + diag(od, -1) + diag(od, 1);
    
end