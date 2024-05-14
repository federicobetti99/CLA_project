function [est, W] = compute_lanczos_estimator(T, V, G)
% This function computes a Lanczos-based estimator diag(inv(M)),
% with M spd of size n x n, given the orthogonal matrix V, 
% the tridigonal matrix T and the preconditioner G
% 
% Inputs:
%   T: tridigonal matrix obtained from Lanczos iterations
%   V: matrix of size n x k whose columns are an orthogonal basis for the
%      k-th Krylov subspace
%   G: preconditioner possibly used during the Lanczos iterates, can also not be passed
% 
% Outputs:
%   est: an estimate for diag(inv(A))
%   W: approximation of the Cholesky factor of inv(A)

    L = chol(T, 'lower'); % compute Cholesky factor of T
    W = (L \ V')';  % cheap solve as L is lower triangular with zeros below the first subdiagonal 
    if nargin == 3
        W = G' \ W; % apply the preconditioner if it is passed
    end
    est = sum(W .^ 2, 2); % compute estimator as diag(W W'), no need to assemble it explicitly

end

