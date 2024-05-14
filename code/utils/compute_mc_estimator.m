function ests = compute_mc_estimator(M, N, W)
% This function computes a Monte Carlo estimator for diag(inv(M)) for a
% number of random Rademacher vectors ranging from 1 to N
% 
% Inputs:
%   M: input matrix of size n x n
%   N: number of Monte Carlo samples
%   W: Lanczos factor for variance reduction, can also not be passed
% 
% Outputs:
%   ests: a matrix n x N with N estimates for diag(inv(M))

   n = size(M, 1);
   ests = zeros(n, N+1);
   L = ichol(M, struct('type', 'ict', 'droptol', 1e-3)); % compute once and for all incomplete Cholesky factorization
   Z = randsrc(n, N); % sample from Rademacher distribution
   for l = 1 : N
       z = Z(:, l);
       [y, ~, ~, ~, ~] = pcg(M, z, [], [], L, L'); % solve with PCG
       if nargin == 3
           y = y - W * W' * z;  % needed for Lanczos-MC estimator
       end
       ests(:, l+1) = ests(:, l) + y .* z;  % update estimate
   end
   ests = ests(:, 2:end) ./ (1:N);

end