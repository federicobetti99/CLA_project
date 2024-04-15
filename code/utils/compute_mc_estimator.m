function ests = compute_mc_estimator(A, N, W)
% This function computes a Monte Carlo estimator for diag(inv(A))
% 
%   A: input matrix
%   N: number of Monte Carlo samples
%   W: Lanczos factor for variance reduction, can also not be passed
%
   n = size(A, 1);
   ests = zeros(n, N+1);
   L = ichol(A, struct('type', 'ict', 'droptol', 1e-3)); % compute incomplete Cholesky factorization
   Z = randsrc(n, N); % sample from Rademacher distribution
   for l = 1 : N
       z = Z(:, l);
       [y, ~, ~, ~, ~] = pcg(A, z, [], [], L, L'); % solve iteratively linear system involving A
       if nargin == 3
           y = y - W * W' * z;
       end
       ests(:, l+1) = ests(:, l) + y .* z;  % update with current point-wise multiplication
   end
   ests = ests(:, 2:end) ./ (1:N);
end