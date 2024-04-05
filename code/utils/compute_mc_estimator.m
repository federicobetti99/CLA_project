function ests = compute_mc_estimator(A, N, W)
% This function computes a Monte Carlo estimator for diag(inv(A))
% 
%   A: input matrix
%   N: number of Monte Carlo samples
%   W: Lanczos factor for variance reduction, can also not be passed
%
   n = size(A, 1);
   ests = zeros(n, N);
   L = ichol(A, struct('type','ict','droptol', 1e-3)); % compute incomplete Cholesky factorization once and for all
   sum = 0;
   Z = (rand(n, N) < .5) * 2 - 1; % sample from Rademacher distribution
   for l = 1 : N
       z = Z(:, l);
       [y, ~, ~, ~, ~] = pcg(A, z, [], [], L, L'); % solve iteratively linear system involving A
       if nargin == 3
           y = y - W * W' * z;
       end
       sum = sum + y .* z;
       ests(:, l) = sum / l;  % update with current point-wise multiplication
   end
end