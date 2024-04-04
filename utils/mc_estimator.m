function ests = mc_estimator(A, N)
% This function computes a Monte Carlo estimator for diag(inv(A))
% 
%   A: input matrix
%   N: number of Monte Carlo samples    
%
   n = size(A, 1);
   ests = [];
   sum = 0;
   for l = 1 : N
       z = ((rand(1, n) < .5) * 2 - 1)';  % sample from Rademacher distribution
       y = pcg_solve(A, z, 1e-10, 200); % solve iteratively linear system involving A
       sum = sum + y .* z;
       ests = [ests, sum / l];  % update with current point-wise multiplication
   end
end