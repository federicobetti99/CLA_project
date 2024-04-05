function [est, W] = compute_lanczos_estimator(G, T, V)
%COMPUTE_LANCZOS_ESTIMATOR Summary of this function goes here
%   Detailed explanation goes here
    L = chol(T);
    W = (L \ V')';
    W = G' \ W;
    est = power(vecnorm(W, 2, 2), 2);
end

