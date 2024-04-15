function [est, W] = compute_lanczos_estimator(T, V, G)
%COMPUTE_LANCZOS_ESTIMATOR Summary of this function goes here
%   Detailed explanation goes here
    L = chol(T);
    W = (L' \ V')';
    if nargin == 3
        W = G' \ W;
    end
    est = sum(W .^ 2, 2);
end

