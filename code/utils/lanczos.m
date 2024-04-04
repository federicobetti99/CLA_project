function [V, H] = lanczos(A, G, x, k)
%   This function performs k Lanczos iterations on the preconditioned matrix inv(G) A inv(G)'
%   
%   A: original matrix
%   G: preconditioner
%   x: starting vector
%   k: number of Lanczos iterations to be carried out
    factor = inv(G) * A * inv(G)';

    vold = x / norm(x, 2);
    w = factor * vold;
    alpha = dot(vold, w);
    
    V = vold;
    H = alpha;
    r = w - alpha * vold;
    
    for j = 1 : (k-1)
        e = zeros(j);
        e(j) = 1;
        beta = norm(r, 2);
        vnew = r / beta;
        V = [V, vnew];
        Hhat = [H; beta * e'];
        z = factor * vold;
        h = V' * z;
        r = z - V * h;
        H = [Hhat, h];
    end
    
end