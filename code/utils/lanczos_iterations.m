function [H, U] = lanczos_iterations(M, m, G)
% Perform m steps of the Lanczos process, storing the full basis U
% with reorthogonalization on symmetric M with starting vector x.
% M is always assumed to be square.
%
% Inputs:
%   M: input matrix of size n x n
% % % %   m: number of iterations
%   G: preconditioner, can also not be passed
%
% Outputs:
%   H: tridigonal matrix obtained from Lanczos iterations 
%   U: matrix of size n x m whose columns are an orthogonal basis for the
%      m-th Krylov subspace
    
    n = size(M, 1);
    x = rand(n, 1);
    x = x / norm(x);
    
    H = zeros(m+1,m+1);
    U = x;
    beta = 0;
    m = min(n, m);
    
    for j = 1:m
        if nargin == 2
            z = M * U(:,end);
        else
            z =  G' \ U(:, end);
            z = M * z;
            z = G \ z;
        end

        alpha = U(:,end)' * z;
        u = z - U(:,end) * alpha;
        if j > 1
            u = u - U(:,end-1) * beta;
        end
        
        % Reorthognalization
        alphas = U' * u;
        u = u - U * alphas;
        alpha = alpha + alphas(end);
        
        beta = norm(u);
        u = u / beta;
        U = [U, u];
        H(j, j) = alpha;
        H(j+1, j) = beta;
        H(j, j+1) = beta;
    end

end