function [H, U] = lanczos_iterations(A, m, G)
% Perform m steps of the Lanczos process, storing the full basis U
% with reorthogonalization on Hermitian A with starting vector x.
% A is always assumed to be square.
    
    n = size(A, 1);
    x = rand(n, 1);
    x = x / norm(x);
    
    H = zeros(m+1,m+1);
    U = x;
    beta = 0;
    m = min(n, m);
    
    for j = 1:m
        if nargin == 2
            z = A * U(:,end);
        else
            z =  G' \ U(:, end);
            z = A * z;
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