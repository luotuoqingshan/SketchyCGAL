function [v, xi, i] = ApproxMinEvecLanczosSE(M, n, q)
% Approximate minimum eigenvector
% Vanilla Lanczos method

q = min(q, n-1);                    % Iterations < dimension!

if isnumeric(M), M = @(x) M*x; end

aleph = zeros(q,1);                 % Diagonal Lanczos coefs
beth = zeros(q,1);                  % Off-diagonal Lanczos coefs

v = randn(n, 1);                   % First Lanczos vector is random
v = v / norm(v);
vi = v;

% First loop is to find coefficients
for i = 1 : q

    vip1 = M ( vi );			% Apply M to previous Lanczos vector
    aleph(i) = real(vi' * vip1);		% Compute diagonal coefficients
    
    if (i == 1)                     % Lanczos iteration
        vip1 = vip1 - aleph(i) * vi;
    else
        vip1 = vip1 - aleph(i) * vi - beth(i-1) * vim1;
    end
    
    beth(i) = norm( vip1 );            % Compute off-diagonal coefficients
    
    if ( abs(beth(i)) < sqrt(n)*eps ), break; end
    
    vip1 = vip1 / beth(i);        % Normalize
    
    vim1 = vi;  % update
    vi = vip1;
    
end

% i contains number of completed iterations
B = diag(aleph(1:i), 0) + diag(beth(1:(i-1)), +1) + diag(beth(1:(i-1)), -1);
[U, D] = cgal_eig(0.5*(B+B'));
[xi, ind] = min(D);
Uind1 = U(:,ind);

% Second loop is to find compute the vector (on the fly)
aleph = zeros(q,1);                 % Diagonal Lanczos coefs
beth = zeros(q,1);                  % Off-diagonal Lanczos coefs
vi = v;
v = zeros(n,1);
for i = 1 : length(Uind1)

    v = v + vi*Uind1(i);

    vip1 = M ( vi );                 % Apply M to previous Lanczos vector
    aleph(i) = real(vi' * vip1);		% Compute diagonal coefficients
    
    if (i == 1)                     % Lanczos iteration
        vip1 = vip1 - aleph(i) * vi;
    else
        vip1 = vip1 - aleph(i) * vi - beth(i-1) * vim1;
    end
    
    beth(i) = norm( vip1 );    % Compute off-diagonal coefficients
    
    % if ( abs(beth(i)) < sqrt(n)*eps ), break; end
    
    % if i >= numit, warning('numerical error in Lanczos'); break; end
    
    vip1 = vip1 / beth(i);          % Normalize
    
    vim1 = vi;  % update
    vi = vip1;
        
end

i = 2*i; % we looped twice

% Next lines are unnecessary in general, but I observed numerical errors in
% norm(v) at some experiments, so let's normalize it for robustness. 
nv = norm(v);
xi = xi*nv;
v = v/nv;
end

function [V,D] = cgal_eig(X)
% Eig in Lanczos based LMO solver sometimes fall into numerical issues. 
% This function replaces eig with a SVD based solver, in case eig does not
% converge. 
try
    [V,D]       = eig(X,'vector');
catch 
	warning('eig did not work. Using the svd based replacement instead.');
    [V,D,W]     = svd(X);
    D           = diag(D).' .* sign(real(dot(V,W,1)));
    [D,ind]     = sort(D);
    V           = V(:,ind);
end
end