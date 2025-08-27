function W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
% Block-wise computation of reachability Gramian for periodic Sylvester systems
% This implements Algorithm 1 from the paper
%
% Inputs: 
%   A_func - function handle for A(t), returns n x n matrix
%   B_func - function handle for B(t), returns n x n matrix  
%   K_func - function handle for K(t), returns n x m matrix
%   T - period
%   N - number of quadrature nodes (must be odd for Simpson rule)
%
% Output: 
%   W - Reachability Gramian (n^2 x n^2)
%
% Author: M. S. V. D. Sudarsan
% Paper: "Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"

% Validate inputs
if mod(N, 2) == 0
    error('N must be odd for composite Simpson rule');
end

% Get dimensions from initial evaluation
K0 = K_func(0);
[n, m] = size(K0);

% Setup quadrature (composite Simpson rule)
tau = linspace(0, T, N);
w = simpson_weights(N, T);

% Initialize Gramian
W = zeros(n^2, n^2);

fprintf('Computing Gramian with n=%d, m=%d, N=%d nodes...\n', n, m, N);

% Main loop over quadrature nodes
for i = 1:N
    if mod(i, ceil(N/10)) == 0
        fprintf('Processing node %d/%d (%.1f%%)\n', i, N, 100*i/N);
    end
    
    Ki = K_func(tau(i));
    M_i = zeros(n^2, m*n);
    
    % Loop over input columns and basis columns
    for k = 1:m
        zcol = Ki(:, k); % k-th column of K(tau_i)
        
        for j = 1:n
            % Create j-th unit vector
            ej = zeros(n, 1);
            ej(j) = 1;
            
            % Initial condition: Z0 = K_k * e_j^T
            Z0 = zcol * ej.';
            
            % Solve Sylvester ODE from tau_i to T
            sylv_ode = @(t, Z_vec) sylvester_rhs(t, Z_vec, A_func, B_func, n);
            
            % ODE options for accuracy
            opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
            
            % Solve ODE
            [~, Z_sol] = ode45(sylv_ode, [tau(i), T], Z0(:), opts);
            
            % Extract final value and vectorize
            Z_final = reshape(Z_sol(end, :), n, n);
            
            % Store in appropriate column of M_i
            col_idx = (k-1)*n + j;
            M_i(:, col_idx) = Z_final(:);
        end
    end
    
    % Accumulate weighted contribution to Gramian
    W = W + w(i) * (M_i * M_i');
end

fprintf('Gramian computation completed.\n');

% Check properties
sigma_min = min(real(eig(W)));
sigma_max = max(real(eig(W)));
cond_num = sigma_max / sigma_min;

fprintf('Gramian properties:\n');
fprintf('  Minimum eigenvalue: %.6e\n', sigma_min);
fprintf('  Maximum eigenvalue: %.6e\n', sigma_max);
fprintf('  Condition number: %.6e\n', cond_num);

if sigma_min > 1e-10
    fprintf('  System appears to be CONTROLLABLE\n');
else
    fprintf('  System may NOT be controllable (σ_min ≈ 0)\n');
end

end

function dZ_vec = sylvester_rhs(t, Z_vec, A_func, B_func, n)
% Right-hand side for Sylvester ODE: dZ/dt = A(t)*Z + Z*B(t)
% Input: Z_vec is vectorized form of Z
% Output: dZ_vec is vectorized form of dZ/dt

% Reshape vectorized Z back to matrix form
Z = reshape(Z_vec, n, n);

% Compute Sylvester equation RHS
dZ = A_func(t) * Z + Z * B_func(t);

% Return in vectorized form
dZ_vec = dZ(:);
end

function w = simpson_weights(N, T)
% Composite Simpson quadrature weights
% N must be odd, T is the integration interval length

if mod(N, 2) == 0
    error('N must be odd for composite Simpson rule');
end

h = T / (N - 1);
w = zeros(1, N);

% First and last points
w(1) = h/3;
w(N) = h/3;

% Interior points
for i = 2:N-1
    if mod(i-1, 2) == 0  % Even index (odd interior point)
        w(i) = 4*h/3;
    else                 % Odd index (even interior point)  
        w(i) = 2*h/3;
    end
end
end
