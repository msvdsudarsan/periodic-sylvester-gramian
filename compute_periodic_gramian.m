
---

Next: the **MATLAB functions**. Copy each block into a separate file with the exact filename shown above the block (for example, copy the first block into a file named `compute_periodic_gramian.m`).

---

### `compute_periodic_gramian.m`
```matlab
function W = compute_periodic_gramian(A_func, B_func, K_func, T, N)
% COMPUTE_PERIODIC_GRAMIAN Computes reachability Gramian for periodic Sylvester systems
%
% Syntax:
%   W = compute_periodic_gramian(A_func, B_func, K_func, T, N)
%
% Inputs:
%   A_func - Function handle for A(t), returns n x n matrix
%   B_func - Function handle for B(t), returns n x n matrix  
%   K_func - Function handle for K(t), returns n x m matrix
%   T      - Period of the system
%   N      - Number of quadrature nodes (odd number recommended)
%
% Output:
%   W      - Reachability Gramian (n^2 x n^2 matrix)
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com
% Date: August 2025

    % Validate inputs
    if nargin < 5
        error('All five inputs are required');
    end
    
    % Get system dimensions
    K0 = K_func(0);
    [n, m] = size(K0);
    
    % Setup quadrature (composite Simpson's rule)
    if mod(N, 2) == 0
        N = N + 1; % Ensure odd number for Simpson's rule
        warning('N adjusted to %d for Simpson''s rule', N);
    end
    
    tau = linspace(0, T, N);
    w = simpson_weights(N, T);
    
    % Initialize Gramian
    W = zeros(n^2, n^2);
    
    % Progress indicator
    fprintf('Computing Gramian with %d quadrature nodes...\n', N);
    
    % Main computation loop
    for i = 1:N
        if mod(i, 20) == 0
            fprintf('Progress: %d/%d nodes processed\n', i, N);
        end
        
        % Initialize M_i matrix for current quadrature node
        M_i = zeros(n^2, m*n);
        
        % Process each input column
        for k = 1:m
            % Extract k-th column of K(tau_i)
            K_col = K_func(tau(i));
            K_k = K_col(:, k);
            
            % Initial condition for Sylvester ODE: Z0 = K_k * ones(1,n)
            Z0 = K_k * ones(1, n);
            
            % Solve Sylvester ODE from tau(i) to T
            Z_final = solve_sylvester_ode(Z0, tau(i), T, A_func, B_func);
            
            % Store vectorized columns in M_i
            col_start = (k-1)*n + 1;
            col_end = k*n;
            M_i(:, col_start:col_end) = reshape(Z_final, n^2, n);
        end
        
        % Accumulate Gramian contribution
        W = W + w(i) * (M_i * M_i');
    end
    
    fprintf('Gramian computation completed.\n');
    
    % Check properties
    sigma_min = min(real(eig(W)));
    sigma_max = max(real(eig(W)));
    condition_number = sigma_max / sigma_min;
    
    fprintf('Gramian properties:\n');
    fprintf('  Minimum eigenvalue: %.6e\n', sigma_min);
    fprintf('  Condition number: %.2e\n', condition_number);
    
    if sigma_min > 1e-12
        fprintf('  System appears to be controllable.\n');
    else
        fprintf('  System may not be controllable (small min eigenvalue).\n');
    end
end

function Z_final = solve_sylvester_ode(Z0, t0, tf, A_func, B_func)
    % Solve dZ/dt = A(t)*Z + Z*B(t) from t0 to tf
    
    n = size(Z0, 1);
    
    % ODE function handle
    sylv_ode = @(t, z_vec) sylvester_rhs_vec(t, z_vec, A_func, B_func, n);
    
    % Initial condition as vector
    z0_vec = Z0(:);
    
    % Solve ODE with tight tolerances
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
    [~, z_sol] = ode45(sylv_ode, [t0, tf], z0_vec, options);
    
    % Extract final solution and reshape
    z_final_vec = z_sol(end, :)';
    Z_final = reshape(z_final_vec, n, n);
end

function dz_dt = sylvester_rhs_vec(t, z_vec, A_func, B_func, n)
    % Right-hand side for vectorized Sylvester equation
    Z = reshape(z_vec, n, n);
    dZ_dt = A_func(t) * Z + Z * B_func(t);
    dz_dt = dZ_dt(:);
end
