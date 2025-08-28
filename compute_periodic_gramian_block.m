function W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
%COMPUTE_PERIODIC_GRAMIAN_BLOCK Block-wise computation of reachability Gramian
%
% W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
%
% Computes the reachability Gramian for the periodic Sylvester system:
%   dX/dt = A(t)*X + X*B(t) + K(t)*U(t)
%
% Inputs:
%   A_func - function handle for A(t), returns n×n matrix
%   B_func - function handle for B(t), returns n×n matrix  
%   K_func - function handle for K(t), returns n×m matrix
%   T      - period
%   N      - number of quadrature nodes (should be odd for Simpson's rule)
%
% Output:
%   W      - n²×n² reachability Gramian matrix
%
% Algorithm avoids forming n²×n² Kronecker matrices by solving
% mn Sylvester ODEs of size n×n, reducing complexity from O(Nn⁶) to O(Nn³m).

% Get dimensions
K0 = K_func(0);
[n, m] = size(K0);

fprintf('Computing Gramian: n=%d, m=%d, period T=%.3f, N=%d nodes\n', n, m, T, N);

% Quadrature setup (composite Simpson's rule)
if mod(N, 2) == 0
    N = N + 1; % Make N odd for Simpson's rule
    fprintf('Adjusted N to %d (odd) for Simpson''s rule\n', N);
end

tau = linspace(0, T, N);
w = simpson_weights(N, T);

% Initialize Gramian
W = zeros(n^2, n^2);

% ODE solver options
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

% Main loop over quadrature nodes
progress_step = max(1, floor(N/10));
for i = 1:N
    if mod(i-1, progress_step) == 0 || i == N
        fprintf('Processing node %d/%d (%.1f%%)\n', i, N, 100*i/N);
    end
    
    % Skip boundary node at tau = T to avoid ODE integration issues
    if abs(tau(i) - T) < 1e-12
        fprintf('Skipping boundary node at tau(%d) = %.6f ≈ T\n', i, tau(i));
        continue;
    end
    
    Ki = K_func(tau(i));
    M_i = zeros(n^2, m*n);
    
    % Loop over input columns and basis columns
    for k = 1:m
        zcol = Ki(:, k); % n×1 column vector
        
        for j = 1:n
            % Create j-th basis vector
            ej = zeros(n, 1); 
            ej(j) = 1;
            
            % Initial condition: Z0 = zcol * ej^T (rank-1 matrix)
            Z0 = zcol * ej';
            
            % Solve Sylvester ODE: dZ/dt = A(t)*Z + Z*B(t)
            % from tau(i) to T with initial condition Z(tau(i)) = Z0
            sylv_ode = @(t, Z_vec) sylvester_rhs(t, Z_vec, A_func, B_func);
            [~, Z_sol] = ode45(sylv_ode, [tau(i), T], Z0(:), opts);
            
            % Extract final value Z(T) and assign to column
            Z_final = reshape(Z_sol(end, :), n, n);
            col_idx = (k-1)*n + j;
            M_i(:, col_idx) = Z_final(:);
        end
    end
    
    % Accumulate weighted contribution to Gramian
    W = W + w(i) * (M_i * M_i');
end

fprintf('Gramian computation completed successfully!\n');
end

function dZ_vec = sylvester_rhs(t, Z_vec, A_func, B_func)
%SYLVESTER_RHS Right-hand side for Sylvester ODE dZ/dt = A(t)*Z + Z*B(t)
n = round(sqrt(length(Z_vec)));
Z = reshape(Z_vec, n, n);
dZ = A_func(t) * Z + Z * B_func(t);
dZ_vec = dZ(:);
end

function w = simpson_weights(N, T)
%SIMPSON_WEIGHTS Composite Simpson quadrature weights
% N must be odd for composite Simpson rule
if mod(N, 2) == 0
    error('N must be odd for composite Simpson rule');
end

h = T / (N - 1);
w = ones(1, N);

% Simpson's rule weights: 1, 4, 2, 4, 2, ..., 4, 1
w(2:2:N-1) = 4;  % Odd indices (except first and last)
w(3:2:N-2) = 2;  % Even indices (except first and last)
w = w * h / 3;
end
