function W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
% COMPUTE_PERIODIC_GRAMIAN_BLOCK Block-wise computation of reachability Gramian
% 
% This function computes the reachability Gramian for periodic Sylvester 
% matrix systems using a structure-exploiting block propagation algorithm
% that avoids explicit formation of n^2 x n^2 Kronecker matrices.
%
% SYNTAX:
%   W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
%
% INPUTS:
%   A_func - Function handle for A(t) matrix (n x n)
%   B_func - Function handle for B(t) matrix (n x n) 
%   K_func - Function handle for K(t) matrix (n x m)
%   T      - Period length
%   N      - Number of quadrature nodes (must be odd for Simpson rule)
%
% OUTPUT:
%   W      - Reachability Gramian (n^2 x n^2)
%
% ALGORITHM:
%   Uses block-wise propagation to reduce computational complexity from
%   O(N*n^6) to O(N*n^3*m) by solving n*m Sylvester ODEs of size n x n
%   instead of one ODE of size n^2 x n^2.
%
% REFERENCE:
%   M.S.V.D. Sudarsan, "Controllability and Efficient Gramian Computation 
%   for Periodic Sylvester Matrix Systems", Applied Mathematics Letters (2025)
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com

% Input validation
if nargin < 5
    error('All five inputs are required: A_func, B_func, K_func, T, N');
end

if mod(N,2) == 0
    error('N must be odd for composite Simpson rule');
end

% Get system dimensions
K0 = K_func(0);
[n, m] = size(K0);

fprintf('Computing Gramian: n=%d, m=%d, T=%.2f, N=%d\n', n, m, T, N);

% Quadrature setup (composite Simpson)
tau = linspace(0, T, N);
w = simpson_weights(N, T);

% Initialize Gramian
W = zeros(n^2, n^2);

% Main loop over quadrature nodes
for i = 1:N
    % Get K(tau_i)
    Ki = K_func(tau(i));
    
    % Initialize propagation matrix for this node
    M_i = zeros(n^2, m*n);
    
    % Loop over input columns (k = 1 to m)
    for k = 1:m
        % Get k-th column of K(tau_i)
        kcol = Ki(:, k); % n x 1 vector
        
        % Loop over basis directions (j = 1 to n)
        for j = 1:n
            % Create standard basis vector e_j
            ej = zeros(n, 1);
            ej(j) = 1;
            
            % Initial condition: Z0 = kcol * ej^T (n x n matrix)
            Z0 = kcol * ej';
            
            % Solve Sylvester ODE: dZ/dt = A(t)Z + ZB(t)
            % from t = tau(i) to t = T
            sylv_ode = @(t, Z_vec) sylvester_rhs(t, Z_vec, A_func, B_func, n);
            
            % ODE solver options for high accuracy
            opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
            
            % Solve ODE (Z0 must be vectorized for ode45)
            [~, Z_sol] = ode45(sylv_ode, [tau(i), T], Z0(:), opts);
            
            % Extract final value and reshape back to n x n
            Z_final = reshape(Z_sol(end, :), n, n);
            
            % Store vectorized result in appropriate column of M_i
            col_idx = (k-1)*n + j;
            M_i(:, col_idx) = Z_final(:);
        end
    end
    
    % Accumulate contribution to Gramian: W += w_i * M_i * M_i^T
    W = W + w(i) * (M_i * M_i');
    
    % Progress indicator for large computations
    if mod(i, max(1, floor(N/10))) == 0
        fprintf('Progress: %d/%d nodes completed\n', i, N);
    end
end

fprintf('Gramian computation completed successfully\n');

end

function dZ_vec = sylvester_rhs(t, Z_vec, A_func, B_func, n)
% Right-hand side function for Sylvester ODE: dZ/dt = A(t)Z + ZB(t)
% 
% INPUTS:
%   t      - Current time
%   Z_vec  - Vectorized Z matrix (n^2 x 1)
%   A_func - Function handle for A(t)
%   B_func - Function handle for B(t) 
%   n      - Matrix dimension
%
% OUTPUT:
%   dZ_vec - Vectorized derivative dZ/dt

% Reshape vectorized Z back to matrix form
Z = reshape(Z_vec, n, n);

% Compute A(t) and B(t) at current time
At = A_func(t);
Bt = B_func(t);

% Sylvester equation: dZ/dt = A(t)*Z + Z*B(t)
dZ = At * Z + Z * Bt;

% Vectorize result
dZ_vec = dZ(:);

end

function w = simpson_weights(N, T)
% SIMPSON_WEIGHTS Composite Simpson quadrature weights
%
% INPUTS:
%   N - Number of nodes (must be odd)
%   T - Integration interval length [0, T]
%
% OUTPUT:
%   w - Vector of quadrature weights (1 x N)

if mod(N, 2) == 0
    error('N must be odd for composite Simpson rule');
end

% Step size
h = T / (N - 1);

% Initialize weights
w = zeros(1, N);

% Simpson's rule weights: 1, 4, 2, 4, 2, ..., 4, 1
w(1) = h/3;        % First point
w(N) = h/3;        % Last point

for i = 2:N-1
    if mod(i-1, 2) == 0  % Even index (middle of interval)
        w(i) = 4*h/3;
    else                 % Odd index
        w(i) = 2*h/3;
    end
end

end
