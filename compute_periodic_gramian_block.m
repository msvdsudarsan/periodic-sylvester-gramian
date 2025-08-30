function W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
%COMPUTE_PERIODIC_GRAMIAN_BLOCK Block-wise computation of reachability Gramian
%
% Syntax:
%   W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
%
% Inputs:
%   A_func - Function handle for A(t), returns n x n matrix
%   B_func - Function handle for B(t), returns n x n matrix  
%   K_func - Function handle for K(t), returns n x m matrix
%   T      - Period (scalar)
%   N      - Number of quadrature nodes (must be odd for Simpson rule)
%
% Output:
%   W      - Reachability Gramian (n^2 x n^2 matrix)
%
% Description:
%   Computes the reachability Gramian for the periodic Sylvester system:
%   dX/dt = A(t)X + XB(t) + K(t)U(t)
%   
%   Uses block-wise propagation to avoid forming n^2 x n^2 Kronecker matrices.
%   Computational complexity: O(N*n^3*m) instead of O(N*n^6).
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com
% Date: 2025

% Input validation
if mod(N,2) == 0
    error('N must be odd for composite Simpson rule');
end

% Get dimensions from K(0)
K0 = K_func(0);
[n, m] = size(K0);

% Quadrature setup (composite Simpson rule)
tau = linspace(0, T, N);
w = simpson_weights(N, T);

% Initialize Gramian
W = zeros(n^2, n^2);

% Main loop over quadrature nodes
for i = 1:N
    Ki = K_func(tau(i));
    M_i = zeros(n^2, m*n);  % Matrix to store propagated vectors

    % Loop over input columns
    for k = 1:m
        zcol = Ki(:, k);  % k-th column of K(tau(i))
        
        % Loop over basis columns  
        for j = 1:n
            % Create rank-1 initial condition: Z0 = zcol * e_j^T
            ej = zeros(n, 1); 
            ej(j) = 1;
            Z0 = zcol * (ej.');

            % Handle endpoint case (no propagation needed)
            if abs(tau(i) - T) < 1e-10
                Z_final = Z0;
            else
                % Solve Sylvester ODE: dZ/dt = A(t)Z + ZB(t)
                sylv_ode = @(t, Z_vec) sylvester_rhs(t, Z_vec, A_func, B_func);
                opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
                [~, Z_sol] = ode45(sylv_ode, [tau(i), T], Z0(:), opts);
                
                % Reshape final solution
                Z_final = reshape(Z_sol(end, :), n, n);
            end

            % Store vectorized result
            col_idx = (k-1)*n + j;
            M_i(:, col_idx) = Z_final(:);
        end
    end

    % Accumulate Gramian contribution
    W = W + w(i) * (M_i * M_i');
end

end

function dZ_vec = sylvester_rhs(t, Z_vec, A_func, B_func)
%SYLVESTER_RHS Right-hand side for vectorized Sylvester ODE
%
% Computes d/dt vec(Z) = vec(A(t)Z + ZB(t))

n = round(sqrt(length(Z_vec)));
Z = reshape(Z_vec, n, n);

% Compute dZ/dt = A(t)Z + ZB(t)
dZ = A_func(t) * Z + Z * B_func(t);

% Return vectorized form
dZ_vec = dZ(:);
end

function w = simpson_weights(N, T)
%SIMPSON_WEIGHTS Composite Simpson quadrature weights
%
% Inputs:
%   N - Number of nodes (must be odd)
%   T - Integration interval [0, T]
%
% Output:
%   w - Weight vector (1 x N)

if mod(N,2) == 0
    error('N must be odd for composite Simpson rule');
end

h = T / (N - 1);
w = zeros(1, N);

% Simpson's 1/3 rule weights
w(1) = h/3;
w(N) = h/3;

for i = 2:N-1
    if mod(i-1, 2) == 0
        w(i) = 4*h/3;  % Even index (middle points)
    else
        w(i) = 2*h/3;  % Odd index
    end
end
end
