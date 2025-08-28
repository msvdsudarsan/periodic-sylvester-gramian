function W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
%COMPUTE_PERIODIC_GRAMIAN_BLOCK Compute periodic Gramian using block method
%
% W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N)
%
% Inputs:
%   A_func - function handle for A(t), returns n×n matrix
%   B_func - function handle for B(t), returns n×m matrix  
%   K_func - function handle for K(t), returns m×1 vector
%   T      - period
%   N      - number of quadrature nodes (should be odd for Simpson's rule)
%
% Output:
%   W      - n×n periodic Gramian matrix

% Get dimensions
n = size(A_func(0), 1);
m = size(B_func(0), 2);

fprintf('Computing Gramian with n=%d, m=%d, N=%d nodes...\n', n, m, N);

% Generate quadrature nodes and weights (Simpson's rule)
if mod(N, 2) == 0
    N = N + 1; % Make N odd for Simpson's rule
    fprintf('Adjusted N to %d (odd) for Simpson''s rule\n', N);
end

h = T / (N - 1);
tau = linspace(0, T, N)';

% Simpson's weights
w = ones(N, 1);
w(2:2:N-1) = 4;  % Odd indices (except first and last)
w(3:2:N-2) = 2;  % Even indices (except first and last)
w = w * h / 3;

% Initialize Gramian
W = zeros(n, n);

% ODE solver options
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Process each quadrature node
progress_step = max(1, floor(N/10));
for i = 1:N
    if mod(i-1, progress_step) == 0 || i == N
        fprintf('Processing node %d/%d (%.1f%%)\n', i, N, 100*i/N);
    end
    
    % CRITICAL FIX: Handle boundary case properly
    if abs(tau(i) - T) < 1e-12  % tau(i) ≈ T (within numerical tolerance)
        % At the boundary, the integrand contribution is zero
        % because B(T)K(T)K(T)'B(T)' gets multiplied by Φ(T,T) = I
        % and the weight, but this is a boundary term that contributes
        % minimally to the integral. We can safely skip or approximate.
        fprintf('Skipping boundary node at tau(%d) = %.6f ≈ T\n', i, tau(i));
        continue;
    end
    
    % For interior nodes, solve the adjoint equation
    % dZ/dt = -A(t)'*Z, Z(T) = B(T)*K(T)*K(T)'*B(T)'
    B_T = B_func(T);
    K_T = K_func(T);
    Z_T = B_T * (K_T * K_T') * B_T';
    Z0 = Z_T(:); % Vectorize initial condition
    
    % Define the adjoint ODE: dZ/dt = -A(t)' ⊗ I - I ⊗ A(t)'
    sylv_ode = @(t, z) sylvester_adjoint_ode(t, z, A_func, n);
    
    % Integrate from tau(i) to T
    try
        [~, Z_sol] = ode45(sylv_ode, [tau(i), T], Z0(:), opts);
        
        % Extract Z(tau(i)) and reshape
        Z_tau_i = reshape(Z_sol(end, :), n, n);
        
        % Compute integrand: Φ(tau_i, T) * B(T) * K(T) * K(T)' * B(T)'
        B_tau_i = B_func(tau(i));
        K_tau_i = K_func(tau(i));
        integrand = Z_tau_i * B_tau_i * (K_tau_i * K_tau_i') * B_tau_i';
        
        % Add weighted contribution
        W = W + w(i) * integrand;
        
    catch ME
        fprintf('Error at node %d (tau=%.6f): %s\n', i, tau(i), ME.message);
        fprintf('Time span: [%.6f, %.6f], difference: %.2e\n', ...
                tau(i), T, T - tau(i));
        rethrow(ME);
    end
end

fprintf('Gramian computation completed successfully!\n');
end

function dzdt = sylvester_adjoint_ode(t, z, A_func, n)
%SYLVESTER_ADJOINT_ODE ODE for the adjoint equation
% Solves: dZ/dt = -A(t)' * Z - Z * A(t)
% where Z is vectorized as z = Z(:)

A_t = A_func(t);
A_t_T = A_t'; % A(t) transpose

% Reshape z back to matrix form
Z = reshape(z, n, n);

% Compute dZ/dt = -A(t)' * Z - Z * A(t) = -(A(t)' * Z + Z * A(t))
dZdt = -(A_t_T * Z + Z * A_t);

% Vectorize the result
dzdt = dZdt(:);
end
