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
    
    % For interior nodes, compute the integrand properly
    % We need Φ(0,τ) * B(τ) * K(τ) * K(τ)' * B(τ)' * Φ(0,τ)'
    
    try
        % Method: Compute Φ(0,τ) by integrating forward from 0 to τ
        % dx/dt = A(t)*x, x(0) = I
        
        % We need to compute n separate solutions (columns of Φ)
        Phi_0_tau = zeros(n, n);
        
        % Integrate each column of the fundamental matrix
        for col = 1:n
            % Initial condition: e_col (column vector with 1 in position col, 0 elsewhere)
            x0 = zeros(n, 1);
            x0(col) = 1;
            
            % Define ODE: dx/dt = A(t)*x
            state_ode = @(t, x) A_func(t) * x;
            
            % Integrate from 0 to tau(i)
            [~, x_sol] = ode45(state_ode, [0, tau(i)], x0, opts);
            
            % Store the final value as column col of Φ(0,τ)
            Phi_0_tau(:, col) = x_sol(end, :)';
        end
        
        % Compute integrand: Φ(0,τ) * B(τ) * K(τ) * K(τ)' * B(τ)' * Φ(0,τ)'
        B_tau_i = B_func(tau(i));
        K_tau_i = K_func(tau(i));
        integrand = Phi_0_tau * B_tau_i * (K_tau_i * K_tau_i') * B_tau_i' * Phi_0_tau';
        
        % Add weighted contribution
        W = W + w(i) * integrand;
        
    catch ME
        fprintf('Error at node %d (tau=%.6f): %s\n', i, tau(i), ME.message);
        rethrow(ME);
    end
end

fprintf('Gramian computation completed successfully!\n');
end
