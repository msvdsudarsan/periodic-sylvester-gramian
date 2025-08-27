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

    if nargin < 5
        error('All five inputs are required');
    end

    K0 = K_func(0);
    [n, m] = size(K0);

    if mod(N, 2) == 0
        N = N + 1;
        warning('N adjusted to %d for Simpson''s rule', N);
    end
    tau = linspace(0, T, N);
    w   = simpson_weights(N, T);

    W = zeros(n^2, n^2);
    fprintf('Computing Gramian with %d quadrature nodes...\n', N);

    for i = 1:N
        if mod(i, 20) == 0
            fprintf('Progress: %d/%d nodes processed\n', i, N);
        end
        K_at_tau = K_func(tau(i));
        M_i = block_sylvester_propagate(A_func, B_func, K_at_tau, tau(i), T, ...
                                        'RelTol',1e-9,'AbsTol',1e-12,'Solver','ode45');
        W = W + w(i) * (M_i * M_i.');
    end

    fprintf('Gramian computation completed.\n');

    % Symmetrize before diagnostics
    Ws = 0.5*(W+W.');
    ev = eig(Ws);
    sigma_min = min(ev);
    sigma_max = max(ev);
    condition_number = sigma_max / max(sigma_min, eps);

    fprintf('Gramian properties:\n');
    fprintf('  Minimum eigenvalue: %.6e\n', sigma_min);
    fprintf('  Condition number:   %.2e\n', condition_number);

    if sigma_min > 1e-12
        fprintf('  System appears to be controllable.\n');
    else
        fprintf('  System may not be controllable (small min eigenvalue).\n');
    end
end
