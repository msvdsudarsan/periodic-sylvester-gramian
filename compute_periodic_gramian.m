

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
% Email:  msvdsudarsan@gmail.com
% Date:   August 2025

    % -------------------------------
    % Basic input validation
    % -------------------------------
    if nargin < 5
        error('All five inputs are required');
    end

    % Infer dimensions from K(0)
    K0 = K_func(0);
    [n, m] = size(K0);

    % -------------------------------
    % Quadrature setup (Composite Simpson)
    % -------------------------------
    if mod(N, 2) == 0
        N = N + 1; % ensure odd N for Simpson rule
        warning('N adjusted to %d for Simpson''s rule', N);
    end
    tau = linspace(0, T, N);
    w   = simpson_weights(N, T);

    % -------------------------------
    % Initialize Gramian
    % -------------------------------
    W = zeros(n^2, n^2);

    fprintf('Computing Gramian with %d quadrature nodes...\n', N);

    % -------------------------------
    % Main loop over quadrature nodes
    % -------------------------------
    for i = 1:N
        if mod(i, 20) == 0
            fprintf('Progress: %d/%d nodes processed\n', i, N);
        end

        % Preallocate M_i (n^2 x (m*n)) so that:
        % columns  (k-1)*n + j  store vec(Z_final^(k,j))
        % where Z_final^(k,j) solves dZ/dt = A Z + Z B with Z(tau_i)=K(:,k) e_j^T
        M_i = zeros(n^2, m * n);

        % Evaluate K at current quadrature node once
        K_at_tau = K_func(tau(i));

        % Loop over input columns k and basis columns j
        for k = 1:m
            Kk = K_at_tau(:, k); % n x 1

            for j = 1:n
                % Unit vector e_j
                ej = zeros(n,1); ej(j) = 1;

                % Correct initial condition for the j-th column:
                % Z0 = K(:,k) * e_j^T  (n x n; only j-th column equals K(:,k))
                Z0 = Kk * (ej.');

                % Propagate Sylvester ODE from t = tau(i) to t = T
                Z_final = solve_sylvester_ode(Z0, tau(i), T, A_func, B_func);

                % Store vec(Z_final) in the appropriate column of M_i
                col_idx = (k-1)*n + j;
                M_i(:, col_idx) = Z_final(:);
            end
        end

        % Accumulate Gramian contribution: W += w_i * M_i * M_i^T
        W = W + w(i) * (M_i * M_i.');
    end

    fprintf('Gramian computation completed.\n');

    % -------------------------------
    % Report basic spectral diagnostics
    % -------------------------------
    % Enforce symmetry for numerical safety
    Ws = 0.5 * (W + W.');
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

% -------------------------------
% Helpers
% -------------------------------
function Z_final = solve_sylvester_ode(Z0, t0, tf, A_func, B_func)
    % Solve dZ/dt = A(t)*Z + Z*B(t) from t0 to tf with Z(t0) = Z0
    n = size(Z0, 1);
    sylv_ode = @(t, z_vec) sylvester_rhs_vec(t, z_vec, A_func, B_func, n);
    z0_vec   = Z0(:);
    options  = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
    [~, z_sol] = ode45(sylv_ode, [t0, tf], z0_vec, options);
    Z_final = reshape(z_sol(end, :).', n, n);
end

function dz_dt = sylvester_rhs_vec(t, z_vec, A_func, B_func, n)
    % Right-hand side for vectorized Sylvester equation
    Z     = reshape(z_vec, n, n);
    dZ_dt = A_func(t) * Z + Z * B_func(t);
    dz_dt = dZ_dt(:);
end

