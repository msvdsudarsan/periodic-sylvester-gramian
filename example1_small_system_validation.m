%% Example 1: Small System Validation (n=2, m=1, T=2π)
% This reproduces the exact example from Section 6.1 of the paper
%
% Author: M. S. V. D. Sudarsan
% Paper: "Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"

clear; clc; close all;

fprintf('=== EXAMPLE 1: Small System Validation ===\n');
fprintf('System parameters: n=2, m=1, T=2π\n\n');

% System parameters
n = 2;
m = 1;
T = 2*pi;
N = 101;  % Number of quadrature nodes (must be odd)

% Define periodic coefficient functions exactly as in the paper
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];

B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];

K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

fprintf('Coefficient matrices:\n');
fprintf('A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]\n');
fprintf('B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)]\n');
fprintf('K(t) = [1 + 0.2*cos(t); 0.5*sin(t)]\n\n');

% Compute reachability Gramian using block method
fprintf('Computing reachability Gramian using block method...\n');
tic;
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
comp_time = toc;
fprintf('Computation time: %.4f seconds\n\n', comp_time);

% Analyze results
fprintf('=== GRAMIAN ANALYSIS ===\n');
fprintf('Gramian matrix W (4x4):\n');
disp(W);

% Compute eigenvalues
lambda = eig(W);
lambda_real = real(lambda);
sigma_min = min(lambda_real);
sigma_max = max(lambda_real);
cond_num = sigma_max / sigma_min;

fprintf('Eigenvalues of W:\n');
for i = 1:length(lambda)
    if abs(imag(lambda(i))) < 1e-12
        fprintf('  λ_%d = %.6e\n', i, real(lambda(i)));
    else
        fprintf('  λ_%d = %.6e + %.6ei\n', i, real(lambda(i)), imag(lambda(i)));
    end
end

fprintf('\nGramian properties:\n');
fprintf('  σ_min(W) ≈ %.2e\n', sigma_min);
fprintf('  σ_max(W) ≈ %.2e\n', sigma_max);
fprintf('  κ(W) ≈ %.1e\n', cond_num);

% Controllability assessment
fprintf('\n=== CONTROLLABILITY ASSESSMENT ===\n');
if sigma_min > 1e-10
    fprintf('✓ System is CONTROLLABLE (σ_min > 0)\n');
    fprintf('  The system can be steered between arbitrary states.\n');
else
    fprintf('✗ System is NOT controllable (σ_min ≈ 0)\n');
    fprintf('  The system cannot reach all states in the state space.\n');
end

% Compare with paper results
fprintf('\n=== COMPARISON WITH PAPER RESULTS ===\n');
fprintf('Paper reports:\n');
fprintf('  σ_min(W) ≈ 1.25×10^-2\n');
fprintf('  κ(W) ≈ 8.4×10^3\n');
fprintf('Computed results:\n');
fprintf('  σ_min(W) ≈ %.2e\n', sigma_min);
fprintf('  κ(W) ≈ %.1e\n', cond_num);

% Relative errors
if abs(1.25e-2) > eps
    rel_err_sigma = abs(sigma_min - 1.25e-2) / 1.25e-2;
    fprintf('  Relative error in σ_min: %.2e\n', rel_err_sigma);
end

if abs(8.4e3) > eps
    rel_err_cond = abs(cond_num - 8.4e3) / 8.4e3;
    fprintf('  Relative error in κ(W): %.2e\n', rel_err_cond);
end

% Minimum energy control example
fprintf('\n=== MINIMUM ENERGY CONTROL EXAMPLE ===\n');
if sigma_min > 1e-10
    % Define initial and final states
    x0_mat = [1, 0; 0, 1];  % Initial state matrix
    xf_mat = [0, 1; -1, 0]; % Final state matrix
    
    x0 = x0_mat(:);  % Vectorize
    xf = xf_mat(:);  % Vectorize
    
    fprintf('Initial state matrix X(0):\n');
    disp(x0_mat);
    fprintf('Final state matrix X(T):\n');
    disp(xf_mat);
    
    % Compute state transition matrix Φ(T,0)
    Phi_T0 = compute_state_transition(A_func, B_func, T, n);
    
    % Compute minimum energy control coefficient
    target = xf - Phi_T0 * x0;
    control_coeff = W \ target;
    
    fprintf('Control coefficient vector:\n');
    disp(control_coeff);
    
    % Compute control energy
    control_energy = control_coeff' * W * control_coeff;
    fprintf('Minimum control energy: %.6f\n', control_energy);
end

fprintf('\n=== EXAMPLE 1 COMPLETED ===\n');

function Phi = compute_state_transition(A_func, B_func, T, n)
% Compute state transition matrix Φ(T,0) for vectorized system
% This computes the fundamental solution of dX/dt = A(t)X + XB(t)

% Initialize identity matrix
I_n2 = eye(n^2);
Phi = zeros(n^2, n^2);

% Solve n^2 initial value problems
for k = 1:n^2
    % k-th unit vector as initial condition
    x0 = I_n2(:, k);
    X0 = reshape(x0, n, n);
    
    % Solve Sylvester ODE
    sylv_ode = @(t, Z_vec) sylvester_rhs_homogeneous(t, Z_vec, A_func, B_func, n);
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
    [~, X_sol] = ode45(sylv_ode, [0, T], X0(:), opts);
    
    % Extract final value
    Phi(:, k) = X_sol(end, :)';
end
end

function dZ_vec = sylvester_rhs_homogeneous(t, Z_vec, A_func, B_func, n)
% RHS for homogeneous Sylvester equation dZ/dt = A(t)*Z + Z*B(t)
Z = reshape(Z_vec, n, n);
dZ = A_func(t) * Z + Z * B_func(t);
dZ_vec = dZ(:);
end
