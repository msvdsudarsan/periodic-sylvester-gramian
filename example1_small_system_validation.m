function example1_small_system_validation()
%EXAMPLE1_SMALL_SYSTEM_VALIDATION Validation of small system from paper
%
% Reproduces Example 1 from the research paper with corrected parameters.
% System: n=2, m=1, T=2π
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com
% Date: 2025

clc;
fprintf('=== EXAMPLE 1: SMALL SYSTEM VALIDATION ===\n');
fprintf('Reproducing corrected results from research paper\n\n');

% System parameters
n = 2; 
m = 1; 
T = 2*pi;
N = 101;  % Quadrature nodes

% Define system matrices exactly as in paper
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];

% CORRECTED K(t) with proper scaling (0.079 instead of 14.958)
K_func = @(t) 0.079 * [1 + 0.2*cos(t); 0.5*sin(t)];

fprintf('System matrices:\n');
fprintf('A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]\n');
fprintf('B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)]\n');
fprintf('K(t) = 0.079 * [1 + 0.2*cos(t); 0.5*sin(t)]  (CORRECTED)\n\n');

% Display K(t) at key points
fprintf('K(t) values:\n');
fprintf('K(0)   = [%.4f; %.4f]\n', K_func(0));
fprintf('K(π/2) = [%.4f; %.4f]\n', K_func(pi/2));
fprintf('K(π)   = [%.4f; %.4f]\n', K_func(pi));
fprintf('K(3π/2) = [%.4f; %.4f]\n\n', K_func(3*pi/2));

% Compute reachability Gramian
fprintf('Computing reachability Gramian (N=%d nodes)...\n', N);
tic;
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
elapsed_time = toc;

% Analyze Gramian
sigma_vals = svd(W);
sigma_min = min(sigma_vals);
sigma_max = max(sigma_vals);
kappa_W = sigma_max / sigma_min;
rank_W = rank(W, 1e-12);

% Display results
fprintf('\n--- NUMERICAL RESULTS ---\n');
fprintf('Elapsed time: %.4f seconds\n', elapsed_time);
fprintf('Gramian size: %dx%d\n', size(W));
fprintf('Rank of W:   %d (full rank = %d)\n', rank_W, n^2);
fprintf('σ_min(W) = %.6e\n', sigma_min);
fprintf('σ_max(W) = %.6e\n', sigma_max);
fprintf('κ(W)     = %.6e\n', kappa_W);

% Controllability assessment
if sigma_min > 1e-10
    fprintf('✓ System is CONTROLLABLE (σ_min > 0)\n');
else
    fprintf('✗ System is NOT CONTROLLABLE (σ_min ≈ 0)\n');
end

% Check conditioning
if kappa_W < 1e6
    fprintf('✓ Gramian is well-conditioned\n');
elseif kappa_W < 1e10
    fprintf('! Gramian is moderately conditioned\n');
else
    fprintf('✗ Gramian is ill-conditioned\n');
end

% Display full Gramian (for small system)
fprintf('\n--- FULL GRAMIAN MATRIX ---\n');
fprintf('W = \n');
disp(W);

% Eigenvalue analysis
eig_W = eig(W);
fprintf('Eigenvalues of W:\n');
for i = 1:length(eig_W)
    fprintf('  λ_%d = %.6e\n', i, eig_W(i));
end

% Compare with paper values
fprintf('\n--- COMPARISON WITH PAPER ---\n');
expected_sigma_min = 1.071e-02;  % Corrected value
expected_kappa = 2.761e+00;     % Corrected value

sigma_error = abs(sigma_min - expected_sigma_min) / expected_sigma_min;
kappa_error = abs(kappa_W - expected_kappa) / expected_kappa;

fprintf('Expected: σ_min = %.6e, κ = %.6e\n', expected_sigma_min, expected_kappa);
fprintf('Computed: σ_min = %.6e, κ = %.6e\n', sigma_min, kappa_W);
fprintf('Relative errors: σ_min = %.2e, κ = %.2e\n', sigma_error, kappa_error);

if sigma_error < 0.05 && kappa_error < 0.05
    fprintf('✓ EXCELLENT AGREEMENT with corrected paper values!\n');
elseif sigma_error < 0.1 && kappa_error < 0.1
    fprintf('✓ Good agreement with corrected paper values\n');
else
    fprintf('! Some discrepancy with paper values\n');
end

fprintf('\n=== EXAMPLE 1 COMPLETE ===\n');

end
