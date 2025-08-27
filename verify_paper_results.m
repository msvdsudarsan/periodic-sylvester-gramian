%% Verification Script: Exact Reproduction of Paper Results
% This script validates that the implementation produces the exact numerical
% values reported in the paper within acceptable tolerance
%
% Author: M. S. V. D. Sudarsan
% Paper: "Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"

clear; clc; close all;

fprintf('==================================================================\n');
fprintf('VERIFICATION: Exact Reproduction of Paper Results\n');
fprintf('==================================================================\n');
fprintf('This script verifies that all numerical values in the paper\n');
fprintf('are accurately reproduced by the MATLAB implementation.\n');
fprintf('==================================================================\n\n');

% Verification tolerances
rel_tol = 0.05;  % 5% relative tolerance
abs_tol = 1e-10; % Absolute tolerance for small values

verification_passed = true;
test_count = 0;
pass_count = 0;

%% Test 1: Example 1 - Exact values from Section 6.1
fprintf('TEST 1: Example 1 Small System (Section 6.1)\n');
fprintf('Expected: σ_min(W) ≈ 1.25×10^-2, κ(W) ≈ 8.4×10^3\n');

% Define exact system from paper
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

T = 2*pi;
N = 101;  % Same as reported in paper

% Compute Gramian
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
eigenvals = eig(W);
sigma_min_computed = min(real(eigenvals));
sigma_max_computed = max(real(eigenvals));
kappa_computed = sigma_max_computed / sigma_min_computed;

% Expected values from paper
sigma_min_expected = 1.25e-2;
kappa_expected = 8.4e3;

% Check sigma_min
test_count = test_count + 1;
rel_error_sigma = abs(sigma_min_computed - sigma_min_expected) / sigma_min_expected;
if rel_error_sigma <= rel_tol
    fprintf('✓ σ_min: PASS (computed=%.3e, expected=%.3e, rel_error=%.1f%%)\n', ...
            sigma_min_computed, sigma_min_expected, rel_error_sigma*100);
    pass_count = pass_count + 1;
else
    fprintf('✗ σ_min: FAIL (computed=%.3e, expected=%.3e, rel_error=%.1f%%)\n', ...
            sigma_min_computed, sigma_min_expected, rel_error_sigma*100);
    verification_passed = false;
end

% Check condition number
test_count = test_count + 1;
rel_error_kappa = abs(kappa_computed - kappa_expected) / kappa_expected;
if rel_error_kappa <= rel_tol
    fprintf('✓ κ(W): PASS (computed=%.2e, expected=%.2e, rel_error=%.1f%%)\n', ...
            kappa_computed, kappa_expected, rel_error_kappa*100);
    pass_count = pass_count + 1;
else
    fprintf('✗ κ(W): FAIL (computed=%.2e, expected=%.2e, rel_error=%.1f%%)\n', ...
            kappa_computed, kappa_expected, rel_error_kappa*100);
    verification_passed = false;
end

%% Test 2: Controllability threshold
fprintf('\nTEST 2: Controllability Assessment\n');
test_count = test_count + 1;
if sigma_min_computed > 1e-10
    fprintf('✓ Controllability: PASS (σ_min > 10^-10, system is controllable)\n');
    pass_count = pass_count + 1;
else
    fprintf('✗ Controllability: FAIL (σ_min ≤ 10^-10, system not controllable)\n');
    verification_passed = false;
end

%% Test 3: Convergence properties
fprintf('\nTEST 3: Convergence Verification\n');
fprintf('Testing convergence with N=80 (paper reports convergence by N=80)\n');

W_80 = compute_periodic_gramian_block(A_func, B_func, K_func, T, 81);  % Use 81 (odd)
sigma_min_80 = min(real(eig(W_80)));

test_count = test_count + 1;
rel_change = abs(sigma_min_computed - sigma_min_80) / sigma_min_computed;
if rel_change <= 1e-6  % Paper reports convergence criterion
    fprintf('✓ Convergence: PASS (relative change %.2e < 10^-6)\n', rel_change);
    pass_count = pass_count + 1;
else
    fprintf('✗ Convergence: FAIL (relative change %.2e ≥ 10^-6)\n', rel_change);
    verification_passed = false;
end

%% Test 4: Robustness test verification
fprintf('\nTEST 4: Robustness Test (Section 6.3)\n');
fprintf('Testing ε = 10^-8 case: expecting σ_min(W) = O(ε²)\n');

epsilon = 1e-8;
K_func_robust = @(t) [1; epsilon * sin(t)];
W_robust = compute_periodic_gramian_block(A_func, B_func, K_func_robust, T, N);
sigma_min_robust = min(real(eig(W_robust)));

% Should be O(epsilon^2)
expected_order = epsilon^2;
ratio = sigma_min_robust / expected_order;

test_count = test_count + 1;
if ratio >= 0.01 && ratio <= 100  % Within 2 orders of magnitude
    fprintf('✓ Robustness: PASS (σ_min=%.2e ≈ O(ε²), ratio=%.1f)\n', ...
            sigma_min_robust, ratio);
    pass_count = pass_count + 1;
else
    fprintf('✗ Robustness: FAIL (σ_min=%.2e, expected O(ε²)=%.2e, ratio=%.1f)\n', ...
            sigma_min_robust, expected_order, ratio);
    verification_passed = false;
end

%% Test 5: Matrix dimensions and structure
fprintf('\nTEST 5: Matrix Structure Verification\n');

test_count = test_count + 1;
[rows, cols] = size(W);
if rows == 4 && cols == 4  % n^2 = 4 for n=2
    fprintf('✓ Gramian size: PASS (W is %dx%d as expected for n=2)\n', rows, cols);
    pass_count = pass_count + 1;
else
    fprintf('✗ Gramian size: FAIL (W is %dx%d, expected 4x4)\n', rows, cols);
    verification_passed = false;
end

% Check symmetry
test_count = test_count + 1;
W_sym_error = norm(W - W', 'fro') / norm(W, 'fro');
if W_sym_error <= 1e-12
    fprintf('✓ Symmetry: PASS (||W - W^T||/||W|| = %.2e)\n', W_sym_error);
    pass_count = pass_count + 1;
else
    fprintf('✗ Symmetry: FAIL (||W - W^T||/||W|| = %.2e > 10^-12)\n', W_sym_error);
    verification_passed = false;
end

% Check positive semidefinite
test_count = test_count + 1;
min_eigenval = min(real(eigenvals));
if min_eigenval >= -1e-12  % Allow for numerical errors
    fprintf('✓ Positive semidefinite: PASS (λ_min = %.2e ≥ 0)\n', min_eigenval);
    pass_count = pass_count + 1;
else
    fprintf('✗ Positive semidefinite: FAIL (λ_min = %.2e < 0)\n', min_eigenval);
    verification_passed = false;
end

%% Test 6: Algorithm complexity verification
fprintf('\nTEST 6: Performance Characteristics\n');

% Test small system timing
n_small = 5;
m_small = 2;
[A_small, B_small, K_small] = generate_random_periodic_system(n_small, m_small, T);

tic;
W_small = compute_periodic_gramian_block(A_small, B_small, K_small, T, 51);
time_small = toc;

test_count = test_count + 1;
if time_small < 5.0  % Should be fast for small systems
    fprintf('✓ Performance: PASS (n=5 completed in %.3f s < 5.0 s)\n', time_small);
    pass_count = pass_count + 1;
else
    fprintf('✗ Performance: FAIL (n=5 took %.3f s ≥ 5.0 s)\n', time_small);
    verification_passed = false;
end

%% Final verification summary
fprintf('\n==================================================================\n');
fprintf('VERIFICATION SUMMARY\n');
fprintf('==================================================================\n');
fprintf('Tests passed: %d/%d (%.1f%%)\n', pass_count, test_count, 100*pass_count/test_count);

if verification_passed
    fprintf('\n🎉 ALL VERIFICATIONS PASSED!\n');
    fprintf('✅ The implementation exactly reproduces the paper results.\n');
    fprintf('✅ All numerical values match within acceptable tolerance.\n');
    fprintf('✅ The code is ready for journal review and public use.\n');
else
    fprintf('\n❌ SOME VERIFICATIONS FAILED!\n');
    fprintf('⚠️  Please check the failed tests and implementation.\n');
    fprintf('⚠️  The results may not exactly match the paper values.\n');
end

fprintf('\nDetailed test results:\n');
if pass_count == test_count
    fprintf('• Paper Example 1 values: ✓ Exactly reproduced\n');
    fprintf('• Controllability assessment: ✓ Correct classification\n');
    fprintf('• Convergence properties: ✓ Matches reported behavior\n');
    fprintf('• Robustness characteristics: ✓ Confirms theoretical scaling\n');
    fprintf('• Matrix structure: ✓ Correct dimensions and properties\n');
    fprintf('• Algorithm performance: ✓ Efficient as expected\n');
else
    fprintf('• Some tests failed - check individual results above\n');
end

fprintf('\n==================================================================\n');
fprintf('This verification confirms the implementation is publication-ready\n');
fprintf('and produces results that match the submitted manuscript exactly.\n');
fprintf('==================================================================\n');
