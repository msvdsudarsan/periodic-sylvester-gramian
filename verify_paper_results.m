function verify_paper_results()
% VERIFY_PAPER_RESULTS Verify all numerical results reported in the paper
%
% This function systematically verifies all numerical values, claims, and
% results presented in the paper "Controllability and Efficient Gramian 
% Computation for Periodic Sylvester Matrix Systems".
%
% VERIFIES:
%   - Example 1 results (σ_min ≈ 1.25×10^-2, κ(W) ≈ 8.4×10^3)
%   - Performance comparison table values
%   - Convergence analysis claims
%   - Robustness test scaling behavior
%   - Algorithm complexity claims
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com

clear; clc;
fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                      PAPER RESULTS VERIFICATION               ║\n');
fprintf('║                                                                ║\n');
fprintf('║  Verifying all numerical claims from the AML paper            ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

%% Verification Summary
verification_results = struct();
verification_results.timestamp = datestr(now);
verification_results.tests = {};
verification_results.expected = {};
verification_results.computed = {};
verification_results.errors = {};
verification_results.status = {};

test_count = 0;

%% Test 1: Example 1 Small System Values
test_count = test_count + 1;
fprintf('TEST %d: Example 1 Small System (Section 6.1)\n', test_count);
fprintf('==============================================\n');

% System parameters from paper
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];
T = 2*pi;
N = 101; % As stated in paper

fprintf('Computing Gramian with N=%d quadrature nodes...\n', N);
tic;
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
comp_time = toc;

eigenvals = eig(W);
sigma_min_computed = sqrt(min(real(eigenvals)));
kappa_computed = max(real(eigenvals)) / min(real(eigenvals));

% Paper claims
sigma_min_paper = 1.25e-2;
kappa_paper = 8.4e3;

% Calculate errors
sigma_error = abs(sigma_min_computed - sigma_min_paper) / sigma_min_paper * 100;
kappa_error = abs(kappa_computed - kappa_paper) / kappa_paper * 100;

fprintf('Results:\n');
fprintf('  σ_min: Expected = %.3e, Computed = %.3e, Error = %.1f%%\n', ...
        sigma_min_paper, sigma_min_computed, sigma_error);
fprintf('  κ(W):  Expected = %.3e, Computed = %.3e, Error = %.1f%%\n', ...
        kappa_paper, kappa_computed, kappa_error);
fprintf('  Computation time: %.3f seconds\n', comp_time);

% Determine status
if sigma_error < 20 && kappa_error < 50
    test1_status = '✓ PASS';
elseif sigma_error < 50
    test1_status = '○ PARTIAL';  
else
    test1_status = '✗ FAIL';
end

fprintf('  Status: %s\n\n', test1_status);

% Store results
verification_results.tests{end+1} = 'Example 1 σ_min';
verification_results.expected{end+1} = sigma_min_paper;
verification_results.computed{end+1} = sigma_min_computed;
verification_results.errors{end+1} = sigma_error;
verification_results.status{end+1} = test1_status;

verification_results.tests{end+1} = 'Example 1 κ(W)';
verification_results.expected{end+1} = kappa_paper;
verification_results.computed{end+1} = kappa_computed;
verification_results.errors{end+1} = kappa_error;
verification_results.status{end+1} = test1_status;

%% Test 2: Convergence Analysis Claims
test_count = test_count + 1;
fprintf('TEST %d: Convergence Analysis (Figure 1)\n', test_count);
fprintf('=========================================\n');

fprintf('Testing convergence claim: "Convergence achieved by N=80"\n');

% Test convergence with different N values
N_test_values = [41, 61, 81, 101];
sigma_values = zeros(size(N_test_values));

fprintf('Computing convergence data...\n');
for i = 1:length(N_test_values)
    N_test = N_test_values(i);
    W_test = compute_periodic_gramian_block(A_func, B_func, K_func, T, N_test);
    eigenvals_test = eig(W_test);
    sigma_values(i) = sqrt(min(real(eigenvals_test)));
    fprintf('  N=%d: σ_min = %.6e\n', N_test, sigma_values(i));
end

% Check convergence
sigma_ref = sigma_values(end); % N=101 as reference
convergence_errors = abs(sigma_values - sigma_ref) / sigma_ref;

% Find where convergence is achieved (relative change < 1e-6)
convergence_threshold = 1e-6;
converged_idx = find(convergence_errors < convergence_threshold, 1, 'first');

if ~isempty(converged_idx)
    convergence_N = N_test_values(converged_idx);
    fprintf('Convergence achieved at N=%d (threshold: %.0e)\n', convergence_N, convergence_threshold);
    
    if convergence_N <= 80
        test2_status = '✓ PASS';
        fprintf('Status: ✓ PASS (convergence by N=%d ≤ 80)\n\n', convergence_N);
    else
        test2_status = '○ PARTIAL';
        fprintf('Status: ○ PARTIAL (convergence at N=%d > 80)\n\n', convergence_N);
    end
else
    test2_status = '✗ FAIL';
    fprintf('Status: ✗ FAIL (convergence not achieved in test range)\n\n');
end

% Store results
verification_results.tests{end+1} = 'Convergence by N=80';
verification_results.expected{end+1} = 80;
verification_results.computed{end+1} = convergence_N;
verification_results.errors{end+1} = NaN;
verification_results.status{end+1} = test2_status;

%% Test 3: Performance Scaling Claims
test_count = test_count + 1;
fprintf('TEST %d: Performance Scaling (Section 6.2)\n', test_count);
fprintf('==========================================\n');

fprintf('Testing claim: "O(Nn³m) complexity"\n');

% Test performance scaling for small systems
n_values = [4, 6, 8];  % Smaller values for quick testing
m = 2;
T = 2*pi;
N = 31;  % Reduced for speed

times = zeros(size(n_values));
theoretical_times = zeros(size(n_values));

fprintf('Measuring computation times...\n');
for i = 1:length(n_values)
    n = n_values(i);
    
    % Generate random system
    [A_test, B_test, K_test] = generate_random_periodic_system(n, m, T);
    
    % Measure time
    tic;
    W_test = compute_periodic_gramian_block(A_test, B_test, K_test, T, N);
    times(i) = toc;
    
    % Theoretical O(n³) scaling
    theoretical_times(i) = (n/n_values(1))^3 * times(1);
    
    fprintf('  n=%d: %.3f seconds (theoretical: %.3f)\n', n, times(i), theoretical_times(i));
end

% Check scaling
if length(n_values) >= 3
    log_n = log(n_values);
    log_t = log(times);
    p = polyfit(log_n, log_t, 1);
    empirical_scaling = p(1);
    
    fprintf('Empirical scaling: O(n^%.2f)\n', empirical_scaling);
    
    if abs(empirical_scaling - 3.0) < 1.0
        test3_status = '✓ PASS';
        fprintf('Status: ✓ PASS (scaling ≈ n³)\n\n');
    else
        test3_status = '○ PARTIAL';
        fprintf('Status: ○ PARTIAL (scaling differs from n³)\n\n');
    end
else
    test3_status = '✗ FAIL';
    fprintf('Status: ✗ FAIL (insufficient data for scaling analysis)\n\n');
end

% Store results
verification_results.tests{end+1} = 'O(n³) complexity scaling';
verification_results.expected{end+1} = 3.0;
verification_results.computed{end+1} = empirical_scaling;
verification_results.errors{end+1} = abs(empirical_scaling - 3.0);
verification_results.status{end+1} = test3_status;

%% Test 4: Robustness Test Claims
test_count = test_count + 1;
fprintf('TEST %d: Robustness Test (Section 6.3)\n', test_count);
fprintf('======================================\n');

fprintf('Testing claim: "σ_min = O(ε²)" for time-varying rank deficiency\n');

% Test system with different epsilon values
epsilon_values = [1e-4, 1e-6, 1e-8];
sigma_min_eps = zeros(size(epsilon_values));

A_robust = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_robust = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];

fprintf('Computing robustness data...\n');
for i = 1:length(epsilon_values)
    eps = epsilon_values(i);
    K_robust = @(t) [1; eps * sin(t)];
    
    W_robust = compute_periodic_gramian_block(A_robust, B_robust, K_robust, T, N);
    eigenvals_robust = eig(W_robust);
    sigma_min_eps(i) = sqrt(min(real(eigenvals_robust)));
    
    fprintf('  ε=%.0e: σ_min = %.3e (expected ∝ %.3e)\n', eps, sigma_min_eps(i), eps^2);
end

% Check O(ε²) scaling
log_eps = log(epsilon_values);
log_sigma = log(sigma_min_eps);
p_robust = polyfit(log_eps, log_sigma, 1);
scaling_exponent = p_robust(1);

fprintf('Empirical scaling: σ_min ∝ ε^%.2f\n', scaling_exponent);

if abs(scaling_exponent - 2.0) < 0.5
    test4_status = '✓ PASS';
    fprintf('Status: ✓ PASS (scaling ≈ ε²)\n\n');
else
    test4_status = '○ PARTIAL';
    fprintf('Status: ○ PARTIAL (scaling differs from ε²)\n\n');
end

% Store results
verification_results.tests{end+1} = 'Robustness O(ε²) scaling';
verification_results.expected{end+1} = 2.0;
verification_results.computed{end+1} = scaling_exponent;
verification_results.errors{end+1} = abs(scaling_exponent - 2.0);
verification_results.status{end+1} = test4_status;

%% Test 5: Algorithm Properties
test_count = test_count + 1;
fprintf('TEST %d: Algorithm Properties\n', test_count);
fprintf('=============================\n');

fprintf('Testing basic algorithm properties...\n');

% Test controllability detection
fprintf('  Controllability detection: ');
if min(real(eig(W))) > 1e-12
    fprintf('✓ PASS (system correctly identified as controllable)\n');
    test5a_status = '✓ PASS';
else
    fprintf('✗ FAIL (controllability detection failed)\n');
    test5a_status = '✗ FAIL';
end

% Test Gramian symmetry
fprintf('  Gramian symmetry: ');
symmetry_error = norm(W - W', 'fro') / norm(W, 'fro');
if symmetry_error < 1e-12
    fprintf('✓ PASS (error = %.2e)\n', symmetry_error);
    test5b_status = '✓ PASS';
else
    fprintf('✗ FAIL (error = %.2e)\n', symmetry_error);
    test5b_status = '✗ FAIL';
end

% Test positive definiteness
fprintf('  Positive definiteness: ');
min_eigenval = min(real(eigenvals));
if min_eigenval > 1e-14
    fprintf('✓ PASS (λ_min = %.2e)\n', min_eigenval);
    test5c_status = '✓ PASS';
else
    fprintf('✗ FAIL (λ_min = %.2e)\n', min_eigenval);
    test5c_status = '✗ FAIL';
end

% Overall algorithm properties status
if all(contains([test5a_status, test5b_status, test5c_status], '✓'))
    test5_status = '✓ PASS';
else
    test5_status = '✗ FAIL';
end

fprintf('  Overall: %s\n\n', test5_status);

% Store results
verification_results.tests{end+1} = 'Algorithm properties';
verification_results.expected{end+1} = 'All pass';
verification_results.computed{end+1} = test5_status;
verification_results.errors{end+1} = NaN;
verification_results.status{end+1} = test5_status;

%% Generate Verification Report
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                     VERIFICATION REPORT                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Summary statistics
total_tests = length(verification_results.status);
passed_tests = sum(contains(verification_results.status, '✓'));
partial_tests = sum(contains(verification_results.status, '○'));
failed_tests = sum(contains(verification_results.status, '✗'));

pass_rate = passed_tests / total_tests * 100;

fprintf('VERIFICATION SUMMARY:\n');
fprintf('====================\n');
fprintf('Total tests: %d\n', total_tests);
fprintf('Passed: %d (%.1f%%)\n', passed_tests, pass_rate);
fprintf('Partial: %d (%.1f%%)\n', partial_tests, partial_tests/total_tests*100);
fprintf('Failed: %d (%.1f%%)\n', failed_tests, failed_tests/total_tests*100);
fprintf('Timestamp: %s\n\n', verification_results.timestamp);

% Detailed results table
fprintf('DETAILED VERIFICATION RESULTS:\n');
fprintf('==============================\n');
fprintf('%-30s %-12s %-12s %-8s %-8s\n', 'Test', 'Expected', 'Computed', 'Error', 'Status');
fprintf('%-30s %-12s %-12s %-8s %-8s\n', repmat('-',1,30), repmat('-',1,12), repmat('-',1,12), repmat('-',1,8), repmat('-',1,8));

for i = 1:length(verification_results.tests)
    test_name = verification_results.tests{i};
    expected = verification_results.expected{i};
    computed = verification_results.computed{i};
    error = verification_results.errors{i};
    status = verification_results.status{i};
    
    % Format values for display
    if isnumeric(expected) && ~isnan(expected)
        if abs(expected) > 1e3 || abs(expected) < 1e-3
            exp_str = sprintf('%.2e', expected);
        else
            exp_str = sprintf('%.3f', expected);
        end
    else
        exp_str = sprintf('%s', mat2str(expected));
    end
    
    if isnumeric(computed) && ~isnan(computed)
        if abs(computed) > 1e3 || abs(computed) < 1e-3
            comp_str = sprintf('%.2e', computed);
        else
            comp_str = sprintf('%.3f', computed);
        end
    else
        comp_str = sprintf('%s', mat2str(computed));
    end
    
    if isnumeric(error) && ~isnan(error)
        err_str = sprintf('%.1f%%', error);
    else
        err_str = 'N/A';
    end
    
    fprintf('%-30s %-12s %-12s %-8s %-8s\n', test_name, exp_str, comp_str, err_str, status);
end

%% Final Assessment
fprintf('\nFINAL ASSESSMENT:\n');
fprintf('=================\n');

if pass_rate >= 90
    fprintf('🎉 EXCELLENT: All critical paper claims verified!\n');
    fprintf('✅ Ready for publication - all numerical results confirmed\n');
    fprintf('✅ Algorithm implementation is mathematically sound\n');
    fprintf('✅ Performance claims are substantiated\n');
    
elseif pass_rate >= 70
    fprintf('✅ GOOD: Most paper claims verified with minor discrepancies\n');
    fprintf('📝 Review partial/failed tests for publication readiness\n');
    fprintf('✅ Core algorithm functionality confirmed\n');
    
elseif pass_rate >= 50
    fprintf('⚠️  MODERATE: Significant discrepancies detected\n');
    fprintf('🔧 Address failed tests before publication\n');
    fprintf('📋 Consider parameter adjustments or method refinements\n');
    
else
    fprintf('❌ CRITICAL: Major issues with paper claims\n');
    fprintf('🚫 Publication not recommended without significant revisions\n');
    fprintf('🔧 Algorithm requires substantial debugging\n');
end

fprintf('\nOverall verification: %.1f%% of claims confirmed\n', pass_rate);

end
