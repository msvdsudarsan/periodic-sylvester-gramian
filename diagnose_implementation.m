function diagnose_implementation()
% DIAGNOSE_IMPLEMENTATION 
% This script helps identify why your current implementation 
% doesn't match the TEX paper values

clear; clc;
fprintf('\n=== IMPLEMENTATION DIAGNOSTIC ===\n');
fprintf('Identifying discrepancies with TEX paper values\n\n');

%% Step 1: Test System Parameters
fprintf('STEP 1: System Parameter Verification\n');
fprintf('====================================\n');

% Define system EXACTLY as in TEX paper
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

% Check at specific time points
test_times = [0, pi/2, pi, 3*pi/2];

fprintf('System matrices at key time points:\n');
for i = 1:length(test_times)
    t = test_times(i);
    A_t = A_func(t);
    B_t = B_func(t);
    K_t = K_func(t);
    
    fprintf('t = %.2f:\n', t);
    fprintf('  A(t) = [%6.3f, %6.3f; %6.3f, %6.3f]\n', A_t(1,1), A_t(1,2), A_t(2,1), A_t(2,2));
    fprintf('  B(t) = [%6.3f, %6.3f; %6.3f, %6.3f]\n', B_t(1,1), B_t(1,2), B_t(2,1), B_t(2,2));
    fprintf('  K(t) = [%6.3f; %6.3f]\n', K_t(1), K_t(2));
    fprintf('\n');
end

%% Step 2: Test Your Current Implementation
fprintf('STEP 2: Testing Current Implementation\n');
fprintf('====================================\n');

T = 2*pi;
N = 101;

fprintf('Computing with your current implementation...\n');
try
    % Try to call your current function
    W_current = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
    
    if size(W_current, 1) ~= 4 || size(W_current, 2) ~= 4
        fprintf('ERROR: Wrong Gramian dimensions!\n');
        fprintf('Expected: 4x4, Got: %dx%d\n', size(W_current,1), size(W_current,2));
    else
        eigenvals_current = eig(W_current);
        eigenvals_current = sort(real(eigenvals_current), 'descend');
        
        sigma_min_current = sqrt(min(eigenvals_current));
        kappa_current = max(eigenvals_current) / min(eigenvals_current);
        
        fprintf('Your results:\n');
        fprintf('  Eigenvalues: [%.6e, %.6e, %.6e, %.6e]\n', ...
                eigenvals_current(1), eigenvals_current(2), eigenvals_current(3), eigenvals_current(4));
        fprintf('  σ_min = %.6e\n', sigma_min_current);
        fprintf('  κ(W)  = %.6e\n', kappa_current);
    end
    
catch ME
    fprintf('ERROR in current implementation: %s\n', ME.message);
    W_current = [];
end

%% Step 3: Compare with Expected Values
fprintf('\nSTEP 3: Comparison with TEX Paper\n');
fprintf('=================================\n');

sigma_min_paper = 1.25e-2;
kappa_paper = 8.4e3;

fprintf('TEX paper claims:\n');
fprintf('  σ_min = %.6e\n', sigma_min_paper);
fprintf('  κ(W)  = %.6e\n', kappa_paper);

if ~isempty(W_current)
    sigma_error = abs(sigma_min_current - sigma_min_paper) / sigma_min_paper * 100;
    kappa_error = abs(kappa_current - kappa_paper) / kappa_paper * 100;
    
    fprintf('Errors:\n');
    fprintf('  σ_min error: %.1f%%\n', sigma_error);
    fprintf('  κ(W) error:  %.1f%%\n', kappa_error);
    
    if sigma_error > 100 || kappa_error > 100
        fprintf('DIAGNOSIS: MAJOR IMPLEMENTATION ERRORS\n');
        fprintf('Possible causes:\n');
        fprintf('1. Wrong system parameters in your code\n');
        fprintf('2. Incorrect algorithm implementation\n');
        fprintf('3. Dimensional errors in matrix operations\n');
        fprintf('4. Wrong quadrature or ODE solver settings\n');
    elseif sigma_error > 20 || kappa_error > 20
        fprintf('DIAGNOSIS: SIGNIFICANT NUMERICAL ISSUES\n');
        fprintf('Possible causes:\n');
        fprintf('1. Insufficient quadrature accuracy\n');
        fprintf('2. ODE solver tolerance too loose\n');
        fprintf('3. Parameter scaling issues\n');
    else
        fprintf('DIAGNOSIS: ACCEPTABLE AGREEMENT\n');
        fprintf('Implementation appears correct\n');
    end
end

%% Step 4: Test Algorithm Components
fprintf('\nSTEP 4: Component Testing\n');
fprintf('========================\n');

% Test Simpson weights
w = simpson_weights(N, T);
fprintf('Simpson weights: first=%f, last=%f, sum=%.6f\n', w(1), w(end), sum(w));
expected_sum = T;
weight_error = abs(sum(w) - expected_sum) / expected_sum * 100;
if weight_error > 1
    fprintf('ERROR: Simpson weights incorrect (%.1f%% error)\n', weight_error);
else
    fprintf('Simpson weights: OK\n');
end

% Test periodicity
periodicity_error = max([
    norm(A_func(0) - A_func(T), 'fro'), 
    norm(B_func(0) - B_func(T), 'fro'),
    norm(K_func(0) - K_func(T), 'fro')
]);
fprintf('Periodicity error: %.2e\n', periodicity_error);
if periodicity_error > 1e-12
    fprintf('WARNING: System may not be exactly periodic\n');
else
    fprintf('Periodicity: OK\n');
end

%% Step 5: Recommendations
fprintf('\nSTEP 5: Recommendations\n');
fprintf('======================\n');

if ~isempty(W_current)
    if sigma_error > 50 || kappa_error > 50
        fprintf('CRITICAL: Replace your implementation with corrected version\n');
        fprintf('Action: Use the corrected codes I provided\n');
    elseif sigma_error > 15 || kappa_error > 15
        fprintf('MODERATE: Debug your current implementation\n');
        fprintf('Action: Check parameter values and algorithm details\n');
    else
        fprintf('MINOR: Fine-tune numerical parameters\n');
        fprintf('Action: Adjust ODE tolerances or quadrature nodes\n');
    end
else
    fprintf('CRITICAL: Fix basic implementation errors first\n');
    fprintf('Action: Use the corrected compute_periodic_gramian_block.m\n');
end

fprintf('\nDiagnostic completed.\n');

end

function w = simpson_weights(N, T)
% Simpson's rule weights for diagnostic
if mod(N, 2) == 0
    error('N must be odd for Simpson rule');
end

h = T / (N - 1);
w = zeros(1, N);
w(1) = h/3; w(N) = h/3;

for i = 2:N-1
    if mod(i-1, 2) == 0
        w(i) = 4*h/3;
    else
        w(i) = 2*h/3;
    end
end

end
