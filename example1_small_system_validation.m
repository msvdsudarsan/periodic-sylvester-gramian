function example1_small_system_validation()
%EXAMPLE1_SMALL_SYSTEM_VALIDATION Small system validation (Example 1 from paper)
%
% Validates the implementation using the small system from Section 6.1:
% n=2, m=1, T=2π with known expected results:
% - σ_min(W) ≈ 1.25×10⁻²
% - κ(W) ≈ 8.4×10³

fprintf('=== EXAMPLE 1: SMALL SYSTEM VALIDATION ===\n');
fprintf('Testing n=2, m=1 system with T=2π\n\n');

% Define system parameters exactly as in paper
T = 2*pi;
n = 2;
m = 1;

% System matrices as defined in Section 6.1
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

% Display system definition
fprintf('System definition:\n');
fprintf('A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]\n');
fprintf('B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)]\n');
fprintf('K(t) = [1 + 0.2*cos(t); 0.5*sin(t)]\n\n');

% Compute Gramian using N=101 nodes (as mentioned in paper)
N = 101;
fprintf('Computing Gramian with N=%d quadrature nodes...\n', N);

tic;
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
computation_time = toc;

% Analyze results
eigenvals = eig(W);
sigma_min = min(real(eigenvals));
sigma_max = max(real(eigenvals));
kappa = sigma_max / sigma_min;

% Display results
fprintf('\n=== RESULTS ===\n');
fprintf('Computation time: %.3f seconds\n', computation_time);
fprintf('Gramian size: %dx%d\n', size(W, 1), size(W, 2));
fprintf('σ_min(W) = %.6e\n', sigma_min);
fprintf('σ_max(W) = %.6e\n', sigma_max);
fprintf('κ(W) = %.6e\n', kappa);

% Check if system is controllable
if sigma_min > 1e-12
    fprintf('System is CONTROLLABLE (σ_min > 0)\n');
else
    fprintf('System may not be controllable (σ_min ≈ 0)\n');
end

% Compare with expected paper results
fprintf('\n=== PAPER COMPARISON ===\n');
paper_sigma_min = 1.25e-2;
paper_kappa = 8.4e3;

rel_error_sigma = abs(sigma_min - paper_sigma_min) / paper_sigma_min;
rel_error_kappa = abs(kappa - paper_kappa) / paper_kappa;

fprintf('Expected σ_min: %.3e, Computed: %.3e, Relative error: %.1f%%\n', ...
        paper_sigma_min, sigma_min, rel_error_sigma*100);
fprintf('Expected κ(W): %.3e, Computed: %.3e, Relative error: %.1f%%\n', ...
        paper_kappa, kappa, rel_error_kappa*100);

% Validation assessment
tolerance = 0.20; % 20% tolerance for numerical differences
if rel_error_sigma < tolerance && rel_error_kappa < tolerance
    fprintf('\n✓ VALIDATION SUCCESSFUL: Results match paper within %.0f%% tolerance\n', tolerance*100);
    validation_status = 'PASSED';
else
    fprintf('\n⚠ VALIDATION WARNING: Results differ from paper by >%.0f%%\n', tolerance*100);
    validation_status = 'NEEDS REVIEW';
    
    % Additional diagnostics
    fprintf('\nDiagnostic information:\n');
    fprintf('- Check if system definition matches paper exactly\n');
    fprintf('- Consider increasing quadrature nodes N\n');
    fprintf('- Verify ODE solver tolerances\n');
end

% Display Gramian properties
fprintf('\n=== GRAMIAN ANALYSIS ===\n');
fprintf('Gramian is symmetric: %s\n', mat2str(issymmetric(W, 1e-12)));
fprintf('Gramian is positive semidefinite: %s\n', mat2str(all(eigenvals >= -1e-12)));
fprintf('Condition number: %.3e\n', kappa);

% Summary
fprintf('\n=== SUMMARY ===\n');
fprintf('Example 1 validation: %s\n', validation_status);
fprintf('System controllability: %s\n', sigma_min > 1e-12);
fprintf('Gramian computation: %.3f seconds\n', computation_time);

end
