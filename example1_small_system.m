function example1_small_system()
% EXAMPLE1_SMALL_SYSTEM Validation example with n=2, m=1 system
%
% Author: M. S. V. D. Sudarsan

    fprintf('=== Example 1: Small System Validation ===\n');
    
    % System parameters
    n = 2; m = 1; T = 2*pi;
    
    % Define periodic matrices
    A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
    B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
    K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];
    
    % Compute Gramian
    fprintf('Computing Gramian with N=101 nodes...\n');
    tic;
    W = compute_periodic_gramian(A_func, B_func, K_func, T, 101);
    comp_time = toc;
    
    % Analysis
    eigenvals = eig(W);
    sigma_min = min(real(eigenvals));
    sigma_max = max(real(eigenvals));
    cond_num = sigma_max / sigma_min;
    
    fprintf('\nResults:\n');
    fprintf('  Computation time: %.4f seconds\n', comp_time);
    fprintf('  Minimum eigenvalue: %.6e\n', sigma_min);
    fprintf('  Maximum eigenvalue: %.6e\n', sigma_max);
    fprintf('  Condition number: %.2e\n', cond_num);
    
    if sigma_min > 1e-10
        fprintf('  ✓ System is controllable\n');
    else
        fprintf('  ✗ System may not be controllable\n');
    end
    
    % Test convergence
    fprintf('\n--- Convergence Test ---\n');
    N_values = [21, 41, 61, 81, 101];
    sigma_mins = zeros(length(N_values), 1);
    
    for i = 1:length(N_values)
        W_temp = compute_periodic_gramian(A_func, B_func, K_func, T, N_values(i));
        sigma_mins(i) = min(real(eig(W_temp)));
        fprintf('N = %3d: σ_min = %.6e\n', N_values(i), sigma_mins(i));
    end
    
    % Check convergence
    rel_change = abs(sigma_mins(end) - sigma_mins(end-1)) / abs(sigma_mins(end));
    fprintf('Relative change (last two): %.2e\n', rel_change);
    if rel_change < 1e-6
        fprintf('✓ Converged to desired tolerance\n');
    end
end
