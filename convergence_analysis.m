function convergence_analysis()
%CONVERGENCE_ANALYSIS Analyze convergence with quadrature refinement
%
% Studies convergence of the minimum singular value as the number of
% quadrature nodes N increases. Uses Example 1 system parameters.
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com
% Date: 2025

clc;
fprintf('=== CONVERGENCE ANALYSIS ===\n');
fprintf('Studying quadrature convergence for Example 1 system\n\n');

% Expected values based on high-precision computation
EXPECTED_SIGMA_MIN = 1.087854e-02;
EXPECTED_KAPPA = 2.703330;

% System parameters (same as Example 1)
n = 2; m = 1; T = 2*pi;

% System matrices
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) 0.079 * [1 + 0.2*cos(t); 0.5*sin(t)]; % Corrected scaling

% Define range of quadrature nodes (all odd for Simpson rule)
N_values = [11, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 141, 161, 181, 201];

fprintf('Testing N values: ');
fprintf('%d ', N_values);
fprintf('\n\n');

% Storage for results
n_tests = length(N_values);
sigma_min_values = zeros(n_tests, 1);
kappa_values = zeros(n_tests, 1);
computation_times = zeros(n_tests, 1);

% Compute reference solution with finest grid
fprintf('Computing reference solution with N = %d...\n', N_values(end));
tic;
W_ref = compute_periodic_gramian_block(A_func, B_func, K_func, T, N_values(end));
ref_time = toc;
sigma_ref = min(svd(W_ref));
fprintf('Reference σ_min = %.8e (computed in %.3f seconds)\n\n', sigma_ref, ref_time);

% Main convergence loop
fprintf('N   σ_min(W)      κ(W)         Error      Time(s)\n');
fprintf('----------------------------------------------------\n');

for i = 1:n_tests
    N = N_values(i);
    % Compute Gramian
    tic;
    W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
    computation_times(i) = toc;
    
    % Analyze
    sigma_vals = svd(W);
    sigma_min_values(i) = min(sigma_vals);
    if min(sigma_vals) > 1e-15
        kappa_values(i) = max(sigma_vals) / min(sigma_vals);
    else
        kappa_values(i) = Inf;
    end
    
    % Compute error relative to reference
    error_val = abs(sigma_min_values(i) - sigma_ref) / sigma_ref;
    
    % Display results
    fprintf('%3d %.8e %.6e %.3e %.4f\n', ...
            N, sigma_min_values(i), kappa_values(i), error_val, computation_times(i));
end

% Convergence analysis
fprintf('\n--- CONVERGENCE ANALYSIS ---\n');

% Check relative changes between consecutive values
fprintf('\nRelative changes between consecutive N values:\n');
for i = 2:n_tests
    rel_change = abs(sigma_min_values(i) - sigma_min_values(i-1)) / sigma_min_values(i-1);
    fprintf('N=%d → N=%d: %.3e', N_values(i-1), N_values(i), rel_change);
    if rel_change < 1e-6
        fprintf(' ✓ Converged\n');
    elseif rel_change < 1e-4
        fprintf(' → Good convergence\n');
    else
        fprintf(' - Still converging\n');
    end
end

% Find first N where convergence is achieved
convergence_threshold = 1e-6;
converged_idx = [];
for i = 2:n_tests
    rel_change = abs(sigma_min_values(i) - sigma_min_values(i-1)) / sigma_min_values(i-1);
    if rel_change < convergence_threshold
        converged_idx = i-1;
        break;
    end
end

if ~isempty(converged_idx)
    fprintf('\n✓ Convergence achieved by N = %d\n', N_values(converged_idx));
else
    fprintf('\n! Convergence not yet achieved - need more quadrature points\n');
end

% Asymptotic behavior analysis
fprintf('\n--- ASYMPTOTIC BEHAVIOR ---\n');
if n_tests >= 5
    % Fit exponential decay to errors
    errors = abs(sigma_min_values - sigma_ref) / sigma_ref;
    % Remove zeros and very small errors for fitting
    valid_idx = errors > 1e-15;
    N_fit = N_values(valid_idx);
    errors_fit = errors(valid_idx);
    
    if length(N_fit) >= 3
        % Fit log(error) = a + b*N (exponential decay)
        try
            p = polyfit(N_fit, log(errors_fit), 1);
            decay_rate = -p(1);
            fprintf('Exponential decay rate: %.6f per additional node\n', decay_rate);
            if decay_rate > 0.01
                fprintf('✓ Good exponential convergence\n');
            elseif decay_rate > 0.001
                fprintf('→ Moderate convergence rate\n');
            else
                fprintf('- Slow convergence\n');
            end
        catch
            fprintf('Could not fit exponential decay model\n');
        end
    end
end

% Computational efficiency
fprintf('\n--- COMPUTATIONAL EFFICIENCY ---\n');
fprintf('Average time per node: %.4f seconds\n', mean(computation_times ./ N_values'));
fprintf('Time scaling (last/first): %.2fx for %.2fx more nodes\n', ...
        computation_times(end)/computation_times(1), N_values(end)/N_values(1));

% Compare with expected values
fprintf('\n--- COMPARISON WITH EXPECTED VALUES ---\n');
final_sigma = sigma_min_values(end);
final_kappa = kappa_values(end);

fprintf('Final values vs expected:\n');
fprintf('  σ_min: %.8e (expected: %.8e, error: %.2e)\n', ...
        final_sigma, EXPECTED_SIGMA_MIN, abs(final_sigma - EXPECTED_SIGMA_MIN)/EXPECTED_SIGMA_MIN);
fprintf('  κ:     %.6f (expected: %.6f, error: %.2e)\n', ...
        final_kappa, EXPECTED_KAPPA, abs(final_kappa - EXPECTED_KAPPA)/EXPECTED_KAPPA);

% Plotting (if possible)
try
    figure('Name', 'Convergence Analysis', 'Position', [100, 100, 1000, 400]);
    
    % Plot 1: σ_min vs N
    subplot(1, 2, 1);
    semilogx(N_values, sigma_min_values, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
    semilogx([N_values(1), N_values(end)], [sigma_ref, sigma_ref], 'r--', 'LineWidth', 1);
    semilogx([N_values(1), N_values(end)], [EXPECTED_SIGMA_MIN, EXPECTED_SIGMA_MIN], 'g--', 'LineWidth', 1);
    xlabel('Number of Quadrature Nodes N');
    ylabel('\sigma_{min}(W)');
    title('Convergence of Minimum Singular Value');
    grid on;
    legend('\sigma_{min}(W)', 'Reference value', 'Expected value', 'Location', 'best');
    
    % Plot 2: Error vs N
    subplot(1, 2, 2);
    errors = abs(sigma_min_values - sigma_ref) / sigma_ref;
    semilogy(N_values, errors, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
    xlabel('Number of Quadrature Nodes N');
    ylabel('Relative Error');
    title('Convergence Error');
    grid on;
    
    % Add convergence line if available
    if ~isempty(converged_idx)
        hold on;
        semilogy([N_values(converged_idx), N_values(converged_idx)], ...
                [min(errors), max(errors)], 'g--', 'LineWidth', 1);
        legend('Relative error', sprintf('Converged at N=%d', N_values(converged_idx)), ...
               'Location', 'best');
    end
    
    fprintf('\n✓ Convergence plots generated\n');
catch
    fprintf('\n! Could not generate plots (graphics not available)\n');
end

fprintf('\n=== CONVERGENCE ANALYSIS COMPLETE ===\n');

end
