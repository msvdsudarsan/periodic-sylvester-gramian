%% Convergence Analysis with Quadrature Refinement
% This reproduces the convergence study from Section 6.2 and Figure 1 of the paper
% Studies how the minimum singular value converges with increasing quadrature nodes
%
% Author: M. S. V. D. Sudarsan
% Paper: "Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"

clear; clc; close all;

fprintf('=== CONVERGENCE ANALYSIS ===\n');
fprintf('Studying convergence of σ_min with quadrature refinement\n\n');

% System parameters (use Example 1 system for consistency)
n = 2;
m = 1; 
T = 2*pi;

% Define the same system as in Example 1
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

fprintf('Using Example 1 system (n=%d, m=%d, T=%.2f)\n\n', n, m, T);

% Range of quadrature nodes to test (must be odd for Simpson)
N_values = [11, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 141, 161, 181, 201];
N_ref = 201;  % Reference value for "exact" solution

fprintf('Testing N = [%s]\n', num2str(N_values));
fprintf('Using N = %d as reference solution\n\n', N_ref);

% Compute reference solution
fprintf('Computing reference solution with N = %d...\n', N_ref);
W_ref = compute_periodic_gramian_block(A_func, B_func, K_func, T, N_ref);
sigma_min_ref = min(real(eig(W_ref)));
fprintf('Reference σ_min = %.6e\n\n', sigma_min_ref);

% Storage for results
results = zeros(length(N_values), 4);  % [N, sigma_min, abs_error, comp_time]

fprintf('%-6s %-12s %-12s %-12s %-8s\n', 'N', 'σ_min', 'Abs Error', 'Rel Error', 'Time (s)');
fprintf('%s\n', repmat('-', 1, 65));

% Convergence study
for i = 1:length(N_values)
    N = N_values(i);
    
    % Compute Gramian
    tic;
    W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
    comp_time = toc;
    
    % Extract minimum singular value
    sigma_min = min(real(eig(W)));
    
    % Compute errors
    abs_error = abs(sigma_min - sigma_min_ref);
    rel_error = abs_error / abs(sigma_min_ref);
    
    % Store results
    results(i, :) = [N, sigma_min, abs_error, comp_time];
    
    fprintf('%-6d %-12.6e %-12.2e %-12.2e %-8.3f\n', ...
            N, sigma_min, abs_error, rel_error, comp_time);
end

fprintf('\n=== CONVERGENCE ANALYSIS ===\n');

% Find convergence point (relative change < 1e-6)
conv_tol = 1e-6;
conv_idx = find(results(:,3)./abs(sigma_min_ref) < conv_tol, 1, 'first');

if ~isempty(conv_idx)
    fprintf('Convergence achieved at N = %d (relative error < %.0e)\n', ...
            N_values(conv_idx), conv_tol);
else
    fprintf('Convergence tolerance %.0e not achieved in tested range\n', conv_tol);
end

% Estimate convergence rate
if length(N_values) >= 3
    % Fit exponential decay to errors: error ≈ C * exp(-α*N)
    N_fit = N_values(3:end);  % Skip first few points for stability
    errors_fit = results(3:end, 3);
    
    % Log-linear fit
    valid_idx = errors_fit > 0;  % Only positive errors
    if sum(valid_idx) >= 3
        N_fit = N_fit(valid_idx);
        log_errors = log(errors_fit(valid_idx));
        
        % Fit log(error) = log(C) - α*N
        p = polyfit(N_fit, log_errors, 1);
        alpha = -p(1);  % Convergence rate
        C = exp(p(2));  % Pre-factor
        
        fprintf('Estimated convergence rate: α ≈ %.4f\n', alpha);
        fprintf('Error model: |σ_min^(N) - σ_min^(∞)| ≈ %.2e * exp(-%.4f*N)\n', C, alpha);
    end
end

% Create convergence plot (matching Figure 1 from paper)
figure('Position', [100, 100, 800, 600]);

% Main convergence plot
subplot(2, 1, 1);
semilogy(N_values, results(:, 3), 'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
hold on;

% Add exponential fit if available
if exist('alpha', 'var') && exist('C', 'var')
    N_smooth = linspace(min(N_values), max(N_values), 100);
    error_fit = C * exp(-alpha * N_smooth);
    semilogy(N_smooth, error_fit, 'b--', 'LineWidth', 1.5);
    legend('Computed errors', sprintf('Fit: %.2e exp(-%.4f N)', C, alpha), ...
           'Location', 'northeast');
else
    legend('Computed errors', 'Location', 'northeast');
end

grid on;
xlabel('Quadrature Nodes N');
ylabel('|σ_{min}^{(N)} - σ_{min}^{(200)}|');
title('Convergence of Minimum Singular Value');

% Add convergence threshold line
if ~isempty(conv_idx)
    yline(conv_tol * abs(sigma_min_ref), 'g--', 'LineWidth', 1, ...
          'Label', sprintf('Tolerance (%.0e)', conv_tol));
    xline(N_values(conv_idx), 'g--', 'LineWidth', 1, ...
          'Label', sprintf('N = %d', N_values(conv_idx)));
end

% Computation time plot
subplot(2, 1, 2);
plot(N_values, results(:, 4), 'bo-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
grid on;
xlabel('Quadrature Nodes N');
ylabel('Computation Time (s)');
title('Computational Cost vs Accuracy Trade-off');

% Add theoretical O(N) line
N_theory = linspace(min(N_values), max(N_values), 100);
time_theory = results(1, 4) * N_theory / N_values(1);  % Scale from first point
hold on;
plot(N_theory, time_theory, 'r--', 'LineWidth', 1.5);
legend('Computed times', 'O(N) scaling', 'Location', 'northwest');

sgtitle('Convergence Analysis: Quadrature Refinement Study');

% Save convergence data
fprintf('\n=== SAVING RESULTS ===\n');
convergence_data.N_values = N_values;
convergence_data.sigma_min_values = results(:, 2);
convergence_data.absolute_errors = results(:, 3);
convergence_data.computation_times = results(:, 4);
convergence_data.sigma_min_reference = sigma_min_ref;
convergence_data.N_reference = N_ref;

if exist('alpha', 'var')
    convergence_data.convergence_rate = alpha;
    convergence_data.error_prefactor = C;
end

save('convergence_analysis_results.mat', 'convergence_data');
fprintf('Results saved to convergence_analysis_results.mat\n');

% Summary
fprintf('\n=== SUMMARY ===\n');
fprintf('Initial error (N=%d): %.2e\n', N_values(1), results(1,3));
fprintf('Final error (N=%d): %.2e\n', N_values(end), results(end,3));
fprintf('Error reduction: %.1fx\n', results(1,3)/results(end,3));

if ~isempty(conv_idx)
    fprintf('Practical convergence: N ≥ %d\n', N_values(conv_idx));
    fprintf('Recommended N for this system: %d\n', N_values(conv_idx));
else
    fprintf('Recommended N for this system: ≥ %d\n', N_values(end));
end

fprintf('\nNote: This analysis demonstrates exponential convergence\n');
fprintf('of the composite Simpson quadrature rule as expected for\n');
fprintf('smooth periodic integrands.\n');

fprintf('\n=== CONVERGENCE ANALYSIS COMPLETE ===\n');
