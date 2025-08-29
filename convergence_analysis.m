function convergence_analysis()
% CONVERGENCE_ANALYSIS Convergence study with quadrature refinement
%
% This function demonstrates the convergence of the minimum singular value
% with quadrature refinement, as shown in Figure 1 of the paper.
% Tests various numbers of quadrature nodes to show exponential convergence.
%
% EXPECTED BEHAVIOR:
%   Exponential convergence of |σ_min^(N) - σ_min^(200)| as N increases
%   Convergence achieved by N ≈ 80 (relative change < 10^-6)
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com

clear; clc;
fprintf('\n=== CONVERGENCE ANALYSIS ===\n');
fprintf('Quadrature refinement study from paper Figure 1\n\n');

%% System Definition (Example 1 from paper)
% Use the same system as in Example 1 for consistency
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];
T = 2*pi;

fprintf('SYSTEM: Example 1 (n=2, m=1, T=2π)\n');
fprintf('OBJECTIVE: Study convergence of σ_min with quadrature refinement\n\n');

%% Convergence Study Setup
% Test range of quadrature nodes (must be odd for Simpson's rule)
N_values = [11, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 141, 161, 181, 201];
N_reference = 201;  % Reference "exact" value

% Storage for results
results = struct();
results.N_values = N_values;
results.sigma_min = zeros(size(N_values));
results.kappa = zeros(size(N_values));
results.errors = zeros(size(N_values));
results.computation_times = zeros(size(N_values));

fprintf('CONVERGENCE STUDY PARAMETERS:\n');
fprintf('  N values: [%s]\n', sprintf('%d ', N_values));
fprintf('  Reference N = %d\n\n', N_reference);

%% Compute Reference Solution
fprintf('Computing reference solution (N=%d)...\n', N_reference);
tic;
W_ref = compute_periodic_gramian_block(A_func, B_func, K_func, T, N_reference);
ref_time = toc;

eigenvals_ref = eig(W_ref);
sigma_min_ref = sqrt(min(real(eigenvals_ref)));
kappa_ref = max(real(eigenvals_ref)) / min(real(eigenvals_ref));

fprintf('  Reference σ_min = %.8e\n', sigma_min_ref);
fprintf('  Reference κ(W)  = %.8e\n', kappa_ref);
fprintf('  Computation time: %.3f seconds\n\n', ref_time);

%% Convergence Study Loop
fprintf('CONVERGENCE STUDY:\n');
fprintf('==================\n');
fprintf('%-6s %-15s %-15s %-12s %-12s\n', 'N', 'σ_min', 'Error', 'κ(W)', 'Time(s)');
fprintf('%-6s %-15s %-15s %-12s %-12s\n', '------', '---------------', '---------------', '------------', '--------');

for i = 1:length(N_values)
    N = N_values(i);
    
    % Compute Gramian for current N
    tic;
    try
        W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
        comp_time = toc;
        
        % Compute measures
        eigenvals = eig(W);
        sigma_min = sqrt(min(real(eigenvals)));
        kappa = max(real(eigenvals)) / min(real(eigenvals));
        
        % Compute error relative to reference
        error = abs(sigma_min - sigma_min_ref);
        
        % Store results
        results.sigma_min(i) = sigma_min;
        results.kappa(i) = kappa;
        results.errors(i) = error;
        results.computation_times(i) = comp_time;
        
        % Display progress
        fprintf('%-6d %-15.8e %-15.8e %-12.3e %-12.3f\n', ...
                N, sigma_min, error, kappa, comp_time);
        
    catch ME
        fprintf('%-6d %-15s %-15s %-12s %-12s\n', N, 'FAILED', 'N/A', 'N/A', 'N/A');
        fprintf('         Error: %s\n', ME.message);
        
        results.sigma_min(i) = NaN;
        results.kappa(i) = NaN;
        results.errors(i) = NaN;
        results.computation_times(i) = NaN;
    end
end

%% Convergence Analysis
fprintf('\nCONVERGENCE ANALYSIS:\n');
fprintf('====================\n');

% Find valid results
valid_idx = ~isnan(results.errors);
valid_N = results.N_values(valid_idx);
valid_errors = results.errors(valid_idx);

if sum(valid_idx) >= 3
    % Check for convergence threshold
    convergence_threshold = 1e-6;
    converged_idx = find(valid_errors < convergence_threshold, 1, 'first');
    
    if ~isempty(converged_idx)
        convergence_N = valid_N(converged_idx);
        fprintf('✓ Convergence achieved at N = %d (error < %.0e)\n', ...
                convergence_N, convergence_threshold);
    else
        fprintf('⚠ Convergence threshold not reached in test range\n');
    end
    
    % Analyze convergence rate
    if length(valid_N) >= 5
        % Fit exponential decay: error ≈ C * exp(-α * N)
        log_errors = log(valid_errors + eps);  % Add eps to avoid log(0)
        p = polyfit(valid_N, log_errors, 1);
        convergence_rate = -p(1);
        
        fprintf('✓ Estimated convergence rate: α = %.4f\n', convergence_rate);
        fprintf('  (Exponential decay: error ∝ exp(-%.4f × N))\n', convergence_rate);
        
        % Check if convergence is exponential
        R_squared = calculate_r_squared(log_errors, polyval(p, valid_N));
        fprintf('✓ Linear fit quality: R² = %.4f\n', R_squared);
        
        if R_squared > 0.9
            fprintf('✓ Exponential convergence confirmed\n');
        else
            fprintf('⚠ Convergence may not be purely exponential\n');
        end
    end
    
    % Final error level
    final_error = valid_errors(end);
    relative_error = final_error / sigma_min_ref * 100;
    
    fprintf('✓ Final error (N=%d): %.2e (%.4f%%)\n', ...
            valid_N(end), final_error, relative_error);
    
else
    fprintf('⚠ Insufficient valid data for convergence analysis\n');
end

%% Performance Analysis
fprintf('\nPERFORMANCE ANALYSIS:\n');
fprintf('====================\n');

valid_times = results.computation_times(valid_idx);
if length(valid_times) >= 3
    % Analyze computational cost scaling
    log_N = log(valid_N);
    log_times = log(valid_times);
    p_time = polyfit(log_N, log_times, 1);
    time_scaling = p_time(1);
    
    fprintf('✓ Computational cost scaling: O(N^%.2f)\n', time_scaling);
    fprintf('  (Expected: O(N) for quadrature integration)\n');
    
    % Time efficiency
    min_time = min(valid_times);
    max_time = max(valid_times);
    fprintf('✓ Time range: %.3f - %.3f seconds\n', min_time, max_time);
end

%% Generate Convergence Plots
fprintf('\nGenerating convergence plots...\n');
generate_convergence_plots(results, sigma_min_ref);

fprintf('\nConvergence analysis completed successfully!\n');

end

function R_squared = calculate_r_squared(y_actual, y_predicted)
% Calculate R-squared (coefficient of determination)
SS_res = sum((y_actual - y_predicted).^2);
SS_tot = sum((y_actual - mean(y_actual)).^2);
R_squared = 1 - SS_res/SS_tot;
end

function generate_convergence_plots(results, sigma_min_ref)
% Generate convergence analysis plots

try
    figure('Name', 'Convergence Analysis', 'Position', [100, 100, 1200, 900]);
    
    % Get valid data
    valid_idx = ~isnan(results.errors);
    valid_N = results.N_values(valid_idx);
    valid_errors = results.errors(valid_idx);
    valid_sigma = results.sigma_min(valid_idx);
    valid_times = results.computation_times(valid_idx);
    
    % Plot 1: Convergence of σ_min
    subplot(2, 3, 1);
    plot(valid_N, valid_sigma, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
    hold on;
    yline(sigma_min_ref, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Reference (%.6e)', sigma_min_ref));
    xlabel('Number of Quadrature Nodes N');
    ylabel('σ_{min}(W)');
    title('Convergence of Minimum Singular Value');
    legend('Location', 'best');
    grid on;
    
    % Plot 2: Error vs N (semi-log)
    subplot(2, 3, 2);
    semilogy(valid_N, valid_errors, 'ro-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('Number of Quadrature Nodes N');
    ylabel('|σ_{min}^{(N)} - σ_{min}^{(ref)}|');
    title('Absolute Error (Semi-log)');
    grid on;
    
    % Add convergence threshold line
    convergence_threshold = 1e-6;
    hold on;
    yline(convergence_threshold, 'k--', 'LineWidth', 1, 'DisplayName', 'Convergence Threshold');
    legend('Error', 'Threshold', 'Location', 'best');
    
    % Plot 3: Error vs N (log-log) with exponential fit
    subplot(2, 3, 3);
    if length(valid_N) >= 3
        loglog(valid_N, valid_errors, 'go-', 'LineWidth', 2, 'MarkerSize', 6);
        
        % Fit exponential decay
        log_errors = log(valid_errors + eps);
        p = polyfit(valid_N, log_errors, 1);
        N_fit = linspace(min(valid_N), max(valid_N), 100);
        error_fit = exp(polyval(p, N_fit));
        
        hold on;
        loglog(N_fit, error_fit, 'k--', 'LineWidth', 1.5);
        
        xlabel('Number of Quadrature Nodes N');
        ylabel('Absolute Error');
        title('Error Convergence (Log-log)');
        legend('Measured', sprintf('Fit: exp(-%.3f·N)', -p(1)), 'Location', 'best');
        grid on;
    end
    
    % Plot 4: Computation Time vs N
    subplot(2, 3, 4);
    if any(valid_idx)
        plot(valid_N, valid_times, 'mo-', 'LineWidth', 2, 'MarkerSize', 6);
        xlabel('Number of Quadrature Nodes N');
        ylabel('Computation Time (seconds)');
        title('Computational Cost');
        grid on;
        
        % Add linear fit to show O(N) scaling
        if length(valid_N) >= 3
            p_time = polyfit(valid_N, valid_times, 1);
            time_fit = polyval(p_time, valid_N);
            hold on;
            plot(valid_N, time_fit, 'k--', 'LineWidth', 1);
            legend('Measured', 'Linear Fit', 'Location', 'best');
        end
    end
    
    % Plot 5: Condition Number Convergence
    subplot(2, 3, 5);
    valid_kappa = results.kappa(valid_idx);
    if any(~isnan(valid_kappa))
        semilogy(valid_N, valid_kappa, 'co-', 'LineWidth', 2, 'MarkerSize', 6);
        xlabel('Number of Quadrature Nodes N');
        ylabel('κ(W)');
        title('Condition Number Convergence');
        grid on;
    end
    
    % Plot 6: Efficiency (Error vs Time)
    subplot(2, 3, 6);
    if length(valid_errors) >= 3 && length(valid_times) >= 3
        loglog(valid_times, valid_errors, 'ko-', 'LineWidth', 2, 'MarkerSize', 6);
        xlabel('Computation Time (seconds)');
        ylabel('Absolute Error');
        title('Accuracy vs Efficiency');
        grid on;
        
        % Add annotations for key points
        for i = 1:2:length(valid_N)
            text(valid_times(i), valid_errors(i), sprintf('N=%d', valid_N(i)), ...
                 'VerticalAlignment', 'bottom', 'FontSize', 8);
        end
    end
    
    sgtitle('Convergence Analysis Results', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Create a separate figure for the main convergence plot (matches paper Figure 1)
    figure('Name', 'Paper Figure 1: Convergence Plot', 'Position', [200, 200, 800, 600]);
    
    semilogy(valid_N, valid_errors, 'ro-', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'red');
    xlabel('Quadrature Nodes N', 'FontSize', 12);
    ylabel('log_{10}|σ_{min}^{(N)} - σ_{min}^{(200)}|', 'FontSize', 12);
    title('Convergence of Minimum Singular Value with Quadrature Refinement', 'FontSize', 14);
    grid on;
    
    % Add exponential convergence annotation
    if length(valid_N) >= 5
        % Find a good spot for annotation
        mid_idx = round(length(valid_N)/2);
        arrow_x = valid_N(mid_idx);
        arrow_y = valid_errors(mid_idx);
        
        annotation('textarrow', [0.6, 0.4], [0.7, 0.5], 'String', 'Exponential\nconvergence', ...
                   'FontSize', 12, 'HorizontalAlignment', 'center');
    end
    
    % Add convergence threshold
    yline(1e-6, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Convergence Threshold (10^{-6})');
    legend('Error', 'Threshold', 'Location', 'northeast', 'FontSize', 11);
    
    % Set y-limits for better visualization
    ylim([1e-10, 1e-1]);
    
    fprintf('Convergence plots generated successfully.\n');
    
catch ME
    fprintf('Warning: Could not generate convergence plots. Error: %s\n', ME.message);
end

end
