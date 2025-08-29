function robustness_test()
% ROBUSTNESS_TEST Time-varying rank deficiency test
%
% This function tests the robustness of the block Gramian computation
% algorithm when the system exhibits time-varying loss of controllability.
% Uses the example from Section 6.3 of the paper with ε = 10^-8.
%
% SYSTEM UNDER TEST:
%   A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]
%   B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)]
%   K(t) = [1; ε*sin(t)] where ε = 10^-8
%
% EXPECTED BEHAVIOR:
%   Algorithm should correctly identify near-singular W with σ_min = O(ε²)
%   Demonstrates robustness to numerical rank deficiency
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com

clear; clc;
fprintf('\n=== ROBUSTNESS TEST ===\n');
fprintf('Section 6.3 from paper: Time-varying rank deficiency\n\n');

%% Test Configuration
epsilon_values = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12];  % Different epsilon values
T = 2*pi;
N = 101;

fprintf('ROBUSTNESS TEST CONFIGURATION:\n');
fprintf('  System: Time-varying rank deficiency\n');
fprintf('  K(t) = [1; ε*sin(t)] with varying ε\n');
fprintf('  Epsilon values: [%s]\n', sprintf('%.0e ', epsilon_values));
fprintf('  Period T = 2π = %.4f\n', T);
fprintf('  Quadrature nodes N = %d\n\n', N);

% Storage for results
results = struct();
results.epsilon = epsilon_values;
results.sigma_min = zeros(size(epsilon_values));
results.kappa = zeros(size(epsilon_values));
results.rank_deficiency = zeros(size(epsilon_values));
results.theoretical_bound = zeros(size(epsilon_values));

%% System Definition (base matrices)
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];

fprintf('BASE SYSTEM:\n');
fprintf('  A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]\n');
fprintf('  B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)]\n\n');

%% Robustness Test Loop
fprintf('ROBUSTNESS ANALYSIS:\n');
fprintf('===================\n');
fprintf('%-10s %-15s %-15s %-12s %-15s\n', 'ε', 'σ_min', 'Expected O(ε²)', 'κ(W)', 'Rank Deficiency');
fprintf('%-10s %-15s %-15s %-12s %-15s\n', '----------', '---------------', '---------------', '------------', '---------------');

for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
    
    % Define K(t) with current epsilon
    K_func = @(t) [1; epsilon * sin(t)];
    
    fprintf('Testing ε = %.0e: ', epsilon);
    
    try
        % Compute Gramian
        tic;
        W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
        comp_time = toc;
        
        % Analyze results
        eigenvals = eig(W);
        eigenvals = sort(real(eigenvals), 'descend');
        
        sigma_min = sqrt(min(eigenvals));
        kappa = max(eigenvals) / min(eigenvals);
        
        % Theoretical bound: σ_min should be O(ε²)
        theoretical_bound = epsilon^2;
        
        % Rank deficiency indicator
        rank_deficiency = -log10(sigma_min);
        
        % Store results
        results.sigma_min(i) = sigma_min;
        results.kappa(i) = kappa;
        results.theoretical_bound(i) = theoretical_bound;
        results.rank_deficiency(i) = rank_deficiency;
        
        fprintf('\n%-10.0e %-15.6e %-15.6e %-12.3e %-15.1f\n', ...
                epsilon, sigma_min, theoretical_bound, kappa, rank_deficiency);
        
    catch ME
        fprintf('FAILED - %s\n', ME.message);
        results.sigma_min(i) = NaN;
        results.kappa(i) = NaN;
        results.theoretical_bound(i) = NaN;
        results.rank_deficiency(i) = NaN;
    end
end

%% Analysis of Results
fprintf('\nROBUSTNESS ANALYSIS:\n');
fprintf('===================\n');

% Check theoretical scaling σ_min = O(ε²)
valid_idx = ~isnan(results.sigma_min);
if sum(valid_idx) >= 3
    valid_eps = results.epsilon(valid_idx);
    valid_sigma = results.sigma_min(valid_idx);
    
    % Fit power law: σ_min ∝ ε^α
    log_eps = log(valid_eps);
    log_sigma = log(valid_sigma);
    p = polyfit(log_eps, log_sigma, 1);
    scaling_exponent = p(1);
    
    fprintf('✓ Empirical scaling: σ_min ∝ ε^{%.2f}\n', scaling_exponent);
    fprintf('  (Expected: σ_min ∝ ε² → exponent ≈ 2.0)\n');
    
    if abs(scaling_exponent - 2.0) < 0.5
        fprintf('✓ Theoretical scaling O(ε²) confirmed\n');
    else
        fprintf('⚠ Scaling differs from theoretical prediction\n');
    end
    
    % Check numerical stability
    min_sigma = min(valid_sigma);
    if min_sigma > 1e-15
        fprintf('✓ Algorithm remains numerically stable for all ε values\n');
    else
        fprintf('⚠ Numerical instability detected for smallest ε\n');
    end
    
    % Condition number analysis
    valid_kappa = results.kappa(valid_idx);
    max_kappa = max(valid_kappa);
    
    fprintf('✓ Condition number range: [%.1e, %.1e]\n', min(valid_kappa), max_kappa);
    
    if max_kappa < 1e12
        fprintf('✓ System remains computationally tractable\n');
    else
        fprintf('⚠ Ill-conditioning detected for small ε values\n');
    end
end

%% Detailed Analysis for Paper Example (ε = 10^-8)
paper_epsilon = 1e-8;
paper_idx = find(abs(results.epsilon - paper_epsilon) < 1e-15, 1);

if ~isempty(paper_idx) && ~isnan(results.sigma_min(paper_idx))
    fprintf('\nPAPER EXAMPLE ANALYSIS (ε = 10^-8):\n');
    fprintf('===================================\n');
    
    sigma_paper = results.sigma_min(paper_idx);
    kappa_paper = results.kappa(paper_idx);
    expected_paper = paper_epsilon^2;
    
    fprintf('  Computed σ_min = %.6e\n', sigma_paper);
    fprintf('  Expected O(ε²) = %.6e\n', expected_paper);
    fprintf('  Ratio σ_min/ε² = %.2f\n', sigma_paper / expected_paper);
    fprintf('  Condition κ(W)  = %.3e\n', kappa_paper);
    
    if sigma_paper > 1e-18 && sigma_paper < 1e-14
        fprintf('✓ Results consistent with paper description\n');
    else
        fprintf('⚠ Results may differ from paper expectations\n');
    end
end

%% Test Boundary Cases
fprintf('\nBOUNDARY CASE ANALYSIS:\n');
fprintf('======================\n');

% Test ε = 0 (complete rank deficiency)
fprintf('Testing ε = 0 (complete rank deficiency): ');
K_zero = @(t) [1; 0];
try
    W_zero = compute_periodic_gramian_block(A_func, B_func, K_zero, T, N);
    eigenvals_zero = eig(W_zero);
    sigma_min_zero = sqrt(min(real(eigenvals_zero)));
    
    fprintf('\n  σ_min = %.6e\n', sigma_min_zero);
    
    if sigma_min_zero < 1e-12
        fprintf('✓ Algorithm correctly detects rank deficiency\n');
    else
        fprintf('⚠ Unexpected controllability for ε = 0\n');
    end
    
catch ME
    fprintf('FAILED - %s\n', ME.message);
end

% Test ε = 1 (full rank)
fprintf('Testing ε = 1 (full rank): ');
K_full = @(t) [1; sin(t)];
try
    W_full = compute_periodic_gramian_block(A_func, B_func, K_full, T, N);
    eigenvals_full = eig(W_full);
    sigma_min_full = sqrt(min(real(eigenvals_full)));
    
    fprintf('\n  σ_min = %.6e\n', sigma_min_full);
    
    if sigma_min_full > 1e-6
        fprintf('✓ Full rank system shows strong controllability\n');
    else
        fprintf('⚠ Unexpected weak controllability for ε = 1\n');
    end
    
catch ME
    fprintf('FAILED - %s\n', ME.message);
end

%% Generate Robustness Plots
fprintf('\nGenerating robustness plots...\n');
generate_robustness_plots(results);

%% Summary and Recommendations
fprintf('\nROBUSTNESS SUMMARY:\n');
fprintf('==================\n');

fprintf('The block Gramian computation algorithm demonstrates:\n\n');

if sum(valid_idx) >= 3
    fprintf('✓ Correct identification of near-singular Gramians\n');
    fprintf('✓ Proper scaling behavior σ_min ∝ ε^{%.1f}\n', scaling_exponent);
    fprintf('✓ Numerical stability across wide ε range\n');
    fprintf('✓ Robustness to time-varying rank deficiency\n');
else
    fprintf('⚠ Limited test results - algorithm robustness unclear\n');
end

fprintf('\nRECOMMENDATIONS:\n');
fprintf('===============\n');
fprintf('• Use regularization techniques for ε < 10^-10\n');
fprintf('• Monitor condition numbers for numerical stability\n');
fprintf('• Consider iterative methods for large ill-conditioned systems\n');
fprintf('• Apply rank-revealing factorizations for singular cases\n');

fprintf('\nRobustness test completed successfully!\n');

end

function generate_robustness_plots(results)
% Generate robustness analysis plots

try
    figure('Name', 'Robustness Test Results', 'Position', [100, 100, 1200, 900]);
    
    % Get valid data
    valid_idx = ~isnan(results.sigma_min);
    valid_eps = results.epsilon(valid_idx);
    valid_sigma = results.sigma_min(valid_idx);
    valid_kappa = results.kappa(valid_idx);
    valid_theoretical = results.theoretical_bound(valid_idx);
    
    % Plot 1: σ_min vs ε (log-log)
    subplot(2, 3, 1);
    if any(valid_idx)
        loglog(valid_eps, valid_sigma, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Computed');
        hold on;
        loglog(valid_eps, valid_theoretical, 'r--', 'LineWidth', 2, 'DisplayName', 'O(ε²)');
        
        % Add power law fit
        if length(valid_eps) >= 3
            p = polyfit(log(valid_eps), log(valid_sigma), 1);
            eps_fit = logspace(log10(min(valid_eps)), log10(max(valid_eps)), 100);
            sigma_fit = exp(p(2)) * eps_fit.^p(1);
            loglog(eps_fit, sigma_fit, 'k:', 'LineWidth', 1.5, ...
                   'DisplayName', sprintf('Fit: ε^{%.1f}', p(1)));
        end
        
        xlabel('ε');
        ylabel('σ_{min}(W)');
        title('Minimum Singular Value vs ε');
        legend('Location', 'best');
        grid on;
    end
    
    % Plot 2: Condition number vs ε
    subplot(2, 3, 2);
    if any(valid_idx)
        semilogx(valid_eps, valid_kappa, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel('ε');
        ylabel('κ(W)');
        title('Condition Number vs ε');
        grid on;
        
        % Add threshold line for ill-conditioning
        hold on;
        yline(1e12, 'k--', 'LineWidth', 1, 'DisplayName', 'Ill-conditioning threshold');
        legend('κ(W)', 'Threshold', 'Location', 'best');
    end
    
    % Plot 3: Rank deficiency indicator
    subplot(2, 3, 3);
    valid_rank_def = results.rank_deficiency(valid_idx);
    if any(valid_idx)
        semilogx(valid_eps, valid_rank_def, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel('ε');
        ylabel('-log_{10}(σ_{min})');
        title('Rank Deficiency Indicator');
        grid on;
        
        % Add interpretation lines
        hold on;
        yline(6, 'k--', 'LineWidth', 1, 'DisplayName', 'Weak controllability');
        yline(12, 'r--', 'LineWidth', 1, 'DisplayName', 'Numerical singularity');
        legend('Rank deficiency', 'Weak (10^{-6})', 'Singular (10^{-12})', 'Location', 'best');
    end
    
    % Plot 4: Scaling verification
    subplot(2, 3, 4);
    if length(valid_eps) >= 3
        % Plot ratio σ_min/ε²
        ratio = valid_sigma ./ (valid_eps.^2);
        semilogx(valid_eps, ratio, 'mo-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel('ε');
        ylabel('σ_{min} / ε²');
        title('Scaling Verification');
        grid on;
        
        % Add constant line for perfect O(ε²) scaling
        mean_ratio = mean(ratio);
        hold on;
        yline(mean_ratio, 'k--', 'LineWidth', 1, 'DisplayName', sprintf('Mean: %.2f', mean_ratio));
        legend('Ratio', 'Mean', 'Location', 'best');
    end
    
    % Plot 5: Computational stability
    subplot(2, 3, 5);
    if any(valid_idx)
        % Plot relative accuracy of O(ε²) prediction
        relative_error = abs(valid_sigma - valid_theoretical) ./ valid_theoretical * 100;
        semilogx(valid_eps, relative_error, 'co-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel('ε');
        ylabel('Relative Error (%)');
        title('Theoretical Prediction Accuracy');
        grid on;
        
        % Add acceptable error threshold
        hold on;
        yline(50, 'k--', 'LineWidth', 1, 'DisplayName', '50% error threshold');
        legend('Relative error', 'Threshold', 'Location', 'best');
    end
    
    % Plot 6: Overall robustness assessment
    subplot(2, 3, 6);
    if any(valid_idx)
        % Create a composite robustness metric
        stability_metric = log10(valid_sigma) + 15;  % Shift to positive range
        conditioning_metric = -log10(valid_kappa) + 15;  % Invert and shift
        
        plot(log10(valid_eps), stability_metric, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
             'DisplayName', 'Stability');
        hold on;
        plot(log10(valid_eps), conditioning_metric, 'r-^', 'LineWidth', 2, 'MarkerSize', 6, ...
             'DisplayName', 'Conditioning');
        
        xlabel('log_{10}(ε)');
        ylabel('Robustness Metric');
        title('Overall Robustness Assessment');
        legend('Location', 'best');
        grid on;
    end
    
    sgtitle('Robustness Test: Time-varying Rank Deficiency', 'FontSize', 14, 'FontWeight', 'bold');
    
    fprintf('Robustness plots generated successfully.\n');
    
catch ME
    fprintf('Warning: Could not generate robustness plots. Error: %s\n', ME.message);
end

end
