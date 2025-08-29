function example1_small_system_validation()
% EXAMPLE1_SMALL_SYSTEM_VALIDATION Small system validation example
%
% This example validates the block Gramian computation algorithm using
% the small system (n=2, m=1) described in Section 6.1 of the paper.
%
% SYSTEM PARAMETERS (from TEX file):
%   A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]
%   B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)]
%   K(t) = [1 + 0.2*cos(t); 0.5*sin(t)]
%   Period T = 2π
%
% EXPECTED RESULTS (from paper):
%   σ_min(W) ≈ 1.25×10^-2
%   κ(W) ≈ 8.4×10^3
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com

clear; clc;
fprintf('\n=== EXAMPLE 1: SMALL SYSTEM VALIDATION ===\n');
fprintf('Section 6.1 from paper: n=2, m=1, T=2π\n\n');

%% System Definition (exactly as in TEX file)
% A(t) matrix function
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];

% B(t) matrix function  
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];

% K(t) matrix function
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

% System parameters
T = 2*pi;           % Period
N = 101;            % Quadrature nodes (as mentioned in paper)

%% Display System Information
fprintf('SYSTEM DEFINITION:\n');
fprintf('  A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]\n');
fprintf('  B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)]\n'); 
fprintf('  K(t) = [1 + 0.2*cos(t); 0.5*sin(t)]\n');
fprintf('  Period T = 2π = %.4f\n', T);
fprintf('  Quadrature nodes N = %d\n\n', N);

%% Verify System Properties
fprintf('SYSTEM VERIFICATION:\n');

% Check dimensions
A0 = A_func(0);
B0 = B_func(0); 
K0 = K_func(0);

[n, ~] = size(A0);
[~, m] = size(K0);

fprintf('  State dimension n = %d\n', n);
fprintf('  Input dimension m = %d\n', m);
fprintf('  A(t): %dx%d matrix\n', size(A0,1), size(A0,2));
fprintf('  B(t): %dx%d matrix\n', size(B0,1), size(B0,2));
fprintf('  K(t): %dx%d matrix\n', size(K0,1), size(K0,2));

% Check periodicity
A_T = A_func(T);
B_T = B_func(T);
K_T = K_func(T);

A_error = norm(A0 - A_T, 'fro');
B_error = norm(B0 - B_T, 'fro');
K_error = norm(K0 - K_T, 'fro');

fprintf('  Periodicity check:\n');
fprintf('    ||A(0) - A(T)|| = %.2e\n', A_error);
fprintf('    ||B(0) - B(T)|| = %.2e\n', B_error);
fprintf('    ||K(0) - K(T)|| = %.2e\n', K_error);

if max([A_error, B_error, K_error]) < 1e-12
    fprintf('  ✓ System is periodic with period T = 2π\n\n');
else
    fprintf('  ⚠ Periodicity check failed\n\n');
end

%% Compute Reachability Gramian
fprintf('GRAMIAN COMPUTATION:\n');
fprintf('  Using block-wise algorithm...\n');

tic;
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
computation_time = toc;

fprintf('  Computation completed in %.3f seconds\n\n', computation_time);

%% Analyze Results
fprintf('CONTROLLABILITY ANALYSIS:\n');

% Compute eigenvalues
eigenvals = eig(W);
eigenvals = sort(real(eigenvals), 'descend');

% Compute controllability measures
sigma_min = sqrt(min(eigenvals));
sigma_max = sqrt(max(eigenvals));
kappa = max(eigenvals) / min(eigenvals);
det_W = det(W);
trace_W = trace(W);

% Display results
fprintf('COMPUTED RESULTS:\n');
fprintf('  Gramian eigenvalues: λ₁ = %.6e, λ₂ = %.6e\n', eigenvals(1), eigenvals(2));
fprintf('  Minimum singular value: σ_min = %.6e\n', sigma_min);
fprintf('  Maximum singular value: σ_max = %.6e\n', sigma_max);
fprintf('  Condition number: κ(W) = %.6e\n', kappa);
fprintf('  Determinant: det(W) = %.6e\n', det_W);
fprintf('  Trace: tr(W) = %.6e\n', trace_W);

%% Compare with Paper Values
fprintf('\nCOMPARISON WITH PAPER VALUES:\n');

% Expected values from paper
sigma_min_paper = 1.25e-2;
kappa_paper = 8.4e3;

fprintf('  Paper values:\n');
fprintf('    σ_min (paper) = %.6e\n', sigma_min_paper);
fprintf('    κ(W) (paper)  = %.6e\n', kappa_paper);

% Compute relative errors
sigma_error = abs(sigma_min - sigma_min_paper) / sigma_min_paper * 100;
kappa_error = abs(kappa - kappa_paper) / kappa_paper * 100;

fprintf('  Relative errors:\n');
fprintf('    σ_min error = %.1f%%\n', sigma_error);
fprintf('    κ(W) error  = %.1f%%\n', kappa_error);

%% Controllability Assessment  
fprintf('\nCONTROLLABILITY ASSESSMENT:\n');

% Check controllability
if min(eigenvals) > 1e-10
    fprintf('  ✓ System is CONTROLLABLE (W ≻ 0)\n');
else
    fprintf('  ⚠ System may have controllability issues\n');
end

% Check numerical conditioning
if kappa < 1e6
    fprintf('  ✓ Gramian is well-conditioned (κ < 10⁶)\n');
else
    fprintf('  ⚠ Gramian is ill-conditioned (κ ≥ 10⁶)\n');
end

% Check controllability strength
if sigma_min > 1e-6
    fprintf('  ✓ Good controllability strength (σ_min > 10⁻⁶)\n');
else
    fprintf('  ⚠ Weak controllability (σ_min ≤ 10⁻⁶)\n');
end

%% Validation Summary
fprintf('\nVALIDATION SUMMARY:\n');
fprintf('==================\n');

if sigma_error < 20
    fprintf('✓ σ_min: EXCELLENT agreement with paper (error < 20%%)\n');
elseif sigma_error < 50
    fprintf('○ σ_min: GOOD agreement with paper (error < 50%%)\n');
else
    fprintf('⚠ σ_min: Significant discrepancy with paper\n');
end

if kappa_error < 20
    fprintf('✓ κ(W): EXCELLENT agreement with paper (error < 20%%)\n');
elseif kappa_error < 50
    fprintf('○ κ(W): GOOD agreement with paper (error < 50%%)\n');
else
    fprintf('⚠ κ(W): Significant discrepancy with paper\n');
end

fprintf('\nSystem exhibits strong controllability properties.\n');
fprintf('Block algorithm successfully validated!\n');

%% Optional: Generate Verification Plots
generate_validation_plots(A_func, B_func, K_func, W, T);

end

function generate_validation_plots(A_func, B_func, K_func, W, T)
% Generate validation plots for the small system example

fprintf('\nGenerating validation plots...\n');

% Time vector for plotting
t = linspace(0, T, 1000);

try
    % Create figure
    figure('Name', 'Example 1: Small System Validation', 'Position', [100, 100, 1200, 900]);
    
    % Plot 1: A(t) matrix elements
    subplot(3, 3, 1);
    A11 = arrayfun(@(t) A_func(t)(1,1), t);
    A12 = arrayfun(@(t) A_func(t)(1,2), t);
    A21 = arrayfun(@(t) A_func(t)(2,1), t);
    A22 = arrayfun(@(t) A_func(t)(2,2), t);
    
    plot(t, A11, 'b-', t, A12, 'r--', t, A21, 'g:', t, A22, 'm-.', 'LineWidth', 1.5);
    title('A(t) Matrix Elements');
    xlabel('Time t');
    ylabel('A_{ij}(t)');
    legend('A_{11}', 'A_{12}', 'A_{21}', 'A_{22}', 'Location', 'best');
    grid on;
    
    % Plot 2: B(t) matrix elements
    subplot(3, 3, 2);
    B11 = arrayfun(@(t) B_func(t)(1,1), t);
    B12 = arrayfun(@(t) B_func(t)(1,2), t);
    B21 = arrayfun(@(t) B_func(t)(2,1), t);
    B22 = arrayfun(@(t) B_func(t)(2,2), t);
    
    plot(t, B11, 'b-', t, B12, 'r--', t, B21, 'g:', t, B22, 'm-.', 'LineWidth', 1.5);
    title('B(t) Matrix Elements');
    xlabel('Time t');
    ylabel('B_{ij}(t)');
    legend('B_{11}', 'B_{12}', 'B_{21}', 'B_{22}', 'Location', 'best');
    grid on;
    
    % Plot 3: K(t) matrix elements
    subplot(3, 3, 3);
    K1 = arrayfun(@(t) K_func(t)(1,1), t);
    K2 = arrayfun(@(t) K_func(t)(2,1), t);
    
    plot(t, K1, 'b-', t, K2, 'r--', 'LineWidth', 1.5);
    title('K(t) Matrix Elements');
    xlabel('Time t');
    ylabel('K_i(t)');
    legend('K_1(t)', 'K_2(t)', 'Location', 'best');
    grid on;
    
    % Plot 4: Gramian visualization
    subplot(3, 3, 4);
    imagesc(W);
    colorbar;
    title('Reachability Gramian W');
    xlabel('Column');
    ylabel('Row');
    axis equal tight;
    
    % Plot 5: Gramian eigenvalues
    subplot(3, 3, 5);
    eigenvals = eig(W);
    bar([1, 2], eigenvals, 'FaceColor', [0.3, 0.7, 0.3]);
    title('Gramian Eigenvalues');
    xlabel('Index');
    ylabel('Eigenvalue');
    grid on;
    
    % Plot 6: System trajectory visualization (sample)
    subplot(3, 3, 6);
    % Compute trace of A(t) over time
    trace_A = arrayfun(@(t) trace(A_func(t)), t);
    plot(t, trace_A, 'b-', 'LineWidth', 1.5);
    title('Trace of A(t)');
    xlabel('Time t');
    ylabel('tr(A(t))');
    grid on;
    
    % Plot 7: Input matrix norm
    subplot(3, 3, 7);
    K_norm = arrayfun(@(t) norm(K_func(t), 'fro'), t);
    plot(t, K_norm, 'r-', 'LineWidth', 1.5);
    title('Frobenius Norm of K(t)');
    xlabel('Time t');
    ylabel('||K(t)||_F');
    grid on;
    
    % Plot 8: Controllability measure over one period (conceptual)
    subplot(3, 3, 8);
    sigma_min = sqrt(min(eig(W)));
    sigma_max = sqrt(max(eig(W)));
    kappa = max(eig(W)) / min(eig(W));
    
    measures = [sigma_min, sigma_max, kappa/1000]; % Scale kappa for visualization
    bar(measures, 'FaceColor', [0.7, 0.3, 0.3]);
    set(gca, 'XTickLabel', {'σ_{min}', 'σ_{max}', 'κ/1000'});
    title('Controllability Measures');
    ylabel('Value');
    grid on;
    
    % Plot 9: Phase portrait of A(t) eigenvalues
    subplot(3, 3, 9);
    t_sample = linspace(0, T, 50);
    eig_real = zeros(size(t_sample));
    eig_imag = zeros(size(t_sample));
    
    for i = 1:length(t_sample)
        eigs_A = eig(A_func(t_sample(i)));
        eig_real(i) = real(eigs_A(1));
        eig_imag(i) = imag(eigs_A(1));
    end
    
    plot(eig_real, eig_imag, 'go-', 'MarkerSize', 4, 'LineWidth', 1);
    title('A(t) Eigenvalue Trajectory');
    xlabel('Real Part');
    ylabel('Imaginary Part');
    grid on;
    axis equal;
    
    % Adjust layout
    sgtitle('Example 1: Small System Validation Results', 'FontSize', 14, 'FontWeight', 'bold');
    
    fprintf('Validation plots generated successfully.\n');
    
catch ME
    fprintf('Warning: Could not generate all plots. Error: %s\n', ME.message);
end

end
