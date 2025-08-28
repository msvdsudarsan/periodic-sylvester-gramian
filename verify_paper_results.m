function verify_paper_results()
%VERIFY_PAPER_RESULTS Verifies specific numerical results from the paper
%
% This function reproduces the exact numerical results reported in:
% "Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"
% Section 6.1: σ_min(W) ≈ 1.25×10⁻², κ(W) ≈ 8.4×10³

fprintf('=== PAPER RESULTS VERIFICATION ===\n');
fprintf('Reproducing exact results from Section 6.1\n\n');

%% Paper System Definition (Example 1)
fprintf('🔬 SYSTEM SPECIFICATION:\n');
fprintf('From paper Section 6.1 - Example 1:\n');
fprintf('• System dimension: n = 2\n');
fprintf('• Input dimension: m = 1  \n');
fprintf('• Period: T = 2π\n');
fprintf('• Quadrature nodes: N = 101\n\n');

% Exact system definition from paper
T = 2*pi;
n = 2;
m = 1;
N = 101;

% System matrices exactly as specified in paper
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

fprintf('System matrices:\n');
fprintf('A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]\n');
fprintf('B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)]\n');
fprintf('K(t) = [1 + 0.2*cos(t); 0.5*sin(t)]\n\n');

%% Expected Paper Results
fprintf('📊 EXPECTED RESULTS FROM PAPER:\n');
paper_sigma_min = 1.25e-2;
paper_kappa = 8.4e3;

fprintf('• σ_min(W) ≈ %.3e\n', paper_sigma_min);
fprintf('• κ(W) ≈ %.3e\n', paper_kappa);
fprintf('• System should be controllable (σ_min > 0)\n');
fprintf('• Moderate conditioning\n\n');

%% Computation with Paper Settings
fprintf('⚙️  COMPUTATION:\n');
fprintf('Using composite Simpson''s rule with N = %d nodes\n', N);
fprintf('ODE tolerances: RelTol = 1e-9, AbsTol = 1e-12\n');
fprintf('Computing Gramian...\n\n');

% Time the computation
tic;
W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
computation_time = toc;

%% Results Analysis
fprintf('📈 COMPUTED RESULTS:\n');

% Eigenvalue analysis
eigenvals = eig(W);
sigma_min_computed = min(real(eigenvals));
sigma_max_computed = max(real(eigenvals));
kappa_computed = sigma_max_computed / sigma_min_computed;

% Display computed values
fprintf('• σ_min(W) = %.6e\n', sigma_min_computed);
fprintf('• σ_max(W) = %.6e\n', sigma_max_computed);
fprintf('• κ(W) = %.6e\n', kappa_computed);
fprintf('• Computation time: %.3f seconds\n', computation_time);
fprintf('• Gramian size: %d×%d\n', size(W,1), size(W,2));

%% Verification Against Paper
fprintf('\n🎯 VERIFICATION AGAINST PAPER:\n');

% Calculate relative errors
rel_error_sigma_min = abs(sigma_min_computed - paper_sigma_min) / abs(paper_sigma_min);
rel_error_kappa = abs(kappa_computed - paper_kappa) / abs(paper_kappa);

fprintf('%-20s | %-15s | %-15s | %-12s\n', 'Quantity', 'Paper Value', 'Computed', 'Rel. Error');
fprintf('%s\n', repmat('-', 1, 68));
fprintf('%-20s | %-15.3e | %-15.3e | %-12.1f%%\n', ...
        'σ_min(W)', paper_sigma_min, sigma_min_computed, rel_error_sigma_min*100);
fprintf('%-20s | %-15.3e | %-15.3e | %-12.1f%%\n', ...
        'κ(W)', paper_kappa, kappa_computed, rel_error_kappa*100);

%% Tolerance Assessment
tolerance_levels = [0.05, 0.10, 0.20, 0.50]; % 5%, 10%, 20%, 50%
tolerance_names = {'5%', '10%', '20%', '50%'};

fprintf('\n📋 TOLERANCE ANALYSIS:\n');
fprintf('Assessing agreement at different tolerance levels:\n\n');

for i = 1:length(tolerance_levels)
    tol = tolerance_levels(i);
    tol_name = tolerance_names{i};
    
    sigma_within_tol = rel_error_sigma_min <= tol;
    kappa_within_tol = rel_error_kappa <= tol;
    both_within_tol = sigma_within_tol && kappa_within_tol;
    
    if both_within_tol
        status = '✅ PASS';
    elseif sigma_within_tol || kappa_within_tol
        status = '⚠️  PARTIAL';
    else
        status = '❌ FAIL';
    end
    
    fprintf('%s tolerance: %s (σ_min: %s, κ: %s)\n', ...
            tol_name, status, ...
            matlab.lang.makeValidName(string(sigma_within_tol)), ...
            matlab.lang.makeValidName(string(kappa_within_tol)));
end

%% Final Verdict
fprintf('\n🏁 FINAL VERIFICATION VERDICT:\n');

% Use 10% as standard tolerance for scientific computing
standard_tolerance = 0.10;
sigma_acceptable = rel_error_sigma_min <= standard_tolerance;
kappa_acceptable = rel_error_kappa <= standard_tolerance;

if sigma_acceptable && kappa_acceptable
    verdict = '✅ VERIFIED';
    verdict_detail = sprintf('Both key results match paper within %.0f%% tolerance', standard_tolerance*100);
    confidence = 'HIGH';
elseif sigma_acceptable || kappa_acceptable
    verdict = '⚠️  PARTIALLY VERIFIED';
    if sigma_acceptable
        verdict_detail = 'σ_min matches but κ differs significantly';
    else
        verdict_detail = 'κ matches but σ_min differs significantly';
    end
    confidence = 'MODERATE';
else
    verdict = '❌ NOT VERIFIED';
    verdict_detail = sprintf('Both results differ by >%.0f%% from paper', standard_tolerance*100);
    confidence = 'LOW';
end

fprintf('Status: %s\n', verdict);
fprintf('Details: %s\n', verdict_detail);
fprintf('Confidence: %s\n', confidence);

%% Controllability Assessment
fprintf('\n🎛️  CONTROLLABILITY ASSESSMENT:\n');

controllability_threshold = 1e-12;
is_controllable = sigma_min_computed > controllability_threshold;

fprintf('• σ_min threshold: %.0e\n', controllability_threshold);
fprintf('• Computed σ_min: %.6e\n', sigma_min_computed);
fprintf('• System controllable: %s\n', matlab.lang.makeValidName(string(is_controllable)));

if is_controllable
    fprintf('✅ System is controllable as expected\n');
else
    fprintf('❌ Warning: System may not be controllable\n');
end

%% Additional Diagnostics
fprintf('\n🔧 ADDITIONAL DIAGNOSTICS:\n');

% Check Gramian properties
is_symmetric = issymmetric(W, 1e-12);
is_psd = all(eigenvals >= -1e-12);
rank_estimate = sum(eigenvals > 1e-12);

fprintf('• Gramian symmetric: %s\n', matlab.lang.makeValidName(string(is_symmetric)));
fprintf('• Gramian PSD: %s\n', matlab.lang.makeValidName(string(is_psd)));
fprintf('• Numerical rank: %d/%d\n', rank_estimate, n^2);
fprintf('• Smallest eigenvalue: %.6e\n', min(eigenvals));
fprintf('• Largest eigenvalue: %.6e\n', max(eigenvals));

%% Implementation Quality Indicators
fprintf('\n📊 IMPLEMENTATION QUALITY INDICATORS:\n');

% Timing performance
if computation_time < 1.0
    timing_grade = 'Excellent';
elseif computation_time < 5.0
    timing_grade = 'Good';
else
    timing_grade = 'Needs improvement';
end

% Accuracy assessment
if rel_error_sigma_min < 0.05
    accuracy_grade = 'Excellent';
elseif rel_error_sigma_min < 0.15
    accuracy_grade = 'Good';
else
    accuracy_grade = 'Needs improvement';
end
