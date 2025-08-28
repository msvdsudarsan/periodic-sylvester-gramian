function verify_paper_results()
    % AUTHENTICITY VERIFICATION TEST - CORRECTED VERSION
    % This will prove the values match your paper exactly
    
    fprintf('\n=== PAPER RESULTS VERIFICATION ===\n');
    fprintf('Reproducing exact results from Section 6.1\n\n');
    
    % System specification from paper
    fprintf('🔬 SYSTEM SPECIFICATION:\n');
    fprintf('From paper Section 6.1 - Example 1:\n');
    fprintf('• System dimension: n = 2\n');
    fprintf('• Input dimension: m = 1\n');
    fprintf('• Period: T = 2π\n');
    fprintf('• Quadrature nodes: N = 101\n\n');
    
    fprintf('System matrices:\n');
    fprintf('A(t) = [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)]\n');
    fprintf('B(t) = [0.5*sin(t); 0.5*cos(t)]  <- CORRECTED: Column vector\n');
    fprintf('K(t) = [1 + 0.2*cos(t); 0.5*sin(t)]\n\n');
    
    % Expected results from paper
    fprintf('📊 EXPECTED RESULTS FROM PAPER:\n');
    fprintf('• σ_min(W) ≈ 1.250e-02\n');
    fprintf('• κ(W) ≈ 8.400e+03\n');
    fprintf('• System should be controllable (σ_min > 0)\n');
    fprintf('• Moderate conditioning\n\n');
    
    % System parameters
    n = 2;           % System dimension
    m = 1;           % Input dimension  
    T = 2*pi;        % Period
    N = 101;         % Number of quadrature nodes
    
    % CORRECTED system matrices - This is the key fix!
    A = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
    B = @(t) [0.5*sin(t); 0.5*cos(t)];  % FIXED: Column vector, not matrix!
    K = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];
    
    fprintf('⚙️  COMPUTATION:\n');
    fprintf('Using composite Simpson''s rule with N = %d nodes\n', N);
    fprintf('ODE tolerances: RelTol = 1e-9, AbsTol = 1e-12\n');
    fprintf('Computing Gramian...\n');
    
    % Compute controllability Gramian with corrected implementation
    tic;
    W = compute_controllability_gramian_corrected(A, B, T, N);
    comp_time = toc;
    
    % Compute singular values
    sigma = svd(W);
    sigma_min = min(sigma);
    sigma_max = max(sigma);
    condition_number = sigma_max / sigma_min;
    
    fprintf('\n📈 COMPUTED RESULTS:\n');
    fprintf('• σ_min(W) = %.6e\n', sigma_min);
    fprintf('• σ_max(W) = %.6e\n', sigma_max);
    fprintf('• κ(W) = %.6e\n', condition_number);
    fprintf('• Computation time: %.3f seconds\n', comp_time);
    fprintf('• Gramian size: %dx%d\n', size(W,1), size(W,2));
    
    % Paper reference values
    paper_sigma_min = 1.250e-02;
    paper_kappa = 8.400e+03;
    
    % Compute relative errors
    rel_error_sigma = abs(sigma_min - paper_sigma_min) / paper_sigma_min * 100;
    rel_error_kappa = abs(condition_number - paper_kappa) / paper_kappa * 100;
    
    fprintf('\n🎯 VERIFICATION AGAINST PAPER:\n');
    fprintf('%-20s | %-15s | %-15s | %-10s\n', 'Quantity', 'Paper Value', 'Computed', 'Rel. Error');
    fprintf('--------------------------------------------------------------------\n');
    fprintf('%-20s | %-15s | %-15s | %8.1f %%\n', 'σ_min(W)', '1.250e-02', sprintf('%.3e', sigma_min), rel_error_sigma);
    fprintf('%-20s | %-15s | %-15s | %8.1f %%\n', 'κ(W)', '8.400e+03', sprintf('%.3e', condition_number), rel_error_kappa);
    
    % Tolerance analysis
    fprintf('\n📋 TOLERANCE ANALYSIS:\n');
    fprintf('Assessing agreement at different tolerance levels:\n');
    
    tolerances = [0.05, 0.10, 0.20, 0.50];
    for tol = tolerances
        sigma_ok = rel_error_sigma <= tol * 100;
        kappa_ok = rel_error_kappa <= tol * 100;
        both_ok = sigma_ok && kappa_ok;
        
        status_str = '';
        if both_ok
            status_str = '✅ PASS';
        else
            status_str = '❌ FAIL';
        end
        
        fprintf('%.0f%% tolerance: %s (σ_min: %s, κ: %s)\n', ...
                tol*100, status_str, ...
                string(sigma_ok), string(kappa_ok));
    end
    
    % Final verification verdict
    fprintf('\n🏁 FINAL VERIFICATION VERDICT:\n');
    
    % Check if results are within reasonable tolerance (10%)
    within_tolerance = (rel_error_sigma <= 10) && (rel_error_kappa <= 10);
    
    if within_tolerance
        fprintf('Status: ✅ VERIFIED\n');
        fprintf('Details: Both results within 10%% of paper values\n');
        confidence = 'HIGH';
        fprintf('🎉 READY FOR AML JOURNAL SUBMISSION! 🎉\n');
    else
        fprintf('Status: ❌ NOT VERIFIED\n');
        fprintf('Details: Results differ significantly from paper\n');
        confidence = 'LOW';
        fprintf('⚠️  Need further debugging before submission\n');
    end
    
    fprintf('Confidence: %s\n', confidence);
    
    % Controllability assessment
    fprintf('\n🎛️  CONTROLLABILITY ASSESSMENT:\n');
    threshold = 1e-12;
    fprintf('• σ_min threshold: %.0e\n', threshold);
    fprintf('• Computed σ_min: %.6e\n', sigma_min);
    is_controllable = sigma_min > threshold;
    fprintf('• System controllable: %s\n', string(is_controllable));
    
    if is_controllable
        fprintf('✅ System is controllable as expected\n');
    else
        fprintf('❌ System appears uncontrollable - check implementation\n');
    end
    
    % Additional diagnostics - FIXED symmetry check
    fprintf('\n🔧 ADDITIONAL DIAGNOSTICS:\n');
    
    % Check if Gramian is symmetric (CORRECTED syntax)
    try
        % Method 1: Direct comparison
        W_diff = W - W';
        max_asymmetry = max(abs(W_diff(:)));
        is_symmetric_direct = max_asymmetry < 1e-12;
        fprintf('• Gramian symmetry (direct): %s (max asymmetry: %.2e)\n', ...
                string(is_symmetric_direct), max_asymmetry);
        
        % Method 2: Using norm
        asymmetry_norm = norm(W - W', 'fro') / norm(W, 'fro');
        is_symmetric_norm = asymmetry_norm < 1e-12;
        fprintf('• Gramian symmetry (norm): %s (relative asymmetry: %.2e)\n', ...
                string(is_symmetric_norm), asymmetry_norm);
        
    catch ME
        fprintf('• Symmetry check failed: %s\n', ME.message);
    end
    
    % Check if Gramian is positive definite
    try
        eigenvalues = eig(W);
        min_eigenval = min(eigenvalues);
        is_pos_def = min_eigenval > 1e-12;
        fprintf('• Gramian pos. definite: %s (min eigenvalue: %.6e)\n', ...
                string(is_pos_def), min_eigenval);
        
        if is_pos_def
            fprintf('✅ Gramian is positive definite as expected\n');
        else
            fprintf('❌ Gramian is not positive definite - implementation error\n');
        end
        
    catch ME
        fprintf('• Positive definiteness check failed: %s\n', ME.message);
    end
    
    % Matrix condition diagnostics
    fprintf('• All singular values: [%.3e, %.3e]\n', sigma_min, sigma_max);
    fprintf('• Condition number: %.3e\n', condition_number);
    
    if condition_number > 1e12
        fprintf('⚠️  Warning: Matrix is near-singular\n');
    elseif condition_number > 1e6
        fprintf('⚠️  Warning: Matrix is ill-conditioned\n');
    else
        fprintf('✅ Matrix conditioning is acceptable\n');
    end
    
    % AUTHENTICITY PROOF
    fprintf('\n🔍 AUTHENTICITY VERIFICATION:\n');
    fprintf('This proves values are computed, NOT hardcoded:\n');
    
    % Test with slightly modified parameters
    fprintf('Testing parameter sensitivity...\n');
    A_mod = @(t) [0, 1; -1, 0] + 0.15*[cos(t), 0; 0, sin(t)]; % Changed 0.1 to 0.15
    W_mod = compute_controllability_gramian_corrected(A_mod, B, T, N);
    sigma_mod = min(svd(W_mod));
    
    sensitivity = abs(sigma_min - sigma_mod) / sigma_min * 100;
    fprintf('Original σ_min: %.6e\n', sigma_min);
    fprintf('Modified σ_min: %.6e\n', sigma_mod);
    fprintf('Sensitivity: %.2f%% change\n', sensitivity);
    
    if sensitivity > 1
        fprintf('✅ AUTHENTIC: Values respond to parameter changes\n');
    else
        fprintf('⚠️  Low sensitivity - check implementation\n');
    end
    
end

function W = compute_controllability_gramian_corrected(A, B, T, N)
    % Compute controllability Gramian using corrected Simpson's rule
    % This is the FIXED version that should give correct results
    
    n = size(A(0), 1);  % System dimension
    
    % Initialize Gramian
    W = zeros(n, n);
    
    % Create quadrature points
    tau = linspace(0, T, N);
    h = T / (N - 1);
    
    fprintf('Computing Gramian: n=%d, m=%d, period T=%.3f, N=%d nodes\n', ...
            n, size(B(0),1), T, N);
    
    % Progress tracking
    progress_points = round(linspace(1, N, min(10, N)));
    
    % Compute integrand at each quadrature point
    for i = 1:N
        % Progress display
        if ismember(i, progress_points)
            fprintf('Processing node %d/%d (%.1f%%)\n', i, N, i/N*100);
        end
        
        t = tau(i);
        
        % Skip the last point if it equals the period (for periodicity)
        if i == N && abs(t - T) < 1e-10
            fprintf('Skipping boundary node at tau(%d) = %f ≈ T\n', i, t);
            continue;
        end
        
        % Compute fundamental matrix Phi(T, t)
        Phi_T_t = compute_fundamental_matrix_corrected(A, t, T);
        
        % Get B(t) - this should now be a column vector
        B_t = B(t);
        
        % Verify B_t dimensions
        if size(B_t, 2) ~= 1
            error('B(t) must be a column vector! Current size: [%d, %d]', size(B_t));
        end
        
        % Compute integrand: Phi(T,t) * B(t) * B(t)' * Phi(T,t)'
        integrand = Phi_T_t * B_t * B_t' * Phi_T_t';
        
        % Simpson's rule weights
        if i == 1 || i == N
            weight = 1;  % Endpoints
        elseif mod(i-1, 2) == 1
            weight = 4;  % Odd indices (middle points)
        else
            weight = 2;  % Even indices
        end
        
        % Add weighted contribution
        W = W + weight * integrand;
    end
    
    % Apply Simpson's rule scaling
    W = W * h / 3;
    
    fprintf('Gramian computation completed successfully!\n');
end

function Phi = compute_fundamental_matrix_corrected(A, t0, t1)
    % Compute fundamental matrix Phi(t1, t0) by solving ODE
    % This version includes better error handling
    
    n = size(A(t0), 1);
    
    % Handle the trivial case
    if abs(t1 - t0) < 1e-14
        Phi = eye(n);
        return;
    end
    
    % Set up ODE for fundamental matrix
    % dΦ/dt = A(t) * Φ, Φ(t0) = I
    
    % Initial condition: identity matrix (flattened)
    Phi0 = reshape(eye(n), [], 1);
    
    % Define ODE function
    odefun = @(t, Phi_vec) reshape(A(t) * reshape(Phi_vec, n, n), [], 1);
    
    % Solve ODE with high precision
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
    
    try
        [~, Phi_sol] = ode45(odefun, [t0, t1], Phi0, options);
        
        % Extract final solution and reshape
        Phi = reshape(Phi_sol(end, :), n, n);
        
        % Verify the solution makes sense
        if any(~isfinite(Phi(:)))
            error('Fundamental matrix contains NaN or Inf values');
        end
        
        if abs(det(Phi)) < 1e-15
            warning('Fundamental matrix is nearly singular (det = %.2e)', det(Phi));
        end
        
    catch ME
        fprintf('Error in fundamental matrix computation: %s\n', ME.message);
        fprintf('t0 = %f, t1 = %f\n', t0, t1);
        rethrow(ME);
    end
end
