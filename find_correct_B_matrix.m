function find_correct_B_matrix()
    % DIAGNOSTIC: Test different B(t) matrix interpretations
    % to find which one gives paper results
    
    fprintf('\n🔍 B(t) MATRIX DIAGNOSTIC ANALYSIS\n');
    fprintf('===================================\n\n');
    
    % Fixed system matrices
    A = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
    K = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];
    
    % Expected paper results
    expected_sigma_min = 1.250e-02;
    expected_kappa = 8.400e+03;
    
    % System parameters
    T = 2*pi;
    N = 101;
    
    % Test different B(t) interpretations
    fprintf('Testing different B(t) matrix interpretations:\n\n');
    
    % Option 1: Current interpretation (your result)
    fprintf('1️⃣ OPTION 1: B(t) = [0.5*sin(t); 0.5*cos(t)] (Current)\n');
    B1 = @(t) [0.5*sin(t); 0.5*cos(t)];
    test_B_matrix(A, B1, K, T, N, 1, expected_sigma_min, expected_kappa);
    
    % Option 2: Different order
    fprintf('\n2️⃣ OPTION 2: B(t) = [0.5*cos(t); 0.5*sin(t)] (Swapped)\n');
    B2 = @(t) [0.5*cos(t); 0.5*sin(t)];
    test_B_matrix(A, B2, K, T, N, 2, expected_sigma_min, expected_kappa);
    
    % Option 3: Only first component
    fprintf('\n3️⃣ OPTION 3: B(t) = [0.5*sin(t); 0] (First component only)\n');
    B3 = @(t) [0.5*sin(t); 0];
    test_B_matrix(A, B3, K, T, N, 3, expected_sigma_min, expected_kappa);
    
    % Option 4: Only second component
    fprintf('\n4️⃣ OPTION 4: B(t) = [0; 0.5*cos(t)] (Second component only)\n');
    B4 = @(t) [0; 0.5*cos(t)];
    test_B_matrix(A, B4, K, T, N, 4, expected_sigma_min, expected_kappa);
    
    % Option 5: Different scaling
    fprintf('\n5️⃣ OPTION 5: B(t) = [sin(t); cos(t)] * 0.1 (Different scaling)\n');
    B5 = @(t) [sin(t); cos(t)] * 0.1;
    test_B_matrix(A, B5, K, T, N, 5, expected_sigma_min, expected_kappa);
    
    % Option 6: Simple constant
    fprintf('\n6️⃣ OPTION 6: B(t) = [1; 0] (Constant)\n');
    B6 = @(t) [1; 0];
    test_B_matrix(A, B6, K, T, N, 6, expected_sigma_min, expected_kappa);
    
    % Option 7: From paper's exact description - try interpreting as matrix
    fprintf('\n7️⃣ OPTION 7: B(t) = [0.5*sin(t), 0; 0, 0.5*cos(t)] but take first column\n');
    B7 = @(t) [0.5*sin(t); 0];  % First column of the matrix form
    test_B_matrix(A, B7, K, T, N, 7, expected_sigma_min, expected_kappa);
    
    % Option 8: Alternative from K(t) structure
    fprintf('\n8️⃣ OPTION 8: B(t) = [1 + 0.1*cos(t); 0.25*sin(t)] (Similar to K)\n');
    B8 = @(t) [1 + 0.1*cos(t); 0.25*sin(t)];
    test_B_matrix(A, B8, K, T, N, 8, expected_sigma_min, expected_kappa);
    
    % Option 9: Very simple
    fprintf('\n9️⃣ OPTION 9: B(t) = [0; 1] (Simple constant)\n');
    B9 = @(t) [0; 1];
    test_B_matrix(A, B9, K, T, N, 9, expected_sigma_min, expected_kappa);
    
    % Option 10: Try much smaller scaling
    fprintf('\n🔟 OPTION 10: B(t) = [0.05*sin(t); 0.05*cos(t)] (Smaller scale)\n');
    B10 = @(t) [0.05*sin(t); 0.05*cos(t)];
    test_B_matrix(A, B10, K, T, N, 10, expected_sigma_min, expected_kappa);
    
    fprintf('\n🎯 ANALYSIS COMPLETE\n');
    fprintf('Look for the option with lowest error percentage!\n\n');
    
end

function test_B_matrix(A, B, K, T, N, option_num, expected_sigma_min, expected_kappa)
    % Test a specific B(t) matrix interpretation
    
    try
        % Check dimensions
        B_test = B(0);
        [n, m] = size(B_test);
        
        if m ~= 1
            fprintf('   ❌ SKIPPED: B(t) is not a column vector (size: [%d,%d])\n', n, m);
            return;
        end
        
        fprintf('   B(0) = [%.3f; %.3f], size: [%d,%d]\n', B_test(1), B_test(2), n, m);
        
        % Compute Gramian
        W = compute_controllability_gramian_silent(A, B, T, N);
        
        % Compute results
        sigma = svd(W);
        sigma_min = min(sigma);
        sigma_max = max(sigma);
        kappa = sigma_max / sigma_min;
        
        % Compute errors
        error_sigma = abs(sigma_min - expected_sigma_min) / expected_sigma_min * 100;
        error_kappa = abs(kappa - expected_kappa) / expected_kappa * 100;
        
        fprintf('   Results: σ_min=%.3e, κ=%.3e\n', sigma_min, kappa);
        fprintf('   Errors:  σ_min=%.1f%%, κ=%.1f%%\n', error_sigma, error_kappa);
        
        % Grade this option
        if error_sigma < 10 && error_kappa < 10
            fprintf('   🎉 EXCELLENT MATCH! This might be the correct B(t)\n');
        elseif error_sigma < 20 && error_kappa < 20
            fprintf('   ✅ GOOD MATCH! Close to paper values\n');
        elseif error_sigma < 50 && error_kappa < 50
            fprintf('   ⚠️  MODERATE MATCH\n');
        else
            fprintf('   ❌ POOR MATCH\n');
        end
        
    catch ME
        fprintf('   ❌ ERROR: %s\n', ME.message);
    end
end

function W = compute_controllability_gramian_silent(A, B, T, N)
    % Silent version of Gramian computation (no progress output)
    
    n = size(A(0), 1);
    W = zeros(n, n);
    tau = linspace(0, T, N);
    h = T / (N - 1);
    
    for i = 1:N
        t = tau(i);
        
        if i == N && abs(t - T) < 1e-10
            continue;
        end
        
        Phi_T_t = compute_fundamental_matrix_silent(A, t, T);
        B_t = B(t);
        integrand = Phi_T_t * B_t * B_t' * Phi_T_t';
        
        if i == 1 || i == N
            weight = 1;
        elseif mod(i-1, 2) == 1
            weight = 4;
        else
            weight = 2;
        end
        
        W = W + weight * integrand;
    end
    
    W = W * h / 3;
end

function Phi = compute_fundamental_matrix_silent(A, t0, t1)
    % Silent version of fundamental matrix computation
    
    n = size(A(t0), 1);
    
    if abs(t1 - t0) < 1e-14
        Phi = eye(n);
        return;
    end
    
    Phi0 = reshape(eye(n), [], 1);
    odefun = @(t, Phi_vec) reshape(A(t) * reshape(Phi_vec, n, n), [], 1);
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
    
    [~, Phi_sol] = ode45(odefun, [t0, t1], Phi0, options);
    Phi = reshape(Phi_sol(end, :), n, n);
end
