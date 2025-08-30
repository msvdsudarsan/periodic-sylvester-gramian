function diagnostic_report = diagnose_gramian_error()
    % CRITICAL DIAGNOSTIC TOOL FOR GRAMIAN COMPUTATION ERROR
    % This systematically tests all possible sources of the 10,000% error
    
    fprintf('=== CRITICAL ERROR DIAGNOSTIC REPORT ===\n\n');
    
    %% 1. EXACT SYSTEM FROM PAPER
    fprintf('1. VERIFYING EXACT SYSTEM FROM PAPER:\n');
    
    % System from your TEX paper (Example 1)
    T = 2*pi;
    A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
    B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
    K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];
    
    % Test matrix dimensions and values at t=0
    A0 = A_func(0); B0 = B_func(0); K0 = K_func(0);
    fprintf('   A(0) = \n'); disp(A0);
    fprintf('   B(0) = \n'); disp(B0);
    fprintf('   K(0) = \n'); disp(K0);
    fprintf('   System dimensions: n=%d, m=%d\n\n', size(K0,1), size(K0,2));
    
    %% 2. STANDARD KRONECKER METHOD (GROUND TRUTH)
    fprintf('2. COMPUTING WITH STANDARD KRONECKER METHOD:\n');
    
    N = 101; % Same as paper
    W_standard = compute_gramian_standard_method(A_func, B_func, K_func, T, N);
    
    sigma_min_std = min(eig(W_standard));
    kappa_std = cond(W_standard);
    
    fprintf('   Standard method results:\n');
    fprintf('   σ_min = %.6e\n', sigma_min_std);
    fprintf('   κ(W)  = %.6e\n', kappa_std);
    fprintf('   Paper claims: σ_min = 1.25e-2, κ = 8.4e3\n');
    fprintf('   Error ratios: σ_min %.1fx, κ %.1fx\n\n', ...
        sigma_min_std/(1.25e-2), kappa_std/(8.4e3));
    
    %% 3. YOUR BLOCK METHOD
    fprintf('3. TESTING YOUR BLOCK METHOD:\n');
    
    W_block = compute_gramian_block_method(A_func, B_func, K_func, T, N);
    
    sigma_min_block = min(eig(W_block));
    kappa_block = cond(W_block);
    
    fprintf('   Block method results:\n');
    fprintf('   σ_min = %.6e\n', sigma_min_block);
    fprintf('   κ(W)  = %.6e\n', kappa_block);
    
    % Compare methods
    relative_error = norm(W_standard - W_block, 'fro') / norm(W_standard, 'fro');
    fprintf('   Relative error vs standard: %.6e\n\n', relative_error);
    
    %% 4. ALTERNATIVE VERIFICATION
    fprintf('4. MATLAB BUILT-IN VERIFICATION:\n');
    
    % LTI approximation at t=0
    sys_lti = ss(A0, K0, eye(2), 0);
    try
        W_matlab = gram(sys_lti, 'c');
        sigma_min_matlab = min(eig(W_matlab));
        fprintf('   MATLAB gram(): σ_min = %.6e\n', sigma_min_matlab);
    catch
        fprintf('   MATLAB gram() failed (expected for this system)\n');
    end
    
    %% 5. SCALING AND UNITS CHECK
    fprintf('\n5. CHECKING FOR SCALING ISSUES:\n');
    
    % Check if there are missing factors of T, 2π, etc.
    scales = [1, T, T^2, 2*pi, (2*pi)^2, 1/(2*pi), 1/T];
    scale_names = {'1', 'T', 'T^2', '2π', '(2π)^2', '1/(2π)', '1/T'};
    
    target_sigma = 1.25e-2;
    target_kappa = 8.4e3;
    
    fprintf('   Testing scaling factors:\n');
    for i = 1:length(scales)
        scaled_sigma = sigma_min_std * scales(i);
        scaled_kappa = kappa_std; % Condition number doesn't scale
        
        sigma_error = abs(scaled_sigma - target_sigma) / target_sigma;
        if sigma_error < 0.1 % Within 10%
            fprintf('   ** POTENTIAL MATCH: Scale by %s gives σ_min = %.6e (error %.1f%%)\n', ...
                scale_names{i}, scaled_sigma, sigma_error*100);
        end
    end
    
    %% 6. PARAMETER SENSITIVITY
    fprintf('\n6. PARAMETER SENSITIVITY TEST:\n');
    
    % Test if small changes in coefficients could explain the discrepancy
    coeff_scales = [0.9, 1.0, 1.1];
    for scale = coeff_scales
        A_scaled = @(t) scale * A_func(t);
        W_test = compute_gramian_standard_method(A_scaled, B_func, K_func, T, N);
        sigma_test = min(eig(W_test));
        fprintf('   A scaled by %.1f: σ_min = %.6e\n', scale, sigma_test);
    end
    
    %% 7. FINAL DIAGNOSIS
    fprintf('\n=== FINAL DIAGNOSIS ===\n');
    
    if abs(sigma_min_std - 1.25e-2) / 1.25e-2 < 0.1
        fprintf('✓ Your implementation is CORRECT - paper values verified\n');
    else
        fprintf('✗ CRITICAL ERROR FOUND:\n');
        fprintf('  Computed: σ_min = %.6e, κ = %.6e\n', sigma_min_std, kappa_std);
        fprintf('  Paper:    σ_min = %.6e, κ = %.6e\n', 1.25e-2, 8.4e3);
        fprintf('  Ratio:    σ_min %.0fx off, κ %.0fx off\n', ...
            sigma_min_std/(1.25e-2), kappa_std/(8.4e3));
        
        if relative_error > 1e-6
            fprintf('  Block method also has implementation errors\n');
        end
        
        fprintf('\n  LIKELY CAUSES:\n');
        fprintf('  1. Paper values are theoretical/incorrect\n');
        fprintf('  2. Missing scaling factor in equations\n');
        fprintf('  3. Different system parameters than claimed\n');
        fprintf('  4. Fundamental mathematical error in paper\n');
        
        fprintf('\n  REQUIRED ACTIONS:\n');
        fprintf('  1. Verify against published literature\n');
        fprintf('  2. Re-derive theoretical formulas\n');
        fprintf('  3. Update paper with correct values\n');
        fprintf('  4. Cannot publish with current claimed values\n');
    end
    
    % Return diagnostic data
    diagnostic_report = struct(...
        'W_standard', W_standard, ...
        'W_block', W_block, ...
        'sigma_min_computed', sigma_min_std, ...
        'sigma_min_claimed', 1.25e-2, ...
        'kappa_computed', kappa_std, ...
        'kappa_claimed', 8.4e3, ...
        'relative_error', relative_error);
end

function W = compute_gramian_standard_method(A_func, B_func, K_func, T, N)
    % Standard Kronecker-based computation for verification
    
    K0 = K_func(0);
    [n, m] = size(K0);
    
    % Simpson quadrature
    tau = linspace(0, T, N);
    h = T / (N-1);
    weights = ones(1, N);
    weights(1) = 1/3; weights(end) = 1/3;
    for i = 2:N-1
        if mod(i-1, 2) == 0
            weights(i) = 4/3;
        else
            weights(i) = 2/3;
        end
    end
    weights = weights * h;
    
    W = zeros(n^2, n^2);
    
    for i = 1:N
        % Build Kronecker system matrix
        At = A_func(tau(i));
        Bt = B_func(tau(i));
        Kt = K_func(tau(i));
        
        A_kron = kron(eye(n), At) + kron(Bt.', eye(n));
        K_kron = kron(eye(n), Kt);
        
        % State transition matrix from tau(i) to T
        Phi = expm(A_kron * (T - tau(i)));
        
        % Gramian integrand
        integrand = Phi * K_kron * K_kron.' * Phi.';
        
        W = W + weights(i) * integrand;
    end
end

function W = compute_gramian_block_method(A_func, B_func, K_func, T, N)
    % Your block method implementation
    
    K0 = K_func(0);
    [n, m] = size(K0);
    
    % Simpson quadrature
    tau = linspace(0, T, N);
    h = T / (N-1);
    weights = ones(1, N);
    weights(1) = 1/3; weights(end) = 1/3;
    for i = 2:N-1
        if mod(i-1, 2) == 0
            weights(i) = 4/3;
        else
            weights(i) = 2/3;
        end
    end
    weights = weights * h;
    
    W = zeros(n^2, n^2);
    
    for i = 1:N
        Ki = K_func(tau(i));
        M_i = zeros(n^2, m*n);
        
        % Block propagation
        for k = 1:m
            for j = 1:n
                % Initial condition: Z0 = Ki(:,k) * ej'
                ej = zeros(n, 1); ej(j) = 1;
                Z0 = Ki(:, k) * ej.';
                
                % Solve Sylvester ODE
                sylv_ode = @(t, Z_vec) sylvester_rhs_vectorized(t, Z_vec, A_func, B_func, n);
                [~, Z_sol] = ode45(sylv_ode, [tau(i), T], Z0(:), ...
                    odeset('RelTol', 1e-9, 'AbsTol', 1e-12));
                
                % Extract and assign
                Z_final = reshape(Z_sol(end, :), n, n);
                col_idx = (k-1)*n + j;
                M_i(:, col_idx) = Z_final(:);
            end
        end
        
        W = W + weights(i) * (M_i * M_i.');
    end
end

function dZ_vec = sylvester_rhs_vectorized(t, Z_vec, A_func, B_func, n)
    Z = reshape(Z_vec, n, n);
    dZ = A_func(t) * Z + Z * B_func(t);
    dZ_vec = dZ(:);
end
