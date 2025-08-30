function salvage_attempt()
    % EMERGENCY SALVAGE: Create a mathematically valid system
    % This could save your paper methodology (but not the current results)
    
    fprintf('=== EMERGENCY SALVAGE ATTEMPT ===\n\n');
    fprintf('Creating STABLE, CONTROLLABLE system to test your methodology...\n\n');
    
    %% CORRECTED STABLE SYSTEM
    fprintf('1. CORRECTED STABLE SYSTEM DESIGN:\n');
    
    T = 2*pi;
    
    % CORRECTED: Ensure stability with damping
    A_func = @(t) [-1, 1; -1, -2] + 0.1*[cos(t), 0; 0, sin(t)];  % Stable base
    B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];                % Same coupling  
    K_func = @(t) [1 + 0.2*cos(t), 0.5*sin(t); 0.5*sin(t), 1 + 0.2*cos(t)]; % Full rank
    
    % Verify stability
    fprintf('   Stability verification:\n');
    t_test = linspace(0, T, 10);
    max_real = -inf;
    
    for i = 1:length(t_test)
        At = A_func(t_test(i));
        Bt = B_func(t_test(i));
        A_kron = kron(eye(2), At) + kron(Bt.', eye(2));
        eigs_t = eig(A_kron);
        max_real = max(max_real, max(real(eigs_t)));
        
        fprintf('     t=%.2f: max Re(λ) = %.3f\n', t_test(i), max(real(eigs_t)));
    end
    
    if max_real < 0
        fprintf('   ✓ SYSTEM IS STABLE (max Re(λ) = %.3f < 0)\n', max_real);
        system_stable = true;
    else
        fprintf('   ✗ System still unstable\n');
        system_stable = false;
    end
    
    % Verify controllability
    fprintf('\n   Controllability verification:\n');
    controllable = true;
    for i = 1:5
        t = (i-1)*T/4;
        Kt = K_func(t);
        rank_K = rank(Kt, 1e-10);
        fprintf('     t=%.2f: rank(K) = %d/%d\n', t, rank_K, size(Kt,1));
        if rank_K < size(Kt,1)
            controllable = false;
        end
    end
    
    if controllable
        fprintf('   ✓ INPUT MATRIX HAS FULL RANK\n');
    else
        fprintf('   ✗ Input matrix still rank deficient\n');
    end
    
    %% COMPUTE CORRECTED GRAMIAN
    if system_stable && controllable
        fprintf('\n2. COMPUTING CORRECTED GRAMIAN:\n');
        
        % Standard method
        W_corrected = compute_gramian_corrected(A_func, B_func, K_func, T, 101);
        
        sigma_min = min(eig(W_corrected));
        kappa = cond(W_corrected);
        
        fprintf('   Corrected system results:\n');
        fprintf('   σ_min = %.6e\n', sigma_min);
        fprintf('   κ(W)  = %.6e\n', kappa);
        fprintf('   System is controllable: %s\n', sigma_min > 1e-10);
        
        % Test block method on corrected system
        fprintf('\n3. TESTING BLOCK METHOD ON CORRECTED SYSTEM:\n');
        
        W_block_corrected = compute_gramian_block_corrected_system(A_func, B_func, K_func, T, 101);
        
        sigma_block = min(eig(W_block_corrected));
        relative_error = norm(W_corrected - W_block_corrected, 'fro') / norm(W_corrected, 'fro');
        
        fprintf('   Block method σ_min = %.6e\n', sigma_block);
        fprintf('   Relative error = %.6e\n', relative_error);
        
        if relative_error < 1e-8
            fprintf('   ✓ Block method works on corrected system\n');
            method_works = true;
        else
            fprintf('   ✗ Block method still has errors\n');
            method_works = false;
        end
        
    else
        fprintf('\n   Cannot compute Gramian - system is invalid\n');
        method_works = false;
        sigma_min = NaN;
        kappa = NaN;
    end
    
    %% SALVAGE ASSESSMENT
    fprintf('\n=== SALVAGE ASSESSMENT ===\n');
    
    if system_stable && controllable && method_works
        fprintf('✓ METHODOLOGY IS SALVAGEABLE\n');
        fprintf('   - Block method works on properly designed systems\n');
        fprintf('   - Your algorithmic approach has merit\n');
        fprintf('   - Problem was the example system choice\n');
        
        fprintf('\n   CORRECTED PAPER STRATEGY:\n');
        fprintf('   1. Replace unstable system with corrected stable system\n');
        fprintf('   2. Update all numerical values with corrected results\n');
        fprintf('   3. Add stability verification to methodology\n');
        fprintf('   4. Acknowledge the correction in revision\n');
        
        fprintf('\n   NEW PAPER VALUES:\n');
        fprintf('   σ_min = %.6e (replaces your 1.25e-2)\n', sigma_min);
        fprintf('   κ(W)  = %.6e (replaces your 8.4e3)\n', kappa);
        
    else
        fprintf('✗ METHODOLOGY IS FUNDAMENTALLY FLAWED\n');
        fprintf('   - Cannot create valid test systems\n');
        fprintf('   - Block method implementation broken\n');
        fprintf('   - Consider abandoning this approach\n');
    end
    
    %% ALTERNATIVE STABLE SYSTEMS
    fprintf('\n4. ALTERNATIVE STABLE SYSTEM EXAMPLES:\n');
    fprintf('   If corrected system above fails, try these:\n\n');
    
    % Example 1: Simple damped oscillator
    fprintf('   EXAMPLE 1: Damped harmonic oscillator\n');
    fprintf('   A(t) = [-0.5, 1; -1, -0.5] + 0.1*[cos(t), 0; 0, sin(t)]\n');
    fprintf('   B(t) = 0.1*[sin(t), 0; 0, cos(t)]\n');
    fprintf('   K(t) = [1; 1] (constant full-rank input)\n\n');
    
    % Example 2: Decoupled stable system  
    fprintf('   EXAMPLE 2: Decoupled stable system\n');
    fprintf('   A(t) = [-2, 0; 0, -3] + 0.2*[cos(2*t), sin(t); -sin(t), cos(2*t)]\n');
    fprintf('   B(t) = zeros(2,2) (no coupling)\n');
    fprintf('   K(t) = [1+0.1*cos(t); 1+0.1*sin(t)]\n\n');
    
    %% FINAL RECOMMENDATION
    fprintf('=== FINAL EMERGENCY RECOMMENDATION ===\n');
    
    if system_stable && controllable && method_works
        fprintf('PROCEED WITH SALVAGE:\n');
        fprintf('1. Use corrected stable system in paper\n');
        fprintf('2. Update all numerical results\n');
        fprintf('3. Add mathematical validity checks\n');
        fprintf('4. Submit as major revision\n');
    else
        fprintf('ABANDON CURRENT APPROACH:\n');
        fprintf('1. Your original system is mathematically invalid\n');
        fprintf('2. Block method may have fundamental flaws\n');  
        fprintf('3. Consider different research direction\n');
        fprintf('4. Consult control theory expert before proceeding\n');
    end
    
    fprintf('\nTIME REMAINING: Your academic reputation depends on handling this correctly.\n');
    
end

function W = compute_gramian_corrected(A_func, B_func, K_func, T, N)
    % Standard method for corrected system
    
    K0 = K_func(0);
    [n, m] = size(K0);
    
    % Simpson quadrature
    if mod(N,2) == 0, N = N+1; end
    tau = linspace(0, T, N);
    h = T/(N-1);
    
    weights = ones(1,N) * h/3;
    weights(1) = h/3; weights(end) = h/3;
    for i = 2:N-1
        if mod(i-1,2) == 0
            weights(i) = 4*h/3;
        else  
            weights(i) = 2*h/3;
        end
    end
    
    W = zeros(n^2, n^2);
    
    for i = 1:N
        At = A_func(tau(i));
        Bt = B_func(tau(i)); 
        Kt = K_func(tau(i));
        
        A_kron = kron(eye(n), At) + kron(Bt.', eye(n));
        K_kron = kron(eye(n), Kt);
        
        % State transition matrix
        Phi = expm(A_kron * (T - tau(i)));
        
        % Gramian integrand  
        integrand = Phi * K_kron * K_kron.' * Phi.';
        W = W + weights(i) * integrand;
    end
end

function W = compute_gramian_block_corrected_system(A_func, B_func, K_func, T, N)
    % Block method for corrected system
    
    K0 = K_func(0);
    [n, m] = size(K0);
    
    % Simpson quadrature
    if mod(N,2) == 0, N = N+1; end
    tau = linspace(0, T, N);
    h = T/(N-1);
    
    weights = ones(1,N) * h/3;
    weights(1) = h/3; weights(end) = h/3;
    for i = 2:N-1
        if mod(i-1,2) == 0
            weights(i) = 4*h/3;
        else
            weights(i) = 2*h/3;
        end
    end
    
    W = zeros(n^2, n^2);
    
    for i = 1:N
        Ki = K_func(tau(i));
        M_i = zeros(n^2, m*n);
        
        for k = 1:m
            for j = 1:n
                ej = zeros(n,1); ej(j) = 1;
                Z0 = Ki(:, k) * ej.';
                
                if abs(tau(i) - T) < 1e-12
                    Z_final = Z0;
                else
                    % Solve Sylvester ODE
                    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-15);
                    [~, Z_sol] = ode45(@(t, Z_vec) sylv_corrected_rhs(t, Z_vec, A_func, B_func, n), ...
                        [tau(i), T], Z0(:), opts);
                    Z_final = reshape(Z_sol(end,:), n, n);
                end
                
                col_idx = (k-1)*n + j;  
                M_i(:, col_idx) = Z_final(:);
            end
        end
        
        W = W + weights(i) * (M_i * M_i.');
    end
end

function dZ = sylv_corrected_rhs(t, Z_vec, A_func, B_func, n)
    Z = reshape(Z_vec, n, n);
    At = A_func(t);
    Bt = B_func(t);
    dZ = At * Z + Z * Bt;
    dZ = dZ(:);
end
