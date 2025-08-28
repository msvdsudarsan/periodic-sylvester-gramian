%% COMPREHENSIVE SYSTEM DIAGNOSTIC
% This will check EVERYTHING: A(t), B(t), dimensions, and computation methods
function complete_system_diagnostic()
    clear; clc;
    
    fprintf('\n🔍 COMPREHENSIVE SYSTEM DIAGNOSTIC\n');
    fprintf('=====================================\n\n');
    
    % Paper target values
    sigma_min_paper = 1.25e-02;
    kappa_paper = 8.4e+03;
    
    fprintf('📊 Paper Target Values:\n');
    fprintf('   σ_min = %.3e, κ(W) = %.3e\n\n', sigma_min_paper, kappa_paper);
    
    % Test different system interpretations
    test_cases = {
        'Original A(t), B(t) from paper', @test_original_system;
        'Alternative A(t) interpretation', @test_alternative_A;
        'Different time parametrization', @test_time_param;
        'Matrix scaling investigation', @test_scaling;
        'Single input verification', @test_single_input;
    };
    
    best_error = inf;
    best_case = '';
    
    for i = 1:length(test_cases)
        fprintf('🧪 TEST %d: %s\n', i, test_cases{i,1});
        fprintf('----------------------------------------\n');
        
        try
            [sigma_min, kappa, details] = test_cases{i,2}();
            
            % Calculate errors
            sigma_error = abs(sigma_min - sigma_min_paper) / sigma_min_paper * 100;
            kappa_error = abs(kappa - kappa_paper) / kappa_paper * 100;
            total_error = sigma_error + kappa_error;
            
            fprintf('   Results: σ_min=%.3e, κ=%.3e\n', sigma_min, kappa);
            fprintf('   Errors:  σ_min=%.1f%%, κ=%.1f%%\n', sigma_error, kappa_error);
            fprintf('   Details: %s\n', details);
            
            if total_error < best_error
                best_error = total_error;
                best_case = test_cases{i,1};
            end
            
            if sigma_error < 5 && kappa_error < 5
                fprintf('   🎉 EXCELLENT MATCH!\n');
            elseif sigma_error < 20 && kappa_error < 20
                fprintf('   ✅ GOOD MATCH\n');
            else
                fprintf('   ❌ POOR MATCH\n');
            end
            
        catch ME
            fprintf('   ❌ ERROR: %s\n', ME.message);
        end
        
        fprintf('\n');
    end
    
    fprintf('🏆 BEST RESULT: %s (Total Error: %.1f%%)\n\n', best_case, best_error);
    
    % Additional dimension analysis
    fprintf('🔍 DIMENSION ANALYSIS\n');
    fprintf('====================\n');
    analyze_dimensions();
end

function [sigma_min, kappa, details] = test_original_system()
    % Original system as we understand it
    A = @(t) [0, 1; -1-0.1*cos(t), -0.25*sin(t)];
    B = @(t) [0.5*sin(t); 0.5*cos(t)];
    
    [sigma_min, kappa] = compute_controllability_measures(A, B, [0, 2*pi]);
    
    % Check dimensions
    A_test = A(0);
    B_test = B(0);
    details = sprintf('A: %dx%d, B: %dx%d', size(A_test,1), size(A_test,2), size(B_test,1), size(B_test,2));
end

function [sigma_min, kappa, details] = test_alternative_A()
    % Try alternative A(t) matrix interpretation
    A = @(t) [0, 1; -1, -0.1*cos(t)-0.25*sin(t)]; % Combined damping term
    B = @(t) [0.5*sin(t); 0.5*cos(t)];
    
    [sigma_min, kappa] = compute_controllability_measures(A, B, [0, 2*pi]);
    details = 'Combined damping in A(t)';
end

function [sigma_min, kappa, details] = test_time_param()
    % Try different time parametrization (maybe it's not 0 to 2π?)
    A = @(t) [0, 1; -1-0.1*cos(t), -0.25*sin(t)];
    B = @(t) [0.5*sin(t); 0.5*cos(t)];
    
    [sigma_min, kappa] = compute_controllability_measures(A, B, [0, pi]); % Half period
    details = 'Time interval [0, π] instead of [0, 2π]';
end

function [sigma_min, kappa, details] = test_scaling()
    % Try different scaling factors
    A = @(t) [0, 1; -1-0.1*cos(t), -0.25*sin(t)];
    B = @(t) [0.05*sin(t); 0.05*cos(t)]; % Much smaller B
    
    [sigma_min, kappa] = compute_controllability_measures(A, B, [0, 2*pi]);
    details = 'B(t) scaled by 0.1';
end

function [sigma_min, kappa, details] = test_single_input()
    % Force single input system
    A = @(t) [0, 1; -1-0.1*cos(t), -0.25*sin(t)];
    B = @(t) [0.5*sin(t)]; % Only first component - this should give error
    
    try
        [sigma_min, kappa] = compute_controllability_measures(A, B, [0, 2*pi]);
        details = 'Single input B(t) = [0.5*sin(t)]';
    catch
        % If it fails, try making it 2x1
        B = @(t) [0.5*sin(t); 0];
        [sigma_min, kappa] = compute_controllability_measures(A, B, [0, 2*pi]);
        details = 'B(t) = [0.5*sin(t); 0]';
    end
end

function analyze_dimensions()
    % Analyze what dimensions we should expect
    A = @(t) [0, 1; -1-0.1*cos(t), -0.25*sin(t)];
    B = @(t) [0.5*sin(t); 0.5*cos(t)];
    
    A_test = A(0);
    B_test = B(0);
    
    n = size(A_test, 1); % Number of states
    m = size(B_test, 2); % Number of inputs
    
    fprintf('System dimensions:\n');
    fprintf('   n (states) = %d\n', n);
    fprintf('   m (inputs) = %d\n', m);
    fprintf('   A(t): %dx%d matrix\n', size(A_test,1), size(A_test,2));
    fprintf('   B(t): %dx%d matrix\n', size(B_test,1), size(B_test,2));
    
    if m == 1
        fprintf('   ✅ Single input system (m=1) as expected\n');
    else
        fprintf('   ⚠️  Multi-input system (m=%d) - might be the issue!\n', m);
    end
    
    % Check if B(t) should be interpreted differently
    fprintf('\n🤔 B(t) Interpretation Analysis:\n');
    fprintf('Current B(0) = [%.3f; %.3f] - this gives m=1 input\n', B_test(1), B_test(2));
    fprintf('Paper might mean: B(t) = [0.5*sin(t), 0.5*cos(t)] - this would give m=2 inputs\n');
    
    % Test the m=2 interpretation
    fprintf('\n🧪 Testing B(t) as [0.5*sin(t), 0.5*cos(t)] (2 inputs):\n');
    B_multi = @(t) [0.5*sin(t), 0.5*cos(t); 0, 0]; % First row has 2 inputs
    
    try
        [sigma_min, kappa] = compute_controllability_measures(A, B_multi, [0, 2*pi]);
        sigma_error = abs(sigma_min - 1.25e-02) / 1.25e-02 * 100;
        kappa_error = abs(kappa - 8.4e+03) / 8.4e+03 * 100;
        
        fprintf('   Results: σ_min=%.3e, κ=%.3e\n', sigma_min, kappa);
        fprintf('   Errors:  σ_min=%.1f%%, κ=%.1f%%\n', sigma_error, kappa_error);
        
        if sigma_error < 10 && kappa_error < 10
            fprintf('   🎉 THIS MIGHT BE IT!\n');
        end
    catch ME
        fprintf('   ❌ Error with 2-input interpretation: %s\n', ME.message);
    end
end

function [sigma_min, kappa] = compute_controllability_measures(A_func, B_func, time_interval)
    % Compute controllability Gramian and measures
    n = size(A_func(0), 1);
    m = size(B_func(0), 2);
    
    % Compute Gramian using ODE solver
    W0 = zeros(n*n, 1);
    
    [~, W_vec] = ode45(@(t,W) gramian_ode(t, W, A_func, B_func, n), time_interval, W0);
    
    W_final = reshape(W_vec(end,:), n, n);
    
    % Compute measures
    eigenvals = eig(W_final);
    sigma_min = sqrt(min(eigenvals));
    kappa = max(eigenvals) / min(eigenvals);
end

function dWdt = gramian_ode(t, W_vec, A_func, B_func, n)
    W = reshape(W_vec, n, n);
    A_t = A_func(t);
    B_t = B_func(t);
    
    dWdt_matrix = A_t * W + W * A_t' + B_t * B_t';
    dWdt = reshape(dWdt_matrix, n*n, 1);
end
