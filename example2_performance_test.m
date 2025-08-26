function example2_performance_test()
% EXAMPLE2_PERFORMANCE_TEST Performance comparison for different dimensions
%
% Author: M. S. V. D. Sudarsan

    fprintf('=== Example 2: Performance Test ===\n');
    
    % Test dimensions
    n_values = [3, 5, 8, 10];
    m = 2;
    T = 2*pi;
    N = 51; % Moderate number of nodes
    
    fprintf('Testing dimensions n = [%s], m = %d\n', ...
            sprintf('%d ', n_values), m);
    fprintf('Period T = %.4f, Quadrature nodes N = %d\n\n', T, N);
    
    results = struct();
    
    for idx = 1:length(n_values)
        n = n_values(idx);
        fprintf('--- Testing n = %d ---\n', n);
        
        % Generate random periodic system
        [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T);
        
        % Time the computation
        tic;
        W = compute_periodic_gramian(A_func, B_func, K_func, T, N);
        comp_time = toc;
        
        % Analyze results
        eigenvals = eig(W);
        sigma_min = min(real(eigenvals));
        sigma_max = max(real(eigenvals));
        
        % Store results
        results(idx).n = n;
        results(idx).time = comp_time;
        results(idx).sigma_min = sigma_min;
        results(idx).condition = sigma_max / sigma_min;
        results(idx).memory_usage = n^4 * 8 / (1024^2); % Approximate MB
        
        fprintf('  Time: %.4f sec, σ_min: %.2e, Condition: %.2e\n', ...
                comp_time, sigma_min, sigma_max/sigma_min);
    end
    
    % Summary table
    fprintf('\n=== Performance Summary ===\n');
    fprintf('n\tTime(s)\tσ_min\t\tCondition\tMemory(MB)\n');
    fprintf('---------------------------------------------------\n');
    for idx = 1:length(results)
        fprintf('%d\t%.3f\t%.2e\t%.2e\t%.1f\n', ...
                results(idx).n, results(idx).time, ...
                results(idx).sigma_min, results(idx).condition, ...
                results(idx).memory_usage);
    end
end
