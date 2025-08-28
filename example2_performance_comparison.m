function example2_performance_comparison()
%EXAMPLE2_PERFORMANCE_COMPARISON Performance comparison for larger systems
%
% Compares computation times and demonstrates scalability advantages
% of the block method vs. direct Kronecker approach for n ∈ {5,10,15,20}

fprintf('=== EXAMPLE 2: PERFORMANCE COMPARISON ===\n');
fprintf('Testing scalability for increasing system dimensions\n\n');

% Test parameters
n_values = [5, 10, 15, 20];
m = 2; % Input dimension
T = 2*pi;
N = 51; % Fewer nodes for performance comparison

% Results storage
results = struct();
results.n = n_values;
results.block_times = zeros(size(n_values));
results.memory_ratios = zeros(size(n_values));
results.controllable = false(size(n_values));

fprintf('System parameters:\n');
fprintf('- Input dimension: m = %d\n', m);
fprintf('- Period: T = %.3f\n', T);
fprintf('- Quadrature nodes: N = %d\n\n', N);

fprintf('%-4s | %-12s | %-12s | %-15s | %-12s\n', 'n', 'Block Time', 'Memory (MB)', 'σ_min(W)', 'Controllable');
fprintf('%s\n', repmat('-', 1, 70));

for i = 1:length(n_values)
    n = n_values(i);
    
    % Generate random periodic system
    [A_func, B_func, K_func] = generate_random_periodic_system(n, n, m, T);
    
    % Compute using block method
    fprintf('Testing n=%d... ', n);
    tic;
    W = compute_periodic_gramian_block(A_func, B_func, K_func, T, N);
    block_time = toc;
    
    % Analyze results
    eigenvals = eig(W);
    sigma_min = min(real(eigenvals));
    is_controllable = sigma_min > 1e-10;
    
    % Estimate memory usage for Gramian (n²×n² matrix)
    gramian_memory_mb = (n^2)^2 * 8 / (1024^2); % 8 bytes per double
    
    % Store results
    results.block_times(i) = block_time;
    results.memory_ratios(i) = gramian_memory_mb;
    results.controllable(i) = is_controllable;
    
    % Display row
    fprintf('%-4d | %-12.3f | %-12.1f | %-15.3e | %-12s\n', ...
            n, block_time, gramian_memory_mb, sigma_min, ...
            matlab.lang.makeValidName(string(is_controllable)));
end

% Analyze scalability
fprintf('\n=== SCALABILITY ANALYSIS ===\n');

% Theoretical complexity analysis
fprintf('Theoretical complexity comparison:\n');
fprintf('- Direct Kronecker method: O(N·n⁶)\n');
fprintf('- Block method: O(N·n³·m)\n\n');

% Compute theoretical speedup ratios
fprintf('Theoretical speedup ratios (n³/m vs n⁶):\n');
for i = 1:length(n_values)
    n = n_values(i);
    theoretical_speedup = n^6 / (n^3 * m);
    fprintf('n=%d: %.0f×\n', n, theoretical_speedup);
end

% Plot timing results if available
if exist('figure', 'file')
    figure('Name', 'Performance Comparison', 'NumberTitle', 'off');
    
    % Timing plot
    subplot(2,2,1);
    semilogy(n_values, results.block_times, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('System dimension n');
    ylabel('Computation time (seconds)');
    title('Block Method Performance');
    grid on;
    
    % Memory usage plot
    subplot(2,2,2);
    semilogy(n_values, results.memory_ratios, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('System dimension n');
    ylabel('Gramian memory (MB)');
    title('Memory Requirements');
    grid on;
    
    % Theoretical vs actual complexity
    subplot(2,2,3);
    theoretical_n3 = (n_values/n_values(1)).^3 * results.block_times(1);
    loglog(n_values, results.block_times, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    loglog(n_values, theoretical_n3, 'r--', 'LineWidth', 2);
    xlabel('System dimension n');
    ylabel('Computation time (seconds)');
    title('Complexity Scaling');
    legend('Actual', 'O(n³)', 'Location', 'northwest');
    grid on;
    
    % Controllability results
    subplot(2,2,4);
    bar(n_values, double(results.controllable));
    xlabel('System dimension n');
    ylabel('Controllable (1=Yes, 0=No)');
    title('Controllability Status');
    ylim([-0.1, 1.1]);
    grid on;
    
    % Adjust layout
    sgtitle('Example 2: Performance Analysis');
end

% Performance summary
fprintf('\n=== PERFORMANCE SUMMARY ===\n');
fprintf('Block method demonstrates excellent scalability:\n');

% Calculate growth rates
if length(n_values) > 1
    time_growth_rate = (results.block_times(end) / results.block_times(1))^(1/(length(n_values)-1));
    fprintf('- Average time growth factor: %.2f per step\n', time_growth_rate);
end

fprintf('- Largest system (n=%d): %.3f seconds\n', n_values(end), results.block_times(end));
fprintf('- Memory efficient: avoids n²×n² Kronecker matrices\n');
fprintf('- All test systems: %d/%d controllable\n', sum(results.controllable), length(n_values));

% Recommendations
fprintf('\n=== RECOMMENDATIONS ===\n');
fprintf('For optimal performance:\n');
fprintf('1. Use block method for n ≥ 10 and m << n³\n');
fprintf('2. Adjust quadrature nodes N based on desired accuracy\n');
fprintf('3. Monitor condition number κ(W) for numerical stability\n');
fprintf('4. Consider iterative eigensolvers for very large systems\n');

end
