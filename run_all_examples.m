%% Master Script: Run All Examples and Tests
% This script reproduces all numerical results from the paper:
% "Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com
% Paper submitted to Applied Mathematics Letters

clear; clc; close all;

fprintf('============================================================\n');
fprintf('PERIODIC SYLVESTER MATRIX SYSTEMS - NUMERICAL VALIDATION\n');
fprintf('============================================================\n');
fprintf('Paper: "Controllability and Efficient Gramian Computation\n');
fprintf('        for Periodic Sylvester Matrix Systems"\n');
fprintf('Author: M. S. V. D. Sudarsan\n');
fprintf('Journal: Applied Mathematics Letters (Submitted)\n');
fprintf('============================================================\n\n');

% Add current directory to path (in case of subdirectories)
addpath(pwd);

% Check for required functions
required_functions = {
    'compute_periodic_gramian_block.m',
    'generate_random_periodic_system.m'
};

fprintf('Checking for required files...\n');
missing_files = {};
for i = 1:length(required_functions)
    if ~exist(required_functions{i}, 'file')
        missing_files{end+1} = required_functions{i};
    else
        fprintf('  ✓ %s\n', required_functions{i});
    end
end

if ~isempty(missing_files)
    fprintf('\nERROR: Missing required files:\n');
    for i = 1:length(missing_files)
        fprintf('  ✗ %s\n', missing_files{i});
    end
    fprintf('Please ensure all files are in the current directory.\n');
    return;
end

fprintf('All required files found.\n\n');

% Menu for user selection
fprintf('Select which examples to run:\n');
fprintf('  1. Example 1: Small System Validation (Section 6.1)\n');
fprintf('  2. Example 2: Performance Comparison (Section 6.2) \n');
fprintf('  3. Convergence Analysis (Section 6.2, Figure 1)\n');
fprintf('  4. Robustness Test (Section 6.3)\n');
fprintf('  5. Run ALL examples\n');
fprintf('  0. Exit\n');

choice = input('Enter your choice (0-5): ');

switch choice
    case 1
        run_example1();
    case 2  
        run_example2();
    case 3
        run_convergence();
    case 4
        run_robustness();
    case 5
        run_all_examples();
    case 0
        fprintf('Exiting...\n');
        return;
    otherwise
        fprintf('Invalid choice. Exiting...\n');
        return;
end

fprintf('\n============================================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('============================================================\n');

%% Individual example functions

function run_example1()
    fprintf('Running Example 1: Small System Validation...\n');
    fprintf('Expected runtime: ~10 seconds\n\n');
    
    try
        example1_small_system_validation;
        fprintf('\n✓ Example 1 completed successfully\n');
    catch ME
        fprintf('\n✗ Example 1 failed: %s\n', ME.message);
    end
    
    input('\nPress Enter to continue...');
end

function run_example2()
    fprintf('Running Example 2: Performance Comparison...\n');
    fprintf('Expected runtime: ~2-5 minutes (depends on system size)\n');
    fprintf('Warning: Large systems (n≥15) may take significant time\n\n');
    
    proceed = input('Continue? (y/n): ', 's');
    if lower(proceed) ~= 'y'
        fprintf('Example 2 skipped.\n');
        return;
    end
    
    try
        example2_performance_comparison;
        fprintf('\n✓ Example 2 completed successfully\n');
    catch ME
        fprintf('\n✗ Example 2 failed: %s\n', ME.message);
    end
    
    input('\nPress Enter to continue...');
end

function run_convergence()
    fprintf('Running Convergence Analysis...\n');
    fprintf('Expected runtime: ~30 seconds\n\n');
    
    try
        convergence_analysis;
        fprintf('\n✓ Convergence analysis completed successfully\n');
    catch ME
        fprintf('\n✗ Convergence analysis failed: %s\n', ME.message);
    end
    
    input('\nPress Enter to continue...');
end

function run_robustness()
    fprintf('Running Robustness Test...\n');
    fprintf('Expected runtime: ~30 seconds\n\n');
    
    try
        robustness_test;
        fprintf('\n✓ Robustness test completed successfully\n');
    catch ME
        fprintf('\n✗ Robustness test failed: %s\n', ME.message);
    end
    
    input('\nPress Enter to continue...');
end

function run_all_examples()
    fprintf('Running ALL examples...\n');
    fprintf('Total expected runtime: ~5-10 minutes\n\n');
    
    proceed = input('This will run all tests. Continue? (y/n): ', 's');
    if lower(proceed) ~= 'y'
        fprintf('All examples skipped.\n');
        return;
    end
    
    % Run all examples in sequence
    fprintf('\n' + string(repmat('=', 1, 60)) + '\n');
    fprintf('STARTING COMPREHENSIVE NUMERICAL VALIDATION\n');
    fprintf(string(repmat('=', 1, 60)) + '\n\n');
    
    example_status = zeros(4, 1);  % Track success/failure
    
    % Example 1
    fprintf('>>> EXAMPLE 1: Small System Validation <<<\n');
    try
        example1_small_system_validation;
        example_status(1) = 1;
        fprintf('✓ Example 1: SUCCESS\n\n');
    catch ME
        fprintf('✗ Example 1: FAILED - %s\n\n', ME.message);
    end
    
    % Example 2
    fprintf('>>> EXAMPLE 2: Performance Comparison <<<\n');
    try
        example2_performance_comparison;
        example_status(2) = 1;
        fprintf('✓ Example 2: SUCCESS\n\n');
    catch ME
        fprintf('✗ Example 2: FAILED - %s\n\n', ME.message);
    end
    
    % Convergence Analysis
    fprintf('>>> CONVERGENCE ANALYSIS <<<\n');
    try
        convergence_analysis;
        example_status(3) = 1;
        fprintf('✓ Convergence Analysis: SUCCESS\n\n');
    catch ME
        fprintf('✗ Convergence Analysis: FAILED - %s\n\n', ME.message);
    end
    
    % Robustness Test
    fprintf('>>> ROBUSTNESS TEST <<<\n');
    try
        robustness_test;
        example_status(4) = 1;
        fprintf('✓ Robustness Test: SUCCESS\n\n');
    catch ME
        fprintf('✗ Robustness Test: FAILED - %s\n\n', ME.message);
    end
    
    % Final summary
    fprintf(string(repmat('=', 1, 60)) + '\n');
    fprintf('FINAL SUMMARY\n');
    fprintf(string(repmat('=', 1, 60)) + '\n');
    
    example_names = {
        'Example 1: Small System Validation',
        'Example 2: Performance Comparison', 
        'Convergence Analysis',
        'Robustness Test'
    };
    
    fprintf('Test Results:\n');
    for i = 1:4
        if example_status(i)
            fprintf('  ✓ %s\n', example_names{i});
        else
            fprintf('  ✗ %s\n', example_names{i});
        end
    end
    
    success_rate = sum(example_status) / length(example_status) * 100;
    fprintf('\nOverall Success Rate: %.0f%% (%d/%d)\n', ...
            success_rate, sum(example_status), length(example_status));
    
    if success_rate == 100
        fprintf('\n🎉 All tests passed! The implementation is validated.\n');
    elseif success_rate >= 75
        fprintf('\n⚠️  Most tests passed. Check failed tests for issues.\n');
    else
        fprintf('\n❌ Multiple tests failed. Please check your setup.\n');
    end
end
