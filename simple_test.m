function simple_test()
% SIMPLE_TEST - Quick test to verify the fixes work

clear; clc;
fprintf('Testing the fixed implementation...\n\n');

% Define system from TEX paper
A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

% Quick test with smaller N
W = compute_periodic_gramian_block(A_func, B_func, K_func, 2*pi, 21);

eigenvals = eig(W);
sigma_min = sqrt(min(real(eigenvals)));
kappa = max(real(eigenvals))/min(real(eigenvals));

fprintf('\nQuick test results:\n');
fprintf('σ_min = %.6e\n', sigma_min);
fprintf('κ(W)  = %.6e\n', kappa);

% Expected from paper
fprintf('\nPaper values:\n');
fprintf('σ_min ≈ 1.25e-02\n');
fprintf('κ(W)  ≈ 8.4e+03\n');

if abs(log10(sigma_min) - log10(1.25e-2)) < 1
    fprintf('\n✓ σ_min is in correct order of magnitude\n');
else
    fprintf('\n✗ σ_min order of magnitude incorrect\n');
end

if abs(log10(kappa) - log10(8.4e3)) < 1
    fprintf('✓ κ(W) is in correct order of magnitude\n');
else
    fprintf('✗ κ(W) order of magnitude incorrect\n');
end

fprintf('\nTest completed.\n');

end
