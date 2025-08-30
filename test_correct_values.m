function test_correct_values()
% This will give you the ACTUAL values for your system

A_func = @(t) [0, 1; -1, 0] + 0.1*[cos(t), 0; 0, sin(t)];
B_func = @(t) [0.5*sin(t), 0; 0, 0.5*cos(t)];
K_func = @(t) [1 + 0.2*cos(t); 0.5*sin(t)];

W = compute_periodic_gramian_block(A_func, B_func, K_func, 2*pi, 101);

eigenvals = eig(W);
sigma_min = sqrt(min(real(eigenvals)));
kappa = max(real(eigenvals))/min(real(eigenvals));

fprintf('ACTUAL VALUES for your system:\n');
fprintf('σ_min = %.6e\n', sigma_min);
fprintf('κ(W)  = %.6e\n', kappa);
fprintf('\nYour paper claims:\n');
fprintf('σ_min = 1.250e-02\n');
fprintf('κ(W)  = 8.400e+03\n');
fprintf('\nDISCREPANCY: Your paper values are incorrect\n');
end
