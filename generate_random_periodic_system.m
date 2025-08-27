function [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T)
% Generate random periodic coefficient matrices for testing
% Ensures the system has nice stability and periodicity properties
%
% Inputs:
%   n - state dimension  
%   m - input dimension
%   T - period
%
% Outputs:
%   A_func, B_func, K_func - function handles for periodic matrices
%
% Author: M. S. V. D. Sudarsan
% Paper: "Controllability and Efficient Gramian Computation for Periodic Sylvester Matrix Systems"

% Random number generator seed for reproducibility
rng(42);  % Fixed seed for consistent results

% Generate base matrices with appropriate scaling
A0 = randn(n, n) * 0.1;                    % Small perturbation around zero
A1 = randn(n, n) * 0.05;                   % Cosine component
A2 = randn(n, n) * 0.05;                   % Sine component

B0 = randn(n, n) * 0.2;                    % Base B matrix
B1 = randn(n, n) * 0.1;                    % Periodic variation

K0 = randn(n, m);                          % Base input matrix
K1 = randn(n, m) * 0.3;                    % Periodic variation

% Ensure A0 has reasonable eigenvalues (stable but not too stable)
[V, D] = eig(A0);
eigenvals = diag(D);
% Keep real parts negative but not too negative
eigenvals = -0.1 - 0.5*abs(real(eigenvals)) + 1i*imag(eigenvals);
A0 = real(V * diag(eigenvals) * inv(V));

% Make B0 reasonably conditioned
[U, S, V] = svd(B0);
S = diag(0.1 + 0.9*diag(S)/max(diag(S)));  % Scale singular values
B0 = U * S * V';

% Scale input matrix for good controllability
K0_norm = norm(K0, 'fro');
if K0_norm > eps
    K0 = K0 / K0_norm * sqrt(n);  % Reasonable scaling
end

% Define periodic functions with fundamental frequency ω = 2π/T
omega = 2*pi/T;

% A(t) = A0 + A1*cos(ωt) + A2*sin(ωt)
A_func = @(t) A0 + A1*cos(omega*t) + A2*sin(omega*t);

% B(t) = B0 + B1*cos(2ωt)  (different frequency for richness)  
B_func = @(t) B0 + B1*cos(2*omega*t);

% K(t) = K0 + K1*sin(ωt)
K_func = @(t) K0 + K1*sin(omega*t);

% Optional: Add some validation
if nargout == 0
    % Test evaluation at a few points
    fprintf('Generated periodic system with n=%d, m=%d, T=%.2f\n', n, m, T);
    
    t_test = [0, T/4, T/2, 3*T/4];
    fprintf('Sample evaluations:\n');
    
    for i = 1:length(t_test)
        t = t_test(i);
        A_val = A_func(t);
        B_val = B_func(t);  
        K_val = K_func(t);
        
        fprintf('t=%.2f: ||A||=%.3f, ||B||=%.3f, ||K||=%.3f\n', ...
                t, norm(A_val,'fro'), norm(B_val,'fro'), norm(K_val,'fro'));
    end
    
    % Verify periodicity
    fprintf('\nPeriodicity check (should be near zero):\n');
    fprintf('||A(T) - A(0)|| = %.2e\n', norm(A_func(T) - A_func(0), 'fro'));
    fprintf('||B(T) - B(0)|| = %.2e\n', norm(B_func(T) - B_func(0), 'fro')); 
    fprintf('||K(T) - K(0)|| = %.2e\n', norm(K_func(T) - K_func(0), 'fro'));
end

end
