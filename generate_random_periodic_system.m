function [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T)
% GENERATE_RANDOM_PERIODIC_SYSTEM Generate random periodic system matrices
%
% This function generates random periodic matrices A(t), B(t), K(t) for
% testing the block Gramian computation algorithm. The matrices are 
% constructed to be smooth, periodic, and ensure numerical stability.
%
% SYNTAX:
%   [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T)
%
% INPUTS:
%   n - State dimension
%   m - Input dimension  
%   T - Period length
%
% OUTPUTS:
%   A_func - Function handle for A(t) matrix (n x n)
%   B_func - Function handle for B(t) matrix (n x n)
%   K_func - Function handle for K(t) matrix (n x m)
%
% DESIGN PRINCIPLES:
%   - A(t) is constructed to be stable on average
%   - B(t) has controlled spectral properties
%   - K(t) ensures controllability for most t
%   - All matrices are smooth and T-periodic
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com

% Input validation
if nargin < 3
    error('Three inputs required: n, m, T');
end

if n < 1 || m < 1 || T <= 0
    error('Dimensions must be positive and T > 0');
end

fprintf('Generating random periodic system: n=%d, m=%d, T=%.2f\n', n, m, T);

%% Generate A(t) matrix
% A(t) = A0 + A1*cos(2πt/T) + A2*sin(2πt/T) + A3*cos(4πt/T) + A4*sin(4πt/T)

% Base matrix A0 - make it stable
A0 = generate_stable_matrix(n);

% Perturbation matrices with controlled magnitudes
perturbation_scale = 0.1;  % Small perturbations to maintain stability
A1 = perturbation_scale * randn(n, n);
A2 = perturbation_scale * randn(n, n);
A3 = perturbation_scale * 0.5 * randn(n, n);  % Smaller higher harmonics
A4 = perturbation_scale * 0.5 * randn(n, n);

% Make perturbations have zero trace to preserve stability better
A1 = A1 - (trace(A1)/n) * eye(n);
A2 = A2 - (trace(A2)/n) * eye(n);
A3 = A3 - (trace(A3)/n) * eye(n);
A4 = A4 - (trace(A4)/n) * eye(n);

% A(t) function handle
A_func = @(t) A0 + A1*cos(2*pi*t/T) + A2*sin(2*pi*t/T) + ...
              A3*cos(4*pi*t/T) + A4*sin(4*pi*t/T);

%% Generate B(t) matrix
% B(t) = B0 + B1*cos(2πt/T) + B2*sin(2πt/T)

% Base matrix B0
B0 = 0.5 * randn(n, n);

% Perturbation matrices
B_perturbation_scale = 0.3;
B1 = B_perturbation_scale * randn(n, n);
B2 = B_perturbation_scale * randn(n, n);

% B(t) function handle
B_func = @(t) B0 + B1*cos(2*pi*t/T) + B2*sin(2*pi*t/T);

%% Generate K(t) matrix
% K(t) = K0 + K1*cos(2πt/T) + K2*sin(2πt/T) + K3*cos(4πt/T)

% Base matrix K0 - ensure it has full column rank
K0 = generate_full_rank_matrix(n, m);

% Perturbation matrices
K_perturbation_scale = 0.2;
K1 = K_perturbation_scale * randn(n, m);
K2 = K_perturbation_scale * randn(n, m);
K3 = K_perturbation_scale * 0.5 * randn(n, m);

% K(t) function handle
K_func = @(t) K0 + K1*cos(2*pi*t/T) + K2*sin(2*pi*t/T) + K3*cos(4*pi*t/T);

%% Verify periodicity
fprintf('Verifying system properties...\n');

% Check periodicity
A_error = norm(A_func(0) - A_func(T), 'fro');
B_error = norm(B_func(0) - B_func(T), 'fro');
K_error = norm(K_func(0) - K_func(T), 'fro');

if max([A_error, B_error, K_error]) < 1e-12
    fprintf('✓ System is periodic (max error: %.2e)\n', max([A_error, B_error, K_error]));
else
    fprintf('⚠ Periodicity check failed (max error: %.2e)\n', max([A_error, B_error, K_error]));
end

% Check stability of A(t) on average
t_sample = linspace(0, T, 100);
max_real_part = -inf;
for i = 1:length(t_sample)
    eigs_A = eig(A_func(t_sample(i)));
    max_real_part = max(max_real_part, max(real(eigs_A)));
end

if max_real_part < 0
    fprintf('✓ A(t) is stable (max Re(λ) = %.3f)\n', max_real_part);
elseif max_real_part < 0.5
    fprintf('○ A(t) is marginally stable (max Re(λ) = %.3f)\n', max_real_part);
else
    fprintf('⚠ A(t) may be unstable (max Re(λ) = %.3f)\n', max_real_part);
end

% Check controllability at t=0
K0_test = K_func(0);
rank_K0 = rank(K0_test);
if rank_K0 == m
    fprintf('✓ K(0) has full column rank (%d)\n', rank_K0);
else
    fprintf('⚠ K(0) is rank deficient (rank = %d, expected %d)\n', rank_K0, m);
end

fprintf('Random periodic system generated successfully.\n');

end

function A_stable = generate_stable_matrix(n)
% Generate a stable matrix (all eigenvalues have negative real parts)

% Method: Generate random matrix and ensure stability
max_attempts = 10;
attempt = 1;

while attempt <= max_attempts
    % Generate random matrix
    A_stable = randn(n, n);
    
    % Make it more likely to be stable
    A_stable = A_stable - (max(real(eig(A_stable))) + 0.5) * eye(n);
    
    % Check stability
    if max(real(eig(A_stable))) < -0.1
        return;
    end
    
    attempt = attempt + 1;
end

% Fallback: construct explicitly stable matrix
fprintf('  Using fallback stable matrix construction...\n');
A_stable = -eye(n) + 0.1 * randn(n, n);

end

function K_full = generate_full_rank_matrix(n, m)
% Generate an n x m matrix with full column rank

if m > n
    error('Cannot generate full rank matrix with m > n');
end

% Generate using QR decomposition for guaranteed full rank
[Q, ~] = qr(randn(n, n));
R_small = triu(randn(m, m));

% Ensure R has non-zero diagonal (full rank)
for i = 1:m
    if abs(R_small(i, i)) < 0.1
        R_small(i, i) = sign(R_small(i, i)) * (0.5 + rand());
    end
end

% Construct full rank matrix
K_full = Q(:, 1:m) * R_small;

% Verify full rank
if rank(K_full) < m
    fprintf('⚠ Fallback: Using orthogonal columns for K matrix\n');
    [Q, ~] = qr(randn(n, m), 0);
    K_full = Q;
end

end
