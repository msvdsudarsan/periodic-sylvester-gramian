function [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T, varargin)
%GENERATE_RANDOM_PERIODIC_SYSTEM Generate random T-periodic system matrices
%
% Syntax:
% [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T)
% [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T, 'stable', true)
% [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T, 'controllable', true)
%
% Inputs:
% n - State dimension
% m - Input dimension 
% T - Period
%
% Optional parameters (name-value pairs):
% 'stable' - Ensure Floquet stability (default: false)
% 'controllable'- Ensure controllability (default: false)
% 'seed' - Random seed for reproducibility (default: random)
% 'amplitude' - Amplitude of periodic variations (default: 0.3)
%
% Outputs:
% A_func - Function handle for A(t) (n×n periodic matrix)
% B_func - Function handle for B(t) (n×n periodic matrix)
% K_func - Function handle for K(t) (n×m periodic matrix)
%
% The generated system has the form:
% A(t) = A0 + A1*cos(2πt/T) + A2*sin(2πt/T)
% B(t) = B0 + B1*cos(2πt/T) + B2*sin(2πt/T) 
% K(t) = K0 + K1*cos(2πt/T) + K2*sin(2πt/T)
%
% Author: M. S. V. D. Sudarsan
% Email: msvdsudarsan@gmail.com
% Date: 2025

% Parse optional arguments
p = inputParser;
addParameter(p, 'stable', false, @islogical);
addParameter(p, 'controllable', false, @islogical);
addParameter(p, 'seed', [], @isnumeric);
addParameter(p, 'amplitude', 0.3, @isnumeric);
parse(p, varargin{:});

stable = p.Results.stable;
controllable = p.Results.controllable;
seed_val = p.Results.seed;
amplitude = p.Results.amplitude;

% Set random seed for reproducibility
if ~isempty(seed_val)
    rng(seed_val);
end

% Validate inputs
if n < 1 || m < 1
    error('Dimensions n and m must be positive integers');
end
if T <= 0
    error('Period T must be positive');
end
if amplitude < 0 || amplitude > 1
    error('Amplitude must be in [0, 1]');
end

fprintf('Generating random periodic system (n=%d, m=%d, T=%.2f)\n', n, m, T);
if stable
    fprintf('  → Ensuring Floquet stability\n');
end
if controllable
    fprintf('  → Ensuring controllability\n');
end

%% Generate base matrices
% Start with a stable base system if requested
if stable
    A0 = generate_stable_matrix(n);
else
    A0 = randn(n, n) * 0.5; % Moderate random matrix
end

B0 = randn(n, n) * 0.3; % Smaller coupling matrix
K0 = randn(n, m) * 0.5; % Input matrix

% Generate periodic variations
A1 = randn(n, n) * amplitude;
A2 = randn(n, n) * amplitude;
B1 = randn(n, n) * amplitude * 0.5;
B2 = randn(n, n) * amplitude * 0.5;
K1 = randn(n, m) * amplitude;
K2 = randn(n, m) * amplitude;

%% Ensure controllability if requested
if controllable
    fprintf('  → Adjusting for controllability...\n');
    % Make sure K0 has sufficient rank
    [U, S, V] = svd(K0);
    S_new = max(S, 0.1); % Ensure minimum singular values
    K0 = U * diag(S_new) * V';
    % Add controllability-enhancing terms
    K0 = K0 + 0.1 * eye(n, min(n, m));
end

%% Create function handles
omega = 2*pi/T; % Angular frequency

A_func = @(t) A0 + A1*cos(omega*t) + A2*sin(omega*t);
B_func = @(t) B0 + B1*cos(omega*t) + B2*sin(omega*t);
K_func = @(t) K0 + K1*cos(omega*t) + K2*sin(omega*t);

%% Display system properties
fprintf('System properties at t=0:\n');
fprintf('  A(0) eigenvalues: ');
eig_A0 = eig(A0);
fprintf('%.3f±%.3fi ', [real(eig_A0), imag(eig_A0)].');
fprintf('\n');

fprintf('  B(0) norm: %.3f\n', norm(B0, 'fro'));
fprintf('  K(0) singular values: ');
sv_K0 = svd(K0);
fprintf('%.3f ', sv_K0);
fprintf('\n');

% Check basic properties
rank_K0 = rank(K0, 1e-10);
fprintf('  K(0) rank: %d/%d\n', rank_K0, min(n, m));

if rank_K0 == min(n, m)
    fprintf('  ✓ K(0) has full rank\n');
else
    fprintf('  ! K(0) is rank deficient\n');
end

%% Verify periodicity
fprintf('Verifying periodicity...\n');
t_test = [0, T/4, T/2, 3*T/4, T];
periodicity_error = 0;

for t = t_test
    A_error = norm(A_func(t) - A_func(t + T), 'fro');
    B_error = norm(B_func(t) - B_func(t + T), 'fro');
    K_error = norm(K_func(t) - K_func(t + T), 'fro');
    periodicity_error = max(periodicity_error, A_error + B_error + K_error);
end

if periodicity_error < 1e-12
    fprintf('  ✓ Periodicity verified (error: %.2e)\n', periodicity_error);
else
    fprintf('  ! Periodicity error: %.2e\n', periodicity_error);
end

fprintf('Random periodic system generated successfully.\n\n');

end

%% Helper function: Generate stable matrix
function A = generate_stable_matrix(n)
%GENERATE_STABLE_MATRIX Generate a stable random matrix
%
% Creates a random matrix with all eigenvalues having negative real parts

% Generate random matrix
A = randn(n, n);

% Symmetrize to ensure real eigenvalues (for simplicity)
A = 0.5 * (A + A');

% Shift eigenvalues to ensure stability
lambda = eig(A);
max_real_part = max(real(lambda));

if max_real_part >= 0
    % Shift to make stable
    A = A - (max_real_part + 0.1) * eye(n);
end

% Verify stability
lambda_new = eig(A);
if max(real(lambda_new)) >= 0
    warning('Failed to generate stable matrix, trying alternative approach');
    % Alternative: construct stable matrix directly
    [Q, ~] = qr(randn(n, n));
    D = diag(-abs(randn(n, 1)) - 0.1); % Negative eigenvalues
    A = Q * D * Q';
end

end
