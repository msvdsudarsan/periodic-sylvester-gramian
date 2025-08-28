function [A_func, B_func, K_func] = generate_random_periodic_system(n_A, n_B, m, T, varargin)
%GENERATE_RANDOM_PERIODIC_SYSTEM Generate random periodic system matrices
%
% [A_func, B_func, K_func] = generate_random_periodic_system(n_A, n_B, m, T)
%
% Generates random periodic matrices for testing the Gramian computation:
%   A(t): n_A × n_A periodic matrix
%   B(t): n_B × n_B periodic matrix  
%   K(t): n_A × m periodic matrix (assuming n_A = n_B for Sylvester systems)
%
% Inputs:
%   n_A    - dimension of A matrix
%   n_B    - dimension of B matrix (typically n_A = n_B for Sylvester systems)
%   m      - number of control inputs
%   T      - period
%
% Optional parameters (name-value pairs):
%   'Seed'          - random seed for reproducibility (default: random)
%   'MaxEigenvalue' - maximum eigenvalue magnitude for stability (default: 2)
%   'Frequencies'   - fundamental frequencies to include (default: [1, 2, 3])
%   'Amplitudes'    - relative amplitudes for periodic terms (default: [0.3, 0.1, 0.05])
%   'Controllable'  - ensure system is controllable (default: true)
%
% Outputs:
%   A_func - function handle A(t) returning n_A × n_A matrix
%   B_func - function handle B(t) returning n_B × n_B matrix
%   K_func - function handle K(t) returning n_A × m matrix

%% Input validation and default parameters
if nargin < 4
    error('At least 4 input arguments required: n_A, n_B, m, T');
end

if n_A ~= n_B
    warning('For Sylvester systems, typically n_A = n_B. Proceeding with n_A = %d, n_B = %d', n_A, n_B);
end

% Parse optional arguments
p = inputParser;
addParameter(p, 'Seed', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
addParameter(p, 'MaxEigenvalue', 2, @(x) isscalar(x) && x > 0);
addParameter(p, 'Frequencies', [1, 2, 3], @(x) isvector(x) && all(x > 0));
addParameter(p, 'Amplitudes', [0.3, 0.1, 0.05], @(x) isvector(x) && all(x >= 0));
addParameter(p, 'Controllable', true, @islogical);
parse(p, varargin{:});

seed = p.Results.Seed;
max_eigenvalue = p.Results.MaxEigenvalue;
frequencies = p.Results.Frequencies;
amplitudes = p.Results.Amplitudes;
ensure_controllable = p.Results.Controllable;

% Set random seed if provided
if ~isempty(seed)
    rng(seed);
    fprintf('Random seed set to %d\n', seed);
end

% Ensure frequencies and amplitudes have same length
if length(amplitudes) < length(frequencies)
    amplitudes = [amplitudes, 0.05 * ones(1, length(frequencies) - length(amplitudes))];
elseif length(amplitudes) > length(frequencies)
    amplitudes = amplitudes(1:length(frequencies));
end

fprintf('Generating random periodic system:\n');
fprintf('• A(t): %d×%d, B(t): %d×%d, K(t): %d×%d\n', n_A, n_A, n_B, n_B, n_A, m);
fprintf('• Period T = %.3f\n', T);
fprintf('• Frequencies: %s\n', mat2str(frequencies));
fprintf('• Amplitudes: %s\n', mat2str(amplitudes, 3));

%% Generate base matrices
% A matrix: start with stable base matrix
A_base = generate_stable_matrix(n_A, max_eigenvalue);

% B matrix: start with stable base matrix  
B_base = generate_stable_matrix(n_B, max_eigenvalue);

% K matrix: random base with controlled norm
K_base = randn(n_A, m);
K_base = K_base / norm(K_base, 'fro') * sqrt(n_A * m); % Normalize

%% Generate periodic components
% For A(t)
A_periodic_terms = cell(length(frequencies), 1);
for i = 1:length(frequencies)
    % Random symmetric matrix for each frequency
    A_rand = randn(n_A, n_A);
    A_sym = (A_rand + A_rand') / 2; % Make symmetric
    A_sym = A_sym / norm(A_sym, 'fro') * amplitudes(i) * max_eigenvalue;
    A_periodic_terms{i} = A_sym;
end

% For B(t)
B_periodic_terms = cell(length(frequencies), 1);
for i = 1:length(frequencies)
    B_rand = randn(n_B, n_B);
    B_sym = (B_rand + B_rand') / 2;
    B_sym = B_sym / norm(B_sym, 'fro') * amplitudes(i) * max_eigenvalue;
    B_periodic_terms{i} = B_sym;
end

% For K(t)
K_periodic_terms = cell(length(frequencies), 1);
for i = 1:length(frequencies)
    K_rand = randn(n_A, m);
    K_rand = K_rand / norm(K_rand, 'fro') * amplitudes(i) * sqrt(n_A * m);
    K_periodic_terms{i} = K_rand;
end

%% Construct function handles
% Fundamental frequency
omega = 2*pi/T;

% A(t) function
A_func = @(t) compute_A(t, A_base, A_periodic_terms, frequencies, omega);

% B(t) function  
B_func = @(t) compute_B(t, B_base, B_periodic_terms, frequencies, omega);

% K(t) function
K_func = @(t) compute_K(t, K_base, K_periodic_terms, frequencies, omega, ensure_controllable);

%% Verification
fprintf('System generation completed.\n');

% Test function handles
try
    A_test = A_func(0);
    B_test = B_func(0);
    K_test = K_func(0);
    
    fprintf('✓ Function handles created successfully\n');
    fprintf('• A(0): %d×%d, eigenvalues in [%.3f, %.3f]\n', ...
            size(A_test,1), size(A_test,2), min(real(eig(A_test))), max(real(eig(A_test))));
    fprintf('• B(0): %d×%d, eigenvalues in [%.3f, %.3f]\n', ...
            size(B_test,1), size(B_test,2), min(real(eig(B_test))), max(real(eig(B_test))));
    fprintf('• K(0): %d×%d, norm = %.3f\n', ...
            size(K_test,1), size(K_test,2), norm(K_test,'fro'));
            
catch ME
    error('Error in generated function handles: %s', ME.message);
end

% Test periodicity
t1 = rand() * T;
t2 = t1 + T;
A_error = norm(A_func(t1) - A_func(t2), 'fro');
B_error = norm(B_func(t1) - B_func(t2), 'fro');
K_error = norm(K_func(t1) - K_func(t2), 'fro');

if max([A_error, B_error, K_error]) < 1e-12
    fprintf('✓ Periodicity verified (error < 1e-12)\n');
else
    warning('Periodicity check failed: max error = %.2e', max([A_error, B_error, K_error]));
end

end

%% Helper functions
function A = generate_stable_matrix(n, max_eigenvalue)
%GENERATE_STABLE_MATRIX Generate a stable matrix with controlled eigenvalues

if n == 1
    A = -0.1 - rand() * 0.5; % Stable scalar
    return;
end

% Generate random eigenvalues with negative real parts for stability
real_parts = -rand(n, 1) * max_eigenvalue; % Negative real parts
imag_parts = (rand(n, 1) - 0.5) * max_eigenvalue; % Small imaginary parts

% Ensure complex eigenvalues come in conjugate pairs
if mod(n, 2) == 1
    imag_parts(end) = 0; % Last eigenvalue is real
end

for i = 2:2:n-1
    imag_parts(i) = -imag_parts(i-1); % Conjugate pair
end

eigenvalues = real_parts + 1i * imag_parts;

% Generate random orthogonal matrix
[Q, ~] = qr(randn(n, n));

% Construct matrix A = Q * diag(eigenvalues) * Q'
A = Q * diag(eigenvalues) * Q';
A = real(A); % Should be real for real eigenvalue problems
end

function A_t = compute_A(t, A_base, periodic_terms, frequencies, omega)
%COMPUTE_A Compute A(t) with periodic components
A_t = A_base;

for i = 1:length(frequencies)
    freq = frequencies(i);
    A_t = A_t + periodic_terms{i} * cos(freq * omega * t);
end
end

function B_t = compute_B(t, B_base, periodic_terms, frequencies, omega)
%COMPUTE_B Compute B(t) with periodic components
B_t = B_base;

for i = 1:length(frequencies)
    freq = frequencies(i);
    B_t = B_t + periodic_terms{i} * sin(freq * omega * t);
end
end

function K_t = compute_K(t, K_base, periodic_terms, frequencies, omega, ensure_controllable)
%COMPUTE_K Compute K(t) with periodic components
K_t = K_base;

for i = 1:length(frequencies)
    freq = frequencies(i);
    % Mix of sin and cos for K(t)
    K_t = K_t + periodic_terms{i} * (cos(freq * omega * t) + 0.5 * sin(freq * omega * t));
end

% Ensure K(t) has sufficient magnitude for controllability
if ensure_controllable
    min_norm = 0.1;
    current_norm = norm(K_t, 'fro');
    if current_norm < min_norm
        K_t = K_t * (min_norm / current_norm);
    end
end
end
