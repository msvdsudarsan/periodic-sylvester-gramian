function [A_func, B_func, K_func] = generate_random_periodic_system(n, m, T)
% Generate random smooth periodic system matrices
%
% Inputs:
%   n - State dimension
%   m - Input dimension  
%   T - Period
%
% Outputs:
%   A_func, B_func, K_func - Function handles for periodic matrices

    % Base matrices
    A0 = randn(n, n); A0 = A0 - A0'; % Make skew-symmetric for stability
    A1 = 0.1 * randn(n, n);
    
    B0 = 0.1 * randn(n, n);
    B1 = 0.1 * randn(n, n);
    
    K0 = randn(n, m);
    K1 = 0.2 * randn(n, m);
    
    % Create smooth periodic functions
    A_func = @(t) A0 + A1 * cos(2*pi*t/T);
    B_func = @(t) B0 + B1 * sin(2*pi*t/T);
    K_func = @(t) K0 + K1 * sin(4*pi*t/T);
end
