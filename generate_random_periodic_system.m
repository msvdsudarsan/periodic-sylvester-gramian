function [A_func,B_func,K_func] = generate_random_periodic_system(n,m,seed)
% Periodic test system with smooth time dependence
if nargin>=3 && ~isempty(seed)
    rng(seed);
end
A0 = randn(n); A0 = 0.5*(A0 - A0');     % skew-symmetric base (stable rotation)
Da = diag(rand(n,1));
Db = diag(rand(n,1));
Kb = randn(n,m);
A_func = @(t) A0 + 0.1*(cos(t)*Da + sin(t)*eye(n));
B_func = @(t) 0.5*(sin(t)*Db + cos(t)*eye(n));
K_func = @(t) Kb + 0.2*cos(t)*Kb;
end
