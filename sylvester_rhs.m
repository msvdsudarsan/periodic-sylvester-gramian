function dZ = sylvester_rhs(t, Z_vec, A_func, B_func)
% Vectorized RHS: dZ/dt = A(t) Z + Z B(t)
n = round(sqrt(numel(Z_vec)));
Z = reshape(Z_vec, n, n);
dZ_mat = A_func(t) * Z + Z * B_func(t);
dZ = dZ_mat(:);
end
