% sylvester_rhs.m
% Placeholder for a separate RHS file if you want to use it outside compute_periodic_gramian.

function dz_dt = sylvester_rhs(t, Z, A_func, B_func)
    % If Z is matrix form:
    dZ_dt = A_func(t) * Z + Z * B_func(t);
    dz_dt = dZ_dt;
end
