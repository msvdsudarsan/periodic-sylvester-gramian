function plot_convergence(Ns, sigmas)
figure; 
semilogy(Ns, abs(sigmas - sigmas(end) + eps), 'o-','LineWidth',1.2);
grid on; xlabel('Quadrature nodes N'); ylabel('|sigma_{min}^{(N)} - sigma_{min}^{(N_{max})}|');
title('Convergence of \sigma_{min}(W) with quadrature refinement');
end
