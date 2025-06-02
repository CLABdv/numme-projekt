% Jacobi matrix for Ic
function J = JacIc(atilde, x0tilde, y0tilde, alphas, eta, omega)
     J = zeros(length(alphas),3);
     for i = 1:length(alphas)
         J(i,1) = eta * cos(omega * (x0tilde .* cos(alphas(i)) + y0tilde .* sin(alphas(i)))); %dF/da
         J(i,2) = -atilde * eta .* sin(omega *(x0tilde .*cos(alphas(i)) + y0tilde .*sin(alphas(i))))*omega .*cos(alphas(i)); %dF/dx0
         J(i,3) = -atilde * eta .* sin(omega *(x0tilde .*cos(alphas(i)) + y0tilde .*sin(alphas(i))))*omega .*sin(alphas(i)); %dF/dy0
     end
end