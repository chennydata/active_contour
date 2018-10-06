function value = DiracDelta(x,epsilon)
% DIRAC Dirac function of x
%    DIRAC( x, epsilon ) Computes the derivative of the heaviside
%    function of x with respect to x. Regularized based on epsilon.

value = ( 1 ./ pi ) .* ( epsilon ./ ( epsilon.^2 + x.^2 ) );

