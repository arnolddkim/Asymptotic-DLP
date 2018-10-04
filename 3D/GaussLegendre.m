%% function [x,w] = GaussLegendre( N )
%
% Compute the Gauss Legendre quadrature points and weights using
% the method given in Spectral Methods in MATLAB by L. N. Trefethen (2000).
%
% Direct implementation by A. D. Kim.

function [x,w] = GaussLegendre( N )

  beta     = 0.5 ./ sqrt( 1.0 - ( 2.0 * (1:N-1) ).^(-2) );
  T        = diag( beta, 1 ) + diag( beta, -1 );
  [ V, D ] = eig( T );
  
  x     = diag( D );
  [x,i] = sort( x' );
  w     = 2 * V(1,i).^2;

return
