%% ComputeDLPmu.m
%
% This code computes the coefficients of the spherical harmonics expansion
% for the density of the double-layer potential on a closed surface using
% the discretized Galerkin method. See Section 5 of Carvalho, Khatri, and
% Kim (2018) for an explanation of this method.
%
% Written by A. D. Kim on 4/9/2018. Modified most recently on 5/4/2018.

function c_coeffs = ComputeDLPmu( n_order )

%% QUADRATURE RULES

N = n_order;
M = 2 * N;

% Gauss-Legendre/trapezoid product rule for theta (z = cos(theta))

[ z, wz ] = GaussLegendre( N );

% transform z quadrature rule for s quadrature rule

s  = 0.5 * pi * ( z' + 1 );
ws = 0.5 * pi * wz';

% trapezoid rules for t

dt   =  pi / N;
t    = -pi : dt : pi - dt;

%% COMPOSITE ANGLE ARRAYS

PHI   = repmat( t', N, 1 );
THETA = repelem( acos( z' ), M );
W_GL  = repelem( wz', M ) * dt;

T     = repmat( t', N, 1 );
S     = repelem( s, M );

%% SPHERICAL HARMONICS MATRIX

[ Y_GL, ~ ] = ComputeSphericalHarmonics( n_order, THETA, PHI );

%% COMPUTE THE DIRICHLET BOUNDARY DATA

[ ~, ~, x, ~, ~ ] = ComputeSurface( 0, 0, THETA, PHI );
[ f, ~, ~, ~ ]    = ComputeHarmonicFunction( x(:,1), x(:,2), x(:,3) );

%% COMPUTE SPHERICAL HARMONICS EXPANSION COEFFICIENTS

f_coeffs = Y_GL' * diag( W_GL ) * f;

%% ALLOCATE MEMORY FOR MATRICES

Kmtx = zeros( N * M, n_order^2 );

%% SET UP THE WAITBAR

wb = waitbar( 0, 'Computing the density. Please wait... ' );

%% COMPUTE THE FIRST INNER PRODUCT: BIE OPERATOR APPLIED TO SPHERICAL HARMONICS

for j = 1 : N * M
    
    % update the waitbar for the impatient user
    
    waitbar( j/(N*M), wb );
    
    % set the values of theta0 and phi0
    
    theta0 = THETA(j);
    phi0   = PHI(j);
   
    % compute ystar
    
    [ ~, ~, ystar, ~, ~ ] = ComputeSurface( 0, 0, theta0, phi0 );

    % compute the surface
    
    [ varTHETA, varPHI, y, nu, J ] = ComputeSurface( theta0, phi0, S, T );
    
    % compute the kernel for the BIE operator
    
    yd = [ ystar(1) - y(:,1), ystar(2) - y(:,2), ystar(3) - y(:,3) ];
    
    ydist = sqrt( yd(:,1).^2 + yd(:,2).^2 + yd(:,3).^2 );
    
    kernel = 0.5 * J .* ( nu(:,1) .* yd(:,1) + nu(:,2) .* yd(:,2) ...
        + nu(:,3) .* yd(:,3) ) ./ ydist.^3;

    % compute the integral operation
    
    k = 0;

    for n = 0 : n_order - 1
        
        for m = -n : n
            
            k = k + 1;
            
            % compute the function to be integrated
            
            ktemp1 = ( Compute_Ynm( n, m, varTHETA, varPHI ) ...
                - Compute_Ynm( n, m, theta0, phi0 ) ) .* kernel .* sin(S);
            
            % compute the integral in t
            
            ktemp2 = sum( reshape( ktemp1, M, N ), 1 ) / M;
                        
            % compute the integral in s
            
            Kmtx(j,k) = ktemp2 * ws;
            
        end
        
    end  
    
end

% close the wait bar

close( wb );

%% COMPUTE THE SECOND INNER PRODUCT

K = Y_GL' * diag( W_GL ) * Kmtx;

%% SOLVE THE GALERKIN SYSTEM OF EQUATIONS

c_coeffs = ( -eye( n_order^2 ) + K ) \ f_coeffs;

return
