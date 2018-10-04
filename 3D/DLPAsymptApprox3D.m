%% DLPAsymptApprox3D.m
%
% This code computes the asymptotic approximation of the close evaluation 
% of the double-layer potential in three dimensions using the method 
% detailed in the paper entitled, "Asymptotic approximation for the close 
% evaluation of double-layer potentials," by C. Carvalho, S. Khatri, and 
% A. D. Kim (2018). This code calls the following codes.
%
% 1. ComputePlane.m -- sets the plane over which the double-layer potential
%    is evaluated.
%
% 2. GaussLegendre.m -- computes the quadrature points and weights for the 
%    Gauss-Legendre quadrature rule based on Trefethen's code from the book,
%    "Spectral Methods in Matlab."
%
% 3. ComputeDLPmu.m -- computes the solution of the boundary integral equation
%    using the Galerkin method. This is the main bottle-neck of this code, so
%    we have saved the results from our computation with N = 48 in the MAT 
%    file, "mushroom_48_density.mat".
%
% 4. ComputeSurface.m -- computes the mesh over the surface in the rotated
%    coordinate system in which (theta0,phi0) gives the coordinates of the
%    north pole.
%
% 5. ComputeSphericalHarmonics.m -- computes a matrix whose columns are the
%    evaluation of the spherical harmonics (in a particular order) on the
%    mesh points provided.
%
% 6. ComputeHarmonicFunction.m -- computes the harmonic function used for
%    the Dirichlet boundary data as well as the exact solution for computing
%    error.

clear;
close all;

%% set the figure parameters

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% set the plane over which to compute the error

plane = 'xz'; %xy, xz, or yz 
Nystar = 50; 

%% compute the grid points on that plane

[ t0, s0 ] = ComputePlane( plane, Nystar );

%% set the order of the Galerkin method and quadrature rule

N = 48;
M = 2*N;

%% compute the Gauss-Legendre quadrature rule in s

[ z, wz ] = GaussLegendre( N ); % Note that -1 < z < 1.
s  = 0.5 * pi * ( z + 1 );      % Map to 0 < s < pi
ws = 0.5 * pi * wz;             % Scale the weights by pi/2

%% compute the periodic trapezoid rule in t

dt = 2 * pi / M;
t  = -pi : dt : pi - dt;

%% compute mesh grid

[ S, T ] = meshgrid( s, t );

% stretch S and T into long column vectors

Svec = S(:);
Tvec = T(:);

%% compute the DLP density spherical harmonics expansion coefficients

% mu_nm = ComputeDLPmu( N ); % uncomment to compute the density
load mushroom_48_density;    % uncomment to use the precomputed density

disp( 'Computed density' );

%% set the values of epsilon 

epsilon = [ 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 5e-1 ];

%% allocate memory for quantities computed

xx     = zeros( length(t0), length(epsilon) );
yy     = zeros( length(t0), length(epsilon) );
zz     = zeros( length(t0), length(epsilon) );
error1 = zeros( length(t0), length(epsilon) );
error2 = zeros( length(t0), length(epsilon) );
order1 = zeros( 1, length(t0) );
order2 = zeros( 1, length(t0) );

%% loop over all evaluation points

for j = 1 : length(t0)

    % set the value of theta0 and phi0
    
    theta0 = s0(j);
    phi0   = t0(j);

    % compute the boundary point and normal at (theta0,phi0)
    
    [ ~, ~, ystar, nustar, ~ ] = ComputeSurface( 0, 0, theta0, phi0 );

    % compute the computational grid
    
    [ theta, phi, y, nu, J ] = ComputeSurface( theta0, phi0, Svec, Tvec );
    
    % compute the spherical harmonics matrix

    Ynm     = ComputeSphericalHarmonics( N, theta, phi );
    Ynmstar = ComputeSphericalHarmonics( N, theta0, phi0 );

    % compute the density

    mu     = Ynm * mu_nm;
    mustar = Ynmstar * mu_nm;

    % compute the kernel for the asymptotic approximation
   
    Yd    = [ ystar(1) - y(:,1), ystar(2) - y(:,2), ystar(3) - y(:,3) ];
    Ydlen = sqrt( Yd(:,1).^2 + Yd(:,2).^2 + Yd(:,3).^2 );
    
    xi1 = ( nu(:,1) .* Yd(:,1) + nu(:,2) .* Yd(:,2) ...
        + nu(:,3) .* Yd(:,3) ) ./ Ydlen;
    
    xi2 = ( nustar(1) .* Yd(:,1) + nustar(2) .* Yd(:,2) ...
        + nustar(3) .* Yd(:,3) ) ./ Ydlen;

    xi3 = ( nustar(1) .* nu(:,1) + nustar(2) .* nu(:,2) ...
        + nustar(3) .* nu(:,3) );
    
    L1kernel = ( 1.5 * xi1 .* xi2 - 0.5 * xi3 ) ...
        ./ Ydlen.^3 .* J .* sin( Svec );

    % compute the integral of the asymptotic approximation with respect to 
    % the azimuthal angle
    
    Fbar = sum( reshape( L1kernel .* ( mu - mustar ), M, N ), 1 ) / M;
    
    % compute the integral of the asymptotic approximation with respect to 
    % polar angle

    U = real( sum( Fbar .* ws ) );
    
    % compute the asymptotic approximation

    fstar  = ComputeHarmonicFunction( ystar(1), ystar(2), ystar(3) );
    asympt = fstar + epsilon * U;

    % compute the target points

    xx(j,:) = ystar(1) - epsilon * nustar(1);
    yy(j,:) = ystar(2) - epsilon * nustar(2);
    zz(j,:) = ystar(3) - epsilon * nustar(3);
    
    % compute the numerical approximation for each epsilon
    
    numerical = zeros( 1, length(epsilon) );
    
    for k = 1 : length( epsilon )
           
        % compute the Yd-vector
        
        Yd1 = xx(j,k) - y(:,1);
        Yd2 = yy(j,k) - y(:,2);
        Yd3 = zz(j,k) - y(:,3);
        
        % compute the length of Yd
        
        Ydlength = sqrt( Yd1.^2 + Yd2.^2 + Yd3.^2 );
            
        % compute the kernel for the double-layer potential
    
        K = 0.5 * J .* sin(Svec) .* ...
            ( ( nu(:,1) .* Yd1  + nu(:,2) .* Yd2  + nu(:,3) .* Yd3 ) ...
            ./ Ydlength.^3 );
         
        % compute integral with respect to the azimuthal angle
    
        Fbar = sum( reshape( K .* ( mu - mustar ), M, N ), 1 ) / M;
        
        % compute the numerical approximation of the double-layer potential

        numerical(k) = real( - mustar + sum( Fbar .* ws ) );

    end
    
    %% compute the exact solution
    
    exact = ComputeHarmonicFunction( xx(j,:), yy(j,:), zz(j,:) );
        
    %% compute the errors

    error1(j,:) = abs( exact - numerical );
    error2(j,:) = abs( exact - asympt );

    %% estimate the order of accuracy

    k = 3:5;
    P1 = polyfit( log( epsilon(k) ), log( error1(j,k) ), 1 );
    order1(j) = P1(1);

    P2 = polyfit( log( epsilon ), log( error2(j,:) ), 1 );
    order2(j) = P2(1);
    
end

%% plot the results

if plane == 'xz'

    indx = 10;

    starget = s0(indx);
    ttarget = 0;

    [ ~, ~, ybdy, ~, ~ ] = ComputeSurface( 0, 0, starget, ttarget );

    figure(1)
    pcolor( xx, zz, log10(error1) );
    shading interp;
    caxis( [ -15 1 ] );
    colorbar;
    axis( [ -2.5 2.5 -2.5 2 ] )
    xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
    ylabel( '$x_{3}$', 'Interpreter', 'LaTeX' );
    title( 'Numerical Approximation', 'Interpreter', 'LaTeX' );
    hold on;
    plot( ybdy(1), ybdy(3), 'rx' );
    text( ybdy(1)+0.1, ybdy(3), '$y_{A}$', 'Interpreter', 'LaTeX', ...
        'FontSize', 14 );
    hold off;
  
    
    figure(2)
    pcolor( xx, zz, log10(error2) );
    shading interp;
    caxis( [ -15 1 ] );
    colorbar;
    axis( [ -2.5 2.5 -2.5 2 ] )
    xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
    ylabel( '$x_{3}$', 'Interpreter', 'LaTeX' );
    title( '$O(\epsilon^{2})$ Approximation', 'Interpreter', 'LaTeX' );
    hold on;
    plot( ybdy(1), ybdy(3), 'rx' );
    text( ybdy(1)+0.1, ybdy(3), '$y_{A}$', 'Interpreter', 'LaTeX', ...
        'FontSize', 14 );
    hold off;
   

    figure(3)
    loglog( epsilon, error1(indx,:), 'o-', epsilon, error2(indx,:), 'x-' );
    ylim( [ 1e-15 1 ]);
    xlabel( '$\epsilon$', 'Interpreter', 'LaTeX' );
    ylabel( 'error at $y_{A} - \epsilon \nu_{A}$', 'Interpreter', 'LaTeX' );
    grid
    legend( {'Numerical', '$O(\epsilon^{2})$'}, ...
        'Interpreter', 'LaTeX', 'Location', 'Northwest' );
   
    
    figure(4)
    plot( [ s0(1:(end-1)/2) 2*pi-s0((end-1)/2+1:end) ], order1,'o-', ...
          [ s0(1:(end-1)/2) 2*pi-s0((end-1)/2+1:end) ], order2 , 'x-');
    xlim([ 0 2*pi ]);
    xlabel( '$s_{0}$ (extended polar angle)', 'Interpreter', 'LaTeX' );
    ylabel( 'estimated order', 'Interpreter', 'LaTeX' )
    grid;
    legend( {'Numerical', '$O(\epsilon^{2})$'}, ...
        'Interpreter', 'LaTeX', 'Location', 'Southeast' );
 
    
elseif plane == 'xy' 
        
    indx = 30;

    starget = pi/2;
    ttarget = t0(indx);

    [ ~, ~, ybdy, ~, ~ ] = ComputeSurface( 0, 0, starget, ttarget );

    figure(1)
    pcolor( xx, yy, log10(error1) );
    shading interp;
    caxis( [ -15 1 ] );
    colorbar;
    axis( [ -2.5 2.5 -4 4] )
    xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
    ylabel( '$x_{2}$', 'Interpreter', 'LaTeX' );
    title( 'Numerical Approximation', 'Interpreter', 'LaTeX' );
    hold on;
    plot( ybdy(1), ybdy(2), 'rx' );
    text( ybdy(1)+0.01, ybdy(2), '$y_{B}$', 'Interpreter', 'LaTeX', ...
        'FontSize', 14 );
    hold off;


    figure(2)
    pcolor( xx, yy, log10(error2) );
    shading interp;
    caxis( [ -15 1 ] );
    colorbar;
    axis( [ -2.5 2.5 -4 4] )
    xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
    ylabel( '$x_{2}$', 'Interpreter', 'LaTeX' );
    title( '$O(\epsilon^{2})$ Approximation', 'Interpreter', 'LaTeX' );
    hold on;
    plot( ybdy(1), ybdy(2), 'rx' );
    text( ybdy(1)+0.01, ybdy(2), '$y_{B}$', 'Interpreter', 'LaTeX', ...
        'FontSize', 14 );
    hold off;


    figure(3)
    loglog( epsilon, error1(indx,:), 'o-', epsilon, error2(indx,:), 'x-' );
    ylim( [ 1e-15 1 ])
    xlabel( '$\epsilon$', 'Interpreter', 'LaTeX' );
    ylabel( 'error at $y_{B} - \epsilon \nu_{B}$', 'Interpreter', 'LaTeX' );
    grid
    legend( {'Numerical', '$O(\epsilon^{2})$'}, ...
        'Interpreter', 'LaTeX', 'Location', 'Northwest' );

    
    figure(4)
    plot( t0, order1, 'o-', t0, order2, 'x-' );
    xlabel( '$t_{0}$ (azimuthal angle)', 'Interpreter', 'LaTeX' );
    ylabel( 'estimated order', 'Interpreter', 'LaTeX' )
    grid;
    legend( {'Numerical', '$O(\epsilon^{2})$'}, ...
        'Interpreter', 'LaTeX', 'Location', 'Southwest' );


end