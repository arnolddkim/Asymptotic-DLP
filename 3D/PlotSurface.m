%% PlotSurface.m
%
% Code to plot the surface used for the computation of the double-layer
% potential.
%

clear;

%% set the figure parameters

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% set the order of the trapezoid rules

N = 48;
M = 2*N;

%% compute the Gauss-Legendre quadrature rule in s

[ s, ws ] = GaussLegendre( N );
s  = 0.5 * pi * ( s + 1 );
ws = 0.5 * pi * ws;

%% compute the periodic trapezoid rule in t

dt = 2 * pi / M;
t  = -pi : dt : pi - dt;

%% compute mesh grid

[ S, T ] = meshgrid( s, t );

% stretch S and T into long column vectors

Svec = S(:);
Tvec = T(:);

%% allocate memory for surface area

SurfaceArea = zeros( 1, M * N );

%% set the value of theta0 and phi0

theta0 = 0;
phi0   = 0;

%% compute the surface and its properties

[ ~, ~, y, ~, ~ ] = ComputeSurface( theta0, phi0, Svec, Tvec );

%% reshape the y vector and impose periodicity in t

y1plot = reshape( y(:,1), M, N ); y1plot = [ y1plot; y1plot(1,:) ];
y2plot = reshape( y(:,2), M, N ); y2plot = [ y2plot; y2plot(1,:) ];
y3plot = reshape( y(:,3), M, N ); y3plot = [ y3plot; y3plot(1,:) ];

figure(1)
surf( y1plot, y2plot, y3plot ); 
axis equal;
xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
ylabel( '$x_{2}$', 'Interpreter', 'LaTeX' );
zlabel( '$x_{3}$', 'Interpreter', 'LaTeX' );


%% compute the slice through the xy-plane

sxy = pi/2;
txy = t;

[ Sxy, Txy ] = meshgrid( sxy, txy );

Sxy = Sxy(:);
Txy = Txy(:);

[ ~, ~, y_xy, ~, ~ ] = ComputeSurface( theta0, phi0, Sxy, Txy );

y_xy = [ y_xy; y_xy(1,:) ];

figure(2)
plot( y_xy(:,1), y_xy(:,2) );
axis equal;
xlim([min(y_xy(:,1))-0.5 max(y_xy(:,1))+0.5]);
ylim([min(y_xy(:,2)) max(y_xy(:,2))]);
xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
ylabel( '$x_{2}$', 'Interpreter', 'LaTeX' );

%% compute the slice through the xz-plane

sxz = s;
txz = [ -pi 0 ];

[ Sxz, Txz ] = meshgrid( sxz, txz );

Sxz = Sxz(:);
Txz = Txz(:);

[ ~, ~, y_xz, ~, ~ ] = ComputeSurface( theta0, phi0, Sxz, Txz );

figure(3)
plot( [ y_xz(1:2:end,1); y_xz(end:-2:2,1); y_xz(1,1) ], ...
      [ y_xz(1:2:end,3); y_xz(end:-2:2,3); y_xz(1,3) ] );
axis equal;
xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
ylabel( '$x_{3}$', 'Interpreter', 'LaTeX' );

