%% DLPAsymptApprox2D.m
%
% This code computes the asymptotic approximation of the close evaluation 
% of the double-layer potential in two dimensions using the method detailed
% in the paper entitled, "Asymptotic approximation for the close evaluation
% of double-layer potentials," by C. Carvalho, S. Khatri, and A. D. Kim
% (2018). This code calls the following codes.
%
% 1. boundary.m -- sets the boundary of the domain
%
% 2. compute_boundary.m -- computes the boundary and its properties on grid
%    points
%
% 3. SolveDLPBIE.m -- computes the numerical solution of the boundary
%    integral equation using the Nystrom method with the periodic trapezoid
%    rule (PTR).

clear;

%% set the figure parameters

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% compute the grid for the parametric boundary curve

N  = 128;
dt = 2 * pi / N;
t  = ( -pi : dt : pi - dt )';

%% compute the FFT grid

omega = fftshift( -N / 2 : N / 2 - 1 )';

%% compute the parametric boundary curve

boundarytype = 'kite'; % choices include circle, ellipse, star, and kite
fbdy = boundary( boundarytype ); 
[ ybdy, normal, Jacobian, kappa, sigma ] = compute_boundary( t, fbdy );

%% set the Dirichlet boundary data

x0   = 1.85;
y0   = 1.65;
dist = sqrt( ( ybdy(:,1) - x0 ).^2 + ( ybdy(:,2) - y0 ).^2 );
f    = - 0.5 / pi * log( dist );

%% compute the density using the Nystrom method

mu = SolveDLPBIE( f, ybdy, normal, Jacobian, kappa );

%% compute the first and second derivatives of the density

mup  = ifft( 1i * omega .* fft( mu ), 'symmetric' );
mupp = ifft(  -omega.^2 .* fft( mu ), 'symmetric' );

%% set the values of epsilon 

epsilon = [ 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 ];

%% allocate memory for quantities computed

xx     = zeros( N + 1, length(epsilon) );
yy     = zeros( N + 1, length(epsilon) );
error1 = zeros( N + 1, length(epsilon) );
error2 = zeros( N + 1, length(epsilon) );
error3 = zeros( N + 1, length(epsilon) );
error4 = zeros( N + 1, length(epsilon) );
Order1 = zeros( 1, N );
Order2 = zeros( 1, N );
Order3 = zeros( 1, N );

%% loop over all t

for j = 1 : N

    % compute the kernel for the asymptotic approximation
    
    yd         = ybdy(j,:) - ybdy;
    ydlength   = sqrt( yd(:,1).^2 + yd(:,2).^2 );
    ydhat      = [ yd(:,1) ./ ydlength yd(:,2) ./ ydlength ];
    ydhat(j,:) = 0;
    
    xi1 = normal(:,1) .* normal(j,1) + normal(:,2) .* normal(j,2);
    xi2 = normal(:,1) .* ydhat(:,1)  + normal(:,2) .* ydhat(:,2);
    xi3 = normal(j,1) .* ydhat(:,1)  + normal(j,2) .* ydhat(:,2);

    K1 = ( 2 * xi2 .* xi3 - xi1 ) ./ ydlength.^2;
    K2 = ( -2 * xi1 .* xi3 + xi2 .* ( 4 * xi3.^2 - 1 ) ) ./ ydlength.^3;
    
    % compute the integrals of the asymptotic approximation
    
    F1    = K1 .* ( mu - mu(j) ) .* Jacobian;
    F1(j) = - 0.5 * mupp(j) / Jacobian(j); % replace with asympt approx          
    U1    = sum( F1 ) / N;
    
    F2    = K2 .* ( mu - mu(j) ) .* Jacobian;
    F2(j) = - 0.25 * kappa(j) * mupp(j) / Jacobian(j); % replace with asympt approx
    U2    = sum( F2 ) / N;    
    
    % compute the local terms in the asymptotic approximation
    
    local = - 0.25 * sigma(j) * mup(j) / Jacobian(j)^4 ...
        + 0.25 * mupp(j) / Jacobian(j)^2;
    
    % compute the asymptotic approximation for each epsilon

    asympt1 = f(j) + epsilon * U1;
    asympt2 = asympt1 + epsilon.^2 * ( U2 + local );

    % compute the target points
    
    xx(j,:)  = ybdy(j,1) - epsilon * normal(j,1);
    yy(j,:)  = ybdy(j,2) - epsilon * normal(j,2);

    % compute the numerical approximation for each epsilon
    
    numerical1 = zeros(1,length(epsilon));
    numerical2 = zeros(1,length(epsilon));
    
    for k = 1 : length(epsilon)
        
        ydiff         = [ xx(j,k) - ybdy(:,1) yy(j,k) - ybdy(:,2) ];
        distance      = sqrt( ydiff(:,1).^2 + ydiff(:,2).^2 );
        costheta      = sum( normal .* ydiff, 2 ) ./ distance;         
        kernel        = 0.5 / pi * costheta ./ distance;              
        numerical1(k) = sum( kernel .* Jacobian .* mu ) * dt;
        numerical2(k) = - mu(j) ...
            + sum( kernel .* Jacobian .* ( mu - mu(j) ) ) * dt;
        
    end
    
    % compute the exact solution

    distance = sqrt( ( xx(j,:) - x0 ).^2 + ( yy(j,:) - y0 ).^2 );
    exact    = -0.5 / pi * log( distance ); 

    % compute the errors

    error1(j,:) = abs( exact - numerical1 );
    error2(j,:) = abs( exact - numerical2 );
    error3(j,:) = abs( exact - asympt1 );
    error4(j,:) = abs( exact - asympt2 );
    
    % estimate the order of accuracy for the asymptotic approximation

    k = 1 : length( epsilon ) - 1;
    P2 = polyfit( log( epsilon(k) ), log( error2(j,k) ), 1 );
    Order1(j) = P2(1);
    
    k  = find( error3(j,:) > 1e-15 );
    P3 = polyfit( log( epsilon(k) ), log( error3(j,k) ), 1 );
    Order2(j) = P3(1);
   
    if boundarytype == 'star'
        k  = find( error4(j,:) > 1e-11 );
    else
        k  = find( error4(j,:) > 1e-15 );
    end
    P4 = polyfit( log( epsilon(k) ), log( error4(j,k) ), 1 );
    Order3(j) = P4(1);

end

%% enforce periodicity

xx(N+1,:)     = xx(1,:);
yy(N+1,:)     = yy(1,:);
error1(N+1,:) = error1(1,:);
error2(N+1,:) = error2(1,:);
error3(N+1,:) = error3(1,:);
error4(N+1,:) = error4(1,:);

%% set the boundary points for plotting errors

i1 = 1 * N / 8 + 1;
i2 = 5 * N / 8 + 1;

%% plot the results

% error made by the PTR

figure(1)
pcolor( xx, yy, log10( error1 ) ); 
caxis( [ -15 1 ] );
shading interp;
colorbar;
xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
ylabel( '$x_{2}$', 'Interpreter', 'LaTeX' );
title( 'PTR Method', 'Interpreter', 'LaTeX' );
% label points y_A and y_B
hold on;

plot( [ ybdy(i1,1), ybdy(i2,1) ], [ ybdy(i1,2), ybdy(i2,2) ], 'rx' );

if boundarytype == 'star'

    text( ybdy(i1,1), ybdy(i1,2)-0.10, '$y_{A}$', ...
        'Interpreter', 'LaTeX', 'FontSize', 14 );

else
    
    text( ybdy(i1,1)+0.05, ybdy(i1,2), '$y_{A}$', ...
        'Interpreter', 'LaTeX', 'FontSize', 14 );

end

text( ybdy(i2,1)+0.05, ybdy(i2,2), '$y_{B}$', ...
    'Interpreter', 'LaTeX', 'FontSize', 14 );

hold off;

% error made by the subtraction method

figure(2)
pcolor( xx, yy, log10( error2 ) ); 
caxis( [ -15 1 ] );
shading interp;
colorbar;
xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
ylabel( '$x_{2}$', 'Interpreter', 'LaTeX' );
title( 'Subtraction Method', 'Interpreter', 'LaTeX' );
% label points y_A and y_B
hold on;

plot( [ ybdy(i1,1), ybdy(i2,1) ], [ ybdy(i1,2), ybdy(i2,2) ], 'rx' );

if boundarytype == 'star'

    text( ybdy(i1,1), ybdy(i1,2)-0.10, '$y_{A}$', ...
        'Interpreter', 'LaTeX', 'FontSize', 14 );

else
    
    text( ybdy(i1,1)+0.05, ybdy(i1,2), '$y_{A}$', ...
        'Interpreter', 'LaTeX', 'FontSize', 14 );

end

text( ybdy(i2,1)+0.05, ybdy(i2,2), '$y_{B}$', ...
    'Interpreter', 'LaTeX', 'FontSize', 14 );

hold off;

% error made by the O(epsilon^2) asymptotic approximation

figure(3)
pcolor( xx, yy, log10( error3 ) ); 
caxis( [ -15 1 ] );
shading interp;
colorbar;
xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
ylabel( '$x_{2}$', 'Interpreter', 'LaTeX' );
title( '$O(\epsilon^{2})$ Approximation', 'Interpreter', 'LaTeX' );
% label points y_A and y_B
hold on;

plot( [ ybdy(i1,1), ybdy(i2,1) ], [ ybdy(i1,2), ybdy(i2,2) ], 'rx' );

if boundarytype == 'star'

    text( ybdy(i1,1), ybdy(i1,2)-0.10, '$y_{A}$', ...
        'Interpreter', 'LaTeX', 'FontSize', 14 );

else
    
    text( ybdy(i1,1)+0.05, ybdy(i1,2), '$y_{A}$', ...
        'Interpreter', 'LaTeX', 'FontSize', 14 );

end

text( ybdy(i2,1)+0.05, ybdy(i2,2), '$y_{B}$', ...
    'Interpreter', 'LaTeX', 'FontSize', 14 );

hold off;

% error made by the O(epsilon^3) asymptotic approximation

figure(4)
pcolor( xx, yy, log10( error4 ) ); 
caxis( [ -15 1 ] );
shading interp;
colorbar;
xlabel( '$x_{1}$', 'Interpreter', 'LaTeX' );
ylabel( '$x_{2}$', 'Interpreter', 'LaTeX' );
title( '$O(\epsilon^{3})$ Approximation', 'Interpreter', 'LaTeX' );
% label points y_A and y_B
hold on;

plot( [ ybdy(i1,1), ybdy(i2,1) ], [ ybdy(i1,2), ybdy(i2,2) ], 'rx' );

if boundarytype == 'star'

    text( ybdy(i1,1), ybdy(i1,2)-0.10, '$y_{A}$', ...
        'Interpreter', 'LaTeX', 'FontSize', 14 );

else
    
    text( ybdy(i1,1)+0.05, ybdy(i1,2), '$y_{A}$', ...
        'Interpreter', 'LaTeX', 'FontSize', 14 );

end

text( ybdy(i2,1)+0.05, ybdy(i2,2), '$y_{B}$', ...
    'Interpreter', 'LaTeX', 'FontSize', 14 );

hold off;

% error at y_A - epsilon * nu_A for all four methods

figure(5)
loglog( epsilon, error1(i1,:), 'ko-', ...
    epsilon, error2(i1,:), 's-',     ...
    epsilon, error3(i1,:)+eps, 'x-', ...
    epsilon, error4(i1,:)+eps, '+-' );
grid;
axis( [ 1e-7 1 1e-20 1e5 ])
xlabel( '$\epsilon$', 'Interpreter', 'LaTeX' );
ylabel( 'error at $y_{A} - \epsilon \nu_{A}$', 'Interpreter', 'LaTeX' );
legend( {'PTR', 'Subtraction', '$O(\epsilon^{2})$', '$O(\epsilon^{3})$'}, 'Interpreter', 'LaTeX', 'Location', 'Southeast' );

% error at y_B - epsilon * nu_B for all four methods

figure(6)
loglog( epsilon, error1(i2,:), 'ko-', ...
    epsilon, error2(i2,:), 's-',     ...
    epsilon, error3(i2,:)+eps, 'x-', ...
    epsilon, error4(i2,:)+eps, '+-' );
grid;
axis( [ 1e-7 1 1e-20 1e5 ])
xlabel( '$\epsilon$', 'Interpreter', 'LaTeX' );
ylabel( 'error at $y_{B} - \epsilon \nu_{B}$', 'Interpreter', 'LaTeX' );
legend( {'PTR', 'Subtraction', '$O(\epsilon^{2})$', '$O(\epsilon^{3})$'}, 'Interpreter', 'LaTeX', 'Location', 'Southeast' );

% estimated order of accuracy for the subtraction method and the asymptotic
% approximations

figure(7)
plot( t, Order1, 'o-', t, Order2, 'x-', t, Order3, '+-' )
axis([ -4 4 0 4.25 ])
grid;
xlabel( '$t$', 'Interpreter', 'LaTeX' );
ylabel( 'estimated order', 'Interpreter', 'LaTeX' );
if boundarytype == 'star'
    legend( {'Subtraction', '$O(\epsilon^{2})$', ...
        '$O(\epsilon^{3})$'}, 'Interpreter', 'LaTeX', ...
        'Location', 'Northwest' );
else
    legend( {'Subtraction', '$O(\epsilon^{2})$', ...
        '$O(\epsilon^{3})$'}, 'Interpreter', 'LaTeX', ...
        'Location', 'Northeast' );
end