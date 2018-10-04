%% ComputeSurface.m
%
% Compute the parametric surface mesh points, their corresponding outward 
% unit normals, and their Jacobians relative to the point at (theta0,phi0).
% The user must supply the radial function and its partial derivatives.
%
% Written by A. D. Kim on 4/20/2018

function [ theta, phi, Y, NU, J ] = ComputeSurface( theta0, phi0, s, t )

[ theta, phi ] = SphereRotation( s, t, theta0, phi0 );

%% compute the vector on the unit sphere and its partial derivatives

x = sin( theta ) .* cos( phi );
y = sin( theta ) .* sin( phi );
z = cos( theta );

x_theta =   cos( theta ) .* cos( phi );
y_theta =   cos( theta ) .* sin( phi );
z_theta = - sin( theta );

x_phi   = - sin( theta ) .* sin( phi );
y_phi   =   sin( theta ) .* cos( phi );
z_phi   = 0 * z;

%% compute the radius function and its partial derivatives

% %Sphere or Ellipsoid
% R   = 1 * ones(size(x));
% R_x = 0 * R;
% R_y = 0 * R;
% R_z = 0 * R;

% Mushroom
R   = 2 - 1 ./ ( 1 + 100 * ( -z + 1 ).^2 );
R_x = 0 * R;
R_y = 0 * R;
R_z = - 200 * ( 1 - z ) ./ ( 1 + 100 * ( 1 - z ).^2 ).^2;

% % Peanut 
% 
% R = sqrt( z.^2 - x.^2 - y.^2 + sqrt(1.1 - 4*(x.^2 + y.^2).*z.^2) ); 
% R_x = (-x - 2.*x.*z.^2./sqrt(1.1 - 4*(x.^2 + y.^2).*z.^2)) ./ R; 
% R_y = (-y - 2.*y.*z.^2./sqrt(1.1 - 4*(x.^2 + y.^2).*z.^2)) ./ R;
% R_z = (z - 2.*z.*(x.^2 + y.^2)./sqrt(1.1 - 4*(x.^2 + y.^2).*z.^2)) ./ R;

%% compute the surface vector and its partial derivatives

% %sphere
% A = 2;
% B = 2;
% C = 2;

% %ellipse
% A = 1.25;
% B = 2.50;
% C = 5.00;

% %ellipse2
% A = 1;
% B = 1.5;
% C = 2;

% % ellipse3
% A = 1;
% B = 2;
% C = 5;

% peanut and mushroom

A = 1; 
B = 2; 
C = 1; 

Y = [ A * R .* x, B * R .* y, C * R .* z ];

Y_x = [ A * R zeros(size(y)) zeros(size(z)) ] ...
    + [ A * R_x .* x B * R_x .* y C * R_x .* z ];

Y_y = [ zeros(size(x)) B * R zeros(size(z)) ] ...
    + [ A * R_y .* x B * R_y .* y C * R_y .* z ];

Y_z = [ zeros(size(x)) zeros(size(y)) C * R ] ...
    + [ A * R_z .* x B * R_z .* y C * R_z .* z ];
    
%% compute the partial derivatives of the surface

Y_theta = [x_theta .* Y_x(:,1) x_theta .* Y_x(:,2) x_theta .* Y_x(:,3)] ...
        + [y_theta .* Y_y(:,1) y_theta .* Y_y(:,2) y_theta .* Y_y(:,3)] ...
        + [z_theta .* Y_z(:,1) z_theta .* Y_z(:,2) z_theta .* Y_z(:,3)]; 
    
Y_phi = [x_phi .* Y_x(:,1) x_phi .* Y_x(:,2) x_phi .* Y_x(:,3)] ...
        + [y_phi .* Y_y(:,1) y_phi .* Y_y(:,2) y_phi .* Y_y(:,3)] ...
        + [z_phi .* Y_z(:,1) z_phi .* Y_z(:,2) z_phi .* Y_z(:,3)];     

%% compute the unit normal vectors

Jvec = cross( Y_theta, Y_phi, 2 ); 
Jlen = sqrt( Jvec(:,1).^2 + Jvec(:,2).^2 + Jvec(:,3).^2 );
NU   = [ Jvec(:,1)./Jlen, Jvec(:,2)./Jlen, Jvec(:,3)./Jlen ];

%% compute the Jacobian  

J = Jlen ./ sin(theta);

%% SPECIAL CASE: When theta == 0 we need to define NU and J differently

indx = find( theta == 0 );

if isempty( indx ) == 0

    % set the normal to be z-hat

    NU(indx,:) = [ 0 0 1 ];

    % set the correct limit for the Jacobian

    J(indx) = sqrt( ( - B * C * R(indx) .* R_x(indx) ).^2 ...
                  + ( - A * C * R(indx) .* R_y(indx) ).^2 ...
                  + (   A * B * R(indx) .* R(indx)   ).^2 );

end

return