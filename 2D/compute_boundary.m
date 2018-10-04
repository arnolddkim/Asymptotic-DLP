%
% -----------------------
% MATLAB: compute_boundary.m
% -----------------------
%

function [ r, normal, Jacobian, kappa, sigma ] = compute_boundary( t, f )

%compute boundary and derivatives of boundary 

r(:,1)   = f{1,1}(t);
r(:,2)   = f{1,2}(t);
rp1(:,1) = f{2,1}(t); 
rp1(:,2) = f{2,2}(t); 
rp2(:,1) = f{3,1}(t);
rp2(:,2) = f{3,2}(t); 

% compute the Jacobian

Jacobian = sqrt( rp1(:,1).^2 + rp1(:,2).^2 );
     
% compute the curvature

kappa = ( rp1(:,1) .* rp2(:,2) - rp1(:,2) .* rp2(:,1) ) ...
    ./ Jacobian.^3;

% compute sigma = dot( y', y'' )

sigma = rp1(:,1) .* rp2(:,1) + rp1(:,2) .* rp2(:,2);

% compute the unit normal vector

normal = [ rp1(:,2) ./ Jacobian -rp1(:,1) ./ Jacobian ];   

return;