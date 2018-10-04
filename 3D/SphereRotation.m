%% SphereRotation.m
%
% This function computes the rotation from the (s,t) angles on the unit
% sphere to the (theta,phi) angles on the unit sphere. See Appendix A of
% Carvalho, Khatri, and Kim (2018) for the mathematical details.
%
% Written by A. D. Kim on 4/9/2018

function [ theta, phi ] = SphereRotation( s, t, theta0, phi0 )

%% compute the auxilliary variables

xi   = cos( theta0 ) * cos( phi0 ) * sin( s ) .* cos( t ) ...
     - sin( phi0 ) * sin( s ) .* sin( t ) ...
     + sin( theta0 ) * cos( phi0 ) * cos( s );

eta  = cos( theta0 ) * sin( phi0 ) * sin( s ) .* cos( t ) ...
     + cos( phi0 ) * sin( s ) .* sin( t ) ...
     + sin( theta0 ) * sin( phi0 ) * cos( s );

zeta = - sin( theta0 ) * sin( s ) .* cos( t ) + cos( theta0 ) * cos( s );

%% compute the rotated angles

theta = atan2( sqrt( xi.^2 + eta.^2 ), zeta );
phi   = atan2( eta, xi );

return