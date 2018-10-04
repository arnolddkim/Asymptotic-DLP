%% ComputeHarmonicFunction.m
%
% Computes the harmonic function and its gradient for use in validating 
% the evaluation of double- and single-layer potentials.
%
% Written by A. D. Kim on 5/4/2018.

function [ u, gradu_x, gradu_y, gradu_z ] = ComputeHarmonicFunction( x, y, z )

x0 = 5;
y0 = 4;
z0 = 3;

% compute the harmonic function

ulen = sqrt( ( x - x0 ).^2 + ( y - y0 ).^2 + ( z - z0 ).^2 );
u    = 1 ./ ulen;

% compute the components of its gradient

gradu_x = - ( x - x0 ) ./ ulen.^3;
gradu_y = - ( y - y0 ) ./ ulen.^3;
gradu_z = - ( z - z0 ) ./ ulen.^3;

return