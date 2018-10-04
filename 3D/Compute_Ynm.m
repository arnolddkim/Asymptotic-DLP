%% Compute Y_nm.m
%
% Evaluate the spherical harmonic function of degree n and order m 
% at polar angle theta and azimuthal angle phi. This code is a modification
% of the one written by Bing Jian that is available at 
%
% https://www.mathworks.com/matlabcentral/fileexchange/15377-real-valued-spherical-harmonics
%
% In particular, the modifications are contained in lines 28 through 30 and
% take into account the case when the order is negative.
%
% Modified by A. D. Kim on 5/4/2018

function Ynm = Compute_Ynm( N, M, THETA, PHI )

MM = abs(M);

Plm = legendre(N,cos(THETA));

if N~=0
  Plm = squeeze(Plm(MM+1,:,:));
end

a1 = ((2*N+1)/(4*pi));
a2 = factorial(N-M)/factorial(N+M);
C = sqrt(a1*a2);

if M < 0
    Plm = (-1)^MM * factorial(N-MM) / factorial(N+MM) * Plm;
end

Ynm = C * Plm.' .* exp( 1i * M * PHI );