%% ComputePlane.m
%
% This function computes a set of grid points in the polar and azimuthal
% angles over a particular plane: xy, yz, xz, etc. over which the 
% double-layer potential is evaluated in DLPAsymptApprox3D.m.
%

function [ t0, s0] = ComputePlane(plane, N )

if strcmp(plane,'xy')

    N = N/2+1; 
    s0 = linspace(pi/2,pi/2,N);
    t0 = linspace(0,pi,N); 
    t0 = [fliplr(-t0(2:end)),t0];
    s0 = [fliplr(s0(2:end)),s0];

elseif strcmp(plane,'yz')
    N = N/2+1; 
    
    s0 = linspace(0,pi,N);
    t0 = 0*s0; 
    
    s0 = [s0,fliplr(s0(1:end-1))];
    t0 = [t0,t0(1:end-1)+pi];
    
elseif strcmp(plane,'xz')
    N = N/2+1; 
    
    s0 = linspace(0,pi,N);
    t0 = 0*s0; % - pi/2; 
    
    s0 = [s0,fliplr(s0(1:end-1))];
    t0 = [t0,t0(1:end-1)+pi];
    
elseif strcmp(plane,'diag1_yz')
    N = N/2+1; 
    
    s0 = linspace(0,pi,N);
    t0 = 0*s0 - pi/8; 
    
    s0 = [s0,fliplr(s0(1:end-1))];
    t0 = [t0,t0(1:end-1)+pi];
    
elseif strcmp(plane,'diag2_yz')
    N = N/2+1; 
    
    s0 = linspace(0,pi,N);
    t0 = 0*s0 - pi/4; 
    
    s0 = [s0,fliplr(s0(1:end-1))];
    t0 = [t0,t0(1:end-1)+pi];
    
end