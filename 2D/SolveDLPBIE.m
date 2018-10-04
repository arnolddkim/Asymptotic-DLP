function mu = SolveDLPBIE( f, ybdy, normal, Jacobian, kappa )
%% function mu = SolveDLPBIE( f, ybdy, normal, Jacobian, kappa )
%
% Solve the boundary integral equation for the density of the double-layer
% potential using a Nystrom method (periodic trapezoid rule).
%

% determine the length of f

N = length( f );

% compute dt

dt = 2 * pi / N;

% allocate memory for the kernel

Kbdy = zeros( N );

% compute the kernel

for i = 1 : N
    
    Kbdy(i,i) = -0.25 * kappa(i) / pi;
    
    for j = 1 : N
        
        if i ~= j

            ydiff     = [ ybdy(i,1) - ybdy(j,1) ybdy(i,2) - ybdy(j,2) ];
            distance  = sqrt( ydiff(:,1).^2 + ydiff(:,2).^2 );
            costheta  = sum( normal(j,:) .* ydiff, 2 ) ./ distance;                           
            Kbdy(i,j) = 0.5 / pi * costheta / distance;
            
        end
        
    end
    
end

% solve the linear system for the density

mu = ( Kbdy * diag( Jacobian ) * dt - 0.5 * eye(N) ) \ f;

return;