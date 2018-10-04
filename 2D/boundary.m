function f = boundary( boundarytype )
%% function f = boundary( boundarytype )
%
% Compute the boundary curve and its derivatives
%

if strcmp( boundarytype, 'circle' )

    a  = 1;
    
    f  = { @(t)  a * cos(t), ...
           @(t)  a * sin(t); 
           @(t) -a * sin(t), ...
           @(t)  a * cos(t);
           @(t) -a * cos(t), ...
           @(t) -a * sin(t) };

elseif strcmp( boundarytype, 'ellipse' )

    a  = 2;
    b  = 1;           
    
    f  = { @(t)  a * cos(t), ...
           @(t)  b * sin(t); 
           @(t) -a * sin(t), ...
           @(t)  b * cos(t);
           @(t) -a * cos(t), ...
           @(t) -b * sin(t) };
    
elseif strcmp( boundarytype, 'star' ) 
    
    a = 0.3;
    w = 5;

    r   = @(t) 1 + a * cos( w * t );
    rp1 = @(t) - w * a * sin( w * t );
    rp2 = @(t) - w * w * a * cos( w * t );

    f  = { @(t) r(t) .* cos(t), ...
           @(t) r(t) .* sin(t);
           @(t) rp1(t) .* cos(t) - r(t) .* sin(t), ...
           @(t) rp1(t) .* sin(t) + r(t) .* cos(t);
           @(t)( rp2(t) - r(t) ) .* cos(t) - 2 * rp1(t) .* sin(t), ...
           @(t)( rp2(t) - r(t) ) .* sin(t) + 2 * rp1(t) .* cos(t) };
        
elseif strcmp( boundarytype, 'kite' )

    f = { @(t) cos(t) + 0.65 * cos(2*t) - 0.65, ...
          @(t) 1.5 * sin(t);
          @(t) - sin(t) - 1.3 * sin(2*t), ...
          @(t) 1.5 * cos(t);
          @(t) - cos(t) - 2.6 * cos(2*t), ...
          @(t) -1.5 * sin(t) };
    
end



