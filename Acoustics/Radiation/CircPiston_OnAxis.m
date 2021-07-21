% ==========================================================================
% The axial sound field radiated by a circular piston.
% Using the exact closed form solution
% --------------------------------------------------------------------------
% INPUT
%   radius  --- The radius of the circular piston
%   k       --- The wavenumber
%   z       --- The coordinate of the field point
% ==========================================================================
function [p, vel] = CircPiston_OnAxis(radius, k, z)

	ka = k*radius;
    za = z/radius;
    p = -2*1i*1.21*343 ...
        .* exp(1i * ka/2 .* (sqrt(1+za.^2)+za)) .*...
        sin(ka/2 .* (sqrt(1+za.^2)-za));

    vel.x = 0;
    vel.y = 0;
	vel.z = exp(1i*k*z) - za .* (1+za.^2).^(-1/2) ...
        .* exp(1i * ka .*sqrt(1+za.^2));
end
