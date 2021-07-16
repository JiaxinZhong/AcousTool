function jd = SphBessel1D_Exact(max_order, z, varargin)

    z_row = z(:).';

    % Pre-calculated closed forms
    jd = cos(z_row)./z_row - sin(z_row)./z_row.^2;
    jd = [jd; (2*z_row.*cos(z_row)+(z_row.^2-2).*sin(z_row))./z_row.^3];
    jd = [jd; (-z_row.*(z_row.^2-9).*cos(z_row) + (4*z_row.^2-9).*sin(z_row))./z_row.^4];

    max_order_allowed = size(jd,1) - 1;
    if (max_order > max_order_allowed)
        error("The maximum order should be no larger than %d\n", max_order_allowed);
    end

    n = (0:max_order).';
    jd = reshape(jd(n+1,:), size(n.*z));
end
