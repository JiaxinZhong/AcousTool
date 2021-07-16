function hd = SphHankel1D_Exact(max_order, z)

    z_row = z(:).';

    % Pre-calculated closed forms
    hd = (1i+z_row)./z_row.^2 .* exp(1i*z_row);
    hd = [hd; (2*1i+2*z_row-1i*z_row.^2)./z_row.^3];
    hd = [hd; -(-9*1i-9*z_row+4*1i*z_row.^2+z_row.^3)./z_row.^4.*exp(1i*z_row)];

    max_order_allowed = size(hd,1) - 1;
    if (max_order > max_order_allowed)
        error("The maximum order should be no larger than %d\n", max_order_allowed);
    end

    n = (0:max_order).';
    hd = reshape(hd(n+1,:), size(n.*z));
end
