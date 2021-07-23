
function res = IntBesselType1(max_order, lower_limit, upper_limit, profile, varargin)

    p = inputParser;
    addParameter(p, 'int_num', 50);
    parse(p, varargin{:});
    ip = p.Results;

    order = (0:max_order).';
    a_row = lower_limit(:).' + 0*lower_limit.*upper_limit.*profile;
    b_row = upper_limit(:).' + 0*lower_limit.*upper_limit.*profile;
    profile_row = profile(:).' + 0*lower_limit.*upper_limit.*profile; 

    res = 0 * order .* a_row;
    gauss = GaussLegendreQuadParam(ip.int_num, 'lower_limit', a_row,...
        'upper_limit', b_row);
    for i = 1:length(order)
        res(i,:) = sum(profile_row .* besselj(order(i), gauss.zero) .* gauss.weight, 1);
    end

    res = reshape(res, size(0*order.*upper_limit.*lower_limit.*profile));
end
