% ==========================================================================
% - j'_n(ka) * h'_n(kr) * int_a^b h_n(xi) / h'n(ka) * P(kd / xi) * xi dxi
% - Dimension: order => kr => xi
% --------------------------------------------------------------------------
% INPUT
% 
% CONCLUSION
%   - Gauss Legendre quadrature is found to be better than the Simpson rule.
% ==========================================================================
function res = IntSphBesselLegendreTypeC(max_order, kr, ka, kd, ...
    int_a, int_b, int_seg_num, varargin)

    p = inputParser;
    addParameter(p, 'method', 'GaussLegendreQuad');
    parse(p, varargin{:});
    ip = p.Results;

    switch ip.method
        case 'SimpsonRule'
            % todo
            % xi = permute(linspace(int_a, int_b, int_seg_num*2+1).', [3,2,1]);
        case 'GaussLegendreQuad'
            param = GaussLegendreQuadParam(int_seg_num, 'lower_limit', int_a,...
                'upper_limit', int_b);
            xi = permute(param.zero, [3,2,1]);
    end

    jh = SphBessel1DTimesSphHankel1D(max_order, ka, kr);
    hh = SphHankel1DDivSphHankel1(max_order, ka, xi);
    leg = LegendrePolynomial(max_order, kd./xi);

    switch ip.method
        case 'SimpsonRule'
            % todo
            % res = SimpsonRule(jh .* leg .* xi, xi(2)-xi(1), 3);
        case 'GaussLegendreQuad'
            res = sum(jh ./ hh .* leg .* xi .* permute(param.weight, [3,2,1]), 3);
    end
end
