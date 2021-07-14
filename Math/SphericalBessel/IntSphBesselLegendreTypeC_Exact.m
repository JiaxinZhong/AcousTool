% ==========================================================================
% - j'_n(ka) * h'_n(kr) * int_a^b h_n(xi) / h'n(ka) * P(kd / xi) * xi dxi
% - Dimension: order => kr => xi
% - Exact
% --------------------------------------------------------------------------
% Calculate using the mathematica code:
%   n = 0;
%   Simplify[FunctionExpand[
%   Integrate[
%   D[SphericalBesselJ[n, ka], ka]/
%   D[SphericalHankelH1[n, ka], ka] SphericalHankelH1[n, x] D[
%    SphericalHankelH1[n, kr], kr] LegendreP[n, kd/x] x, x]]]
% ==========================================================================
function res = IntSphBesselLegendreTypeC_Exact(max_order, kr, ka, kd, ...
    int_a, int_b, varargin)

    res = func1(ka, kr, kd, int_b) - func1(ka, kr, kd, int_a);
    res = [res; func2(ka, kr, kd, int_b) - func2(ka, kr, kd, int_a)];
    res = [res; func3(ka, kr, kd, int_b) - func3(ka, kr, kd, int_a)];

    MAX_ORDER = size(res,1);
    if (max_order > MAX_ORDER)
        error('The maximum order should not be greater than %d.\n', MAX_ORDER);
    end

end

function res = func1(ka, kr, kd, xi)
    res = -(1i+kr) .* (ka.*cos(ka) - sin(ka)) ...
        ./(1i+ka) ./ kr.^2 ...
        .* exp(-1i*(ka-kr-xi));
end

function res = func2(ka, kr, kd, xi)
    res = (1i*exp(-1i*(ka-kr-xi)).*kd.*(-2+2*1i*kr+kr.^2) ...
        .*(2*ka.*cos(ka)+(ka.^2-2).*sin(ka)))...
        ./ ((-2+2*1i*ka+ka.^2).*kr.^3.*xi);
end

function res = func3(ka, kr, kd, xi)
    res = (exp(-1i*(ka-kr-xi)) .* (-9*1i-9*kr+4*1i*kr.^2+kr.^3) ...
        .* (-3*kd.^2*(1i+xi)+xi.^2.*(3*1i+xi)) ...
        .* (ka.*(ka.^2-9).*cos(ka) + (9-4*ka.^2).*sin(ka)))...
        ./ (2*(-9*1i-9*ka+4*1i.*ka.^2+ka.^3).*kr.^4.*xi.^3);
end
