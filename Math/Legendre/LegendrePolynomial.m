% ==============================================================================
% Legendre polynomials
% ------------------------------------------------------------------------------
% INPUT
%	max_degree  --- Maximum degree
%	x			---	argument; must be in dimensions 2, 3, ...
% ------------------------------------------------------------------------------
% OUTPUT
%	P			--- Legendre polynomial
% ==============================================================================
function P = LegendrePolynomial(max_degree, x)

    if (numel(max_degree) > 1)
        error('Requirement not satisfied: numel(max_degree) <= 1\n');
    end
    if size(x, 1) > 1
        error('Requirement not satisfied: size(x, 1) <= 1\n');
	end

	x_row = x(:).';

    degree = (0:max_degree).';
    P = 0 * degree .* x_row;

    % initial value
	P(1, :) = 1;
    if (max_degree > 0)
        P(2, :) = x_row;
    end

    % recurrence relation
	for nn = 2:max_degree
		P(nn+1,:) = ((2*nn-1).*x_row.*P(nn-1+1,:)-(nn-1)*P(nn-2+1,:))/nn;
    end

    P = reshape(P, size(degree.*x));
end
