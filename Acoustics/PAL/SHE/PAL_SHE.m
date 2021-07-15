% Dimension: fp.theta -> fp.r -> l -> m -> n -> Gauss
function prs = PAL_SHE(transducer, ultra, audio, fp, varargin)

	if ~isequal(fp.r, sort(fp.r))
		error('The elements in fp.r must be in correct order!!\n');
    end
    
	p = inputParser;
	% the number of truncated terms
	addParameter(p, 'TruncatedTerms', 30);
    % 'Westervelt' or 'Kuznetsov'
	addParameter(p, 'Equation', 'Westervelt');
    % 'Default' or 'Mute'
	addParameter(p, 'PrintInfo', 'Default');
	%
	addParameter(p, 'isPreCalculateRadial', 0);
	addParameter(p, 'radialData', nan);
	addParameter(p, 'isFocus', 0);
	addParameter(p, 'AngularAperture', 40);
	parse(p, varargin{:});
	ip = p.Results;

    N = ip.TruncatedTerms;
	lmn = (0:N).';
	l = reshape(lmn, 1, 1, N+1);
	m = reshape(lmn, 1, 1, 1, N+1);
	n = reshape(lmn, 1, 1, 1, 1, N+1);

	% see Eq. (13) of Zhong2020JASA
	C_l = (-1).^l .* (4*l+1) / sqrt(pi) .* exp(gammaln(l+1/2) - gammaln(l+1));
	C_m = (-1).^m .* (4*m+1) / sqrt(pi) .* exp(gammaln(m+1/2) - gammaln(m+1));
	C_n = -ultra.vel0_low*conj(ultra.vel0_high)./audio.angfreq.*(4*n+1);
	% Wigner 3j symbol. see Eq. (20) of Zhong2020JASA
	wigner = 0*l.*m.*n;
	for i_l = 1:length(l)
		for j_m = 1:length(m)
			for k_n = 1:length(n)
				wigner(1, 1, i_l, j_m, k_n) = Wigner3j000(2*l(i_l),...
					2*m(j_m), 2*n(k_n));
			end
		end
	end

	% Legendre polynomials in Eq. (21) of Zhong2020JASA
	% leg_2n = permute(...
		% LegendrePolynomial(2*n(:), cos(permute(fp.theta, [5,2,3,4,1]))),...
		% [5,2,3,4,1]);
    leg_2n = permute( ...
        LegendrePolynomial(2*max(n), cos(permute(fp.theta, [5,2,3,4,1]))), [5,2,3,4,1]);
    leg_2n = leg_2n(:,:,:,:,1:2:end);

	%% Calculate the radial componentns
	if ~ip.isPreCalculateRadial
		% the points in inner region
		fp.r_inner = fp.r(fp.r<transducer.radius);
		% the points in outer region
		fp.r_outer = fp.r(fp.r>=transducer.radius);
		% Radial components. Eqs. xx of the submitted manuscript
		F.p = []; F.r = []; F.theta = [];
		if ~isempty(fp.r_inner)
			if ~strcmp(ip.PrintInfo, 'Mute')
				fprintf('============Inner region============\n')
				fprintf('Processing the inner region points...\n')
			end
			F_inner = PAL_SHE_RadialRobust(transducer, ...
				ultra, audio, ...
				fp.r_inner, l, m, n, ...
				'Equation', ip.Equation, ...
				'Region', 'inner', ...
				'PrintInfo', ip.PrintInfo, ...
				'isFocus', ip.isFocus,...
				'AngularAperture', ip.AngularAperture);
			F.p = [F.p, F_inner.p];
			F.r = [F.r, F_inner.r];
			F.theta = [F.theta, F_inner.theta];
		end
		if ~isempty(fp.r_outer)
			if ~strcmp(ip.PrintInfo, 'Mute')
				fprintf('============Outer region============\n')
				fprintf('Processing the outer region points...\n')
			end
			F_outer = PAL_SHE_RadialRobust(transducer, ....
				ultra, audio,...
				fp.r_outer, l, m, n, ...
				'Equation', ip.Equation, ...
				'Region', 'outer', ...
				'PrintInfo', ip.PrintInfo,...
				'isFocus', ip.isFocus,...
				'AngularAperture', ip.AngularAperture);
			F.p = [F.p, F_outer.p];
			F.r = [F.r, F_outer.r];
			F.theta = [F.theta, F_outer.theta];
        end
    else
        F = ip.radialData.F;
	end

	%% Calculate the velocity potentials
    % see Eq. xx of the submitted manuscript
	potential_buf = C_l .* C_m .* C_n .* wigner.^2;
	potential_p_lmn = sum(sum(potential_buf .* F.p, 4), 3) .* leg_2n;
	potential_p = sum(potential_p_lmn, 5);
	switch ip.Equation
		case 'Westervelt'
			potential = potential_p * 1.2;
			prs = 1i * 1.21 * audio.angfreq .* potential;   
		case 'Kuznetsov'
			potential_r_lmn = sum(sum(potential_buf .* F.r, 4), 3) ...
				.* leg_2n;
			potential_r = sum(potential_r_lmn, 5);
			potential_theta_lmn = sum(sum(potential_buf .* F.theta...
				.* (l.*(2*l+1) + m.*(2*m+1) - n.*(2*n+1)), 4), 3) ...
				.* leg_2n ;
			potential_theta = sum(potential_theta_lmn, 5);
			potential = potential_p .* (1.2-1) ...
				+ potential_r + potential_theta;
			
			% sound pressure and velocity of primary waves at field points
			transducer.vel0 = ultra.vel0_low;
			[prs_low, vel_low] = ...
				CircPiston_SHE(...
				ultra.low, transducer, fp, ...
				'is_cal_velocity', 1,...
                'isFocus', ip.isFocus,...
                'AngularAperture', ip.AngularAperture);
			transducer.vel0 = ultra.vel0_high;
			[prs_high, vel_high] = ...
				CircPiston_SHE(...
				ultra.high, transducer, fp, ...
				'is_cal_velocity', 1,...
                'isFocus', ip.isFocus,...
                'AngularAperture', ip.AngularAperture);
            % modifications of the Lagrangian density
			lagrangian = 1.21/2 .* (vel_high.r .* conj(vel_low.r)) ...
				- 1./(2*1.21*343^2) ...
				.* (prs_high .* conj(prs_low));
			prs = 1i*1.21*audio.angfreq .* potential ...
				- lagrangian;
	end
end
