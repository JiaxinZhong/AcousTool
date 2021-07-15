% =========================================================================
% Calculate the second radial component
% -------------------------------------------------------------------------
% 输入
%	r				---	the radial coordinate of the field point;
%						must be in the 2nd dimension
%	l, m, n			--- orders; must be in the 3rd, 4th, and 5th dimensions
%	ultra			---	ultrasounds
%	audio			---	audio sound
% =========================================================================
function F2 = ...
    cal_radial_component1(transducer, ultra, audio, r, l, m, n, varargin)

	%% check the dimensions
	if ~isscalar(r) && (size(r,1) > 1 ...
			|| size(r,3) > 1 || size(r,4) > 1 ...
			|| size(r,5) > 1)
		error('r must be in the dimension 2!\n')
	end
	if ~isscalar(l) && (size(l,1) > 1 ...
			|| size(l,2) > 1 || size(l,4) > 1 ...
			|| size(l,5) > 1)
		error('l must be in the dimension 3!\n')
	end
	if ~isscalar(m) && (size(m,1) > 1 ...
			|| size(m,2) > 1 || size(m,3) > 1 ...
			|| size(m,5) > 1)
		error('m must be in the dimension 4!\n')
	end
	if ~isscalar(n) && (size(n,1) > 1 ...
			|| size(n,2) > 1 || size(n,3) > 1 ...
			|| size(n,4) > 1)
		error('n must be in the dimension 5!\n')
	end

	%% process options
	p = inputParser;
	addParameter(p, 'Equation', 'Westervelt');
	addParameter(p, 'Region', 'outer');
	addParameter(p, 'isOrigin', 0);
	addParameter(p, 'GaussNumber', 100);
	addParameter(p, 'isFocus', 0);
	addParameter(p, 'AngularAperture', 40);
	parse(p, varargin{:});
	ip = p.Results;

	% wavnumbers
	k1 = ultra.high.wavnum;
	k2 = ultra.low.wavnum;
	ka = audio.wavnum;
	a = transducer.radius;
    if ip.isFocus
        focalLength = a * cot(ip.AngularAperture/180*pi/2);
    end

	h_2n_kar = permute(SphHankel1(2*n(:), ka*r), [5,2,3,4,1]);
	switch ip.Region
		case 'inner'
			jbar_2n_kar = reshape(SphBessel1Norm(2*n(:), ka*r).', 1, length(r), 1, 1, length(n));
		case 'outer'
			;
    end
    if ip.isFocus
        J_2l_0_k1a = reshape(cal_sphBessel_int_Gauss_0n_focus(l(end), ...
            0, k1*a, k1*focalLength, 60), 1, 1, length(l));
        J_2m_0_k2a = reshape(cal_sphBessel_int_Gauss_0n_focus(m(end), ...
            0, k2*a, k2*focalLength, 60), 1, 1, 1, length(m));
    else
        J_2l_0_k1a = reshape(cal_sphBessel_int_Gauss_0n(l(end), ...
            0, k1*a, 60), 1, 1, length(l));
        J_2m_0_k2a = reshape(cal_sphBessel_int_Gauss_0n(m(end), ...
            0, k2*a, 60), 1, 1, 1, length(m));
    end
	% 使用GaussLegendre进行积分
	GaussNumber_hJ = 60;
	GaussNumber_jH = 60;
	ip.GaussNumber = 80;
	gauss = GaussLegendreQuadParam(ip.GaussNumber);
	zero = reshape(gauss.zero(1:ip.GaussNumber,1), 1, 1, 1, 1, 1, ip.GaussNumber);
	weight = reshape(gauss.weight(1:ip.GaussNumber,1), 1, 1, 1, 1, 1, ip.GaussNumber);

	switch ip.Region
		case 'inner'
			rv2 = ((a-r) .* zero + (a+r))/2;
			C2 = (ka*r).^(2.*n) ./ ((ka.*rv2).^(2*n+1)) ./ 1i ./ (4*n+1) .* (a-r)/2;
		case 'outer'
			rv2 = ((r-a) .* zero + (r+a))/2;
			C2 = (r-a)/2;
	end
	switch ip.Region
		case 'inner'
            if ip.isFocus
                hJ_2l_rv2 = permute(cal_hJ_norm_focus(2*l(:), k1*rv2, 0, k1*rv2, k1*focalLength,...
                    'GaussNumber', GaussNumber_hJ), [3,2,1,4,5,6]);
                jH_2l_rv2 = permute(cal_jH_norm_focus(2*l(:), k1*rv2, k1*rv2, k1*a, k1*focalLength,...
                    'GaussNumber', GaussNumber_jH), [3,2,1,4,5,6]);
                hJ_2m_rv2 = permute(cal_hJ_norm_focus(2*m(:), k2*rv2, 0, k2*rv2, k2*focalLength,...
                    'GaussNumber', GaussNumber_hJ), [4,2,3,1,5,6]);
                jH_2m_rv2 = permute(cal_jH_norm_focus(2*m(:), k2*rv2, k2*rv2, k2*a, k2*focalLength,...
                    'GaussNumber', GaussNumber_jH), [4,2,3,1,5,6]);

            else
                hJ_2l_rv2 = permute(cal_hJ_norm(2*l(:), k1*rv2, 0, k1*rv2, ...
                    'GaussNumber', GaussNumber_hJ), [3,2,1,4,5,6]);
                jH_2l_rv2 = permute(cal_jH_norm(2*l(:), k1*rv2, k1*rv2, k1*a, ...
                    'GaussNumber', GaussNumber_jH), [3,2,1,4,5,6]);
                hJ_2m_rv2 = permute(cal_hJ_norm(2*m(:), k2*rv2, 0, k2*rv2, ...
                    'GaussNumber', GaussNumber_hJ), [4,2,3,1,5,6]);
                jH_2m_rv2 = permute(cal_jH_norm(2*m(:), k2*rv2, k2*rv2, k2*a, ...
                    'GaussNumber', GaussNumber_jH), [4,2,3,1,5,6]);
            end
			hbar_2n_karv2 = permute(SphHankel1Norm(2*n(:), ka*rv2), [5,2,3,4,1,6]);
			R_2l_rv2 = hJ_2l_rv2 + jH_2l_rv2;
			R_2m_rv2 = hJ_2m_rv2 + jH_2m_rv2;
		case 'outer'
			R_2l_rv2 = permute(SphHankel1(2*l(:), k1*rv2), [3,2,1,4,5,6]) .* J_2l_0_k1a;
			R_2m_rv2 = permute(SphHankel1(2*m(:), k2*rv2), [4,2,3,1,5,6]) .* J_2m_0_k2a;
			j_2n_karv2 = permute(SphBessel1(2*n(:), ka*rv2), [5,2,3,4,1,6]);
	end
	switch ip.Equation
		case 'Westervelt'
			;
		case 'Kuznetsov'
			switch ip.Region
				case 'inner'
                    if ip.isFocus
                        potential_r.R_2l_rv2 = ...
                            permute(cal_hDJ_norm_focus(2*l(:), k1*rv2, 0, k1*rv2, k1*focalLength) ...
                            + cal_jDH_norm_focus(2*l(:), k1*rv2, k1*rv2, k1*a, k1*focalLength), [3,2,1,4,5,6]);
                        potential_r.R_2m_rv2 = ...
                            permute(cal_hDJ_norm_focus(2*m(:), k2*rv2, 0, k2*rv2, k2*focalLength) ...
                            + cal_jDH_norm_focus(2*m(:), k2*rv2, k2*rv2, k2*a, k2*focalLength), [4,2,3,1,5,6]);
                    else
                        potential_r.R_2l_rv2 = ...
                            permute(cal_hDJ_norm(2*l(:), k1*rv2, 0, k1*rv2) ...
                            + cal_jDH_norm(2*l(:), k1*rv2, k1*rv2, k1*a), [3,2,1,4,5,6]);
                        potential_r.R_2m_rv2 = ...
                            permute(cal_hDJ_norm(2*m(:), k2*rv2, 0, k2*rv2) ...
                            + cal_jDH_norm(2*m(:), k2*rv2, k2*rv2, k2*a), [4,2,3,1,5,6]);
                    end
				case 'outer'
					potential_r.R_2l_rv2 = ...
						permute(SphHankel1D(2*l(:), k1*rv2), [3,2,1,4,5,6]) ...
						.* J_2l_0_k1a;
					potential_r.R_2m_rv2 = ...
						permute(SphHankel1D(2*m(:), k2*rv2), [4,2,3,1,5,6]) ...
						.* J_2m_0_k2a;
			end
			if check_infnan(potential_r.R_2l_rv2, 'Mode', 'mute') ...
				|| check_infnan(potential_r.R_2m_rv2, 'Mode', 'mute') 
				error('存在inf或nan!\n');
			end
	end
	switch ip.Region
		case 'inner'
			F2_buf = weight .* C2 .* jbar_2n_kar .* hbar_2n_karv2;
		case 'outer'
			F2_buf = weight .* C2 .* j_2n_karv2 .* h_2n_kar;
	end
	F2.p = sum(F2_buf .* ka^3 .* rv2.^2 .* R_2l_rv2 .* conj(R_2m_rv2), 6);

	switch ip.Equation
		case 'Westervelt'
			;
		case 'Kuznetsov'
			F2.r = sum(F2_buf .* ka^3 .* rv2.^2 ...
				.* potential_r.R_2l_rv2 .* conj(potential_r.R_2m_rv2), 6);
			F2.theta = sum(F2_buf .* ka^3 ./ k1 ./ conj(k2) ...
				.* R_2l_rv2 .* conj(R_2m_rv2), 6);
	end
end
