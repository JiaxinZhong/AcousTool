% =========================================================================
% Calculate the first radial component
% -------------------------------------------------------------------------
% 输入
%	r				---	the radial coordinate of the field point;
%						must be in the 2nd dimension
%	l, m, n			--- orders; must be in the 3rd, 4th, and 5th dimensions
%	ultra			---	ultrasounds
%	audio			---	audio sound
% =========================================================================
function F1 = ...
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
    
	hbar_2n_kar = permute(SphHankel1Norm(2*n(:), ka*r), [5,2,3,4,1]);

	% 使用GaussLegendre进行积分
	GaussNumber_hJ = 60;
	GaussNumber_jH = 60;
	gauss = GaussLegendreQuadParam(ip.GaussNumber);
	zero = reshape(gauss.zero(1:ip.GaussNumber,1), 1, 1, 1, 1, 1, ip.GaussNumber);
	weight = reshape(gauss.weight(1:ip.GaussNumber,1), 1, 1, 1, 1, 1, ip.GaussNumber);

	switch ip.Region
		case 'inner'
			rv1 = ((r-0).*zero+(r-0))/2;
			% Using the normalized spherical Bessel functions, see Eq. (28) of Zhong2020JASA
			C1bar = (ka*rv1).^(2*n)./((ka*r).^(2*n+1)) ./ (4*n+1)/1i .* (r-0)/2;
		case 'outer'
			rv1 = ((a-0).*zero+(a+0))/2;
			C1bar = (ka*rv1).^(2*n)./((ka*r).^(2*n+1)) ./ (4*n+1)/1i .* (a-0)/2;
    end
    
    

	jbar_2n_karv1 = permute(SphBessel1Norm(2*n(:), ka*rv1), [5,2,3,4,1,6]);
    if ip.isFocus
        hJ_2l_rv1 = permute(...
            cal_hJ_norm_focus(2*l(:), k1*rv1, 0, k1*rv1, k1*focalLength,...
            'GaussNumber', GaussNumber_hJ), [3,2,1,4,5,6]);
        jH_2l_rv1 = permute(cal_jH_norm_focus(2*l(:), k1*rv1, k1*rv1, k1*a, k1*focalLength,...
            'GaussNumber', GaussNumber_jH), [3,2,1,4,5,6]);
        hJ_2m_rv1 = permute(cal_hJ_norm_focus(2*m(:), k2*rv1, 0, k2*rv1, k2*focalLength, ...
           'GaussNumber', GaussNumber_hJ), [4,2,3,1,5,6]);
        jH_2m_rv1 = permute(cal_jH_norm_focus(2*m(:), k2*rv1, k2*rv1, k2*a, k2*focalLength,...
           'GaussNumber', GaussNumber_jH), [4,2,3,1,5,6]);
    else
        hJ_2l_rv1 = permute(cal_hJ_norm(2*l(:), k1*rv1, 0, k1*rv1, ...
            'GaussNumber', GaussNumber_hJ), [3,2,1,4,5,6]);
        jH_2l_rv1 = permute(cal_jH_norm(2*l(:), k1*rv1, k1*rv1, k1*a,...
            'GaussNumber', GaussNumber_jH), [3,2,1,4,5,6]);
        hJ_2m_rv1 = permute(cal_hJ_norm(2*m(:), k2*rv1, 0, k2*rv1, ...
           'GaussNumber', GaussNumber_hJ), [4,2,3,1,5,6]);
        jH_2m_rv1 = permute(cal_jH_norm(2*m(:), k2*rv1, k2*rv1, k2*a, ...
           'GaussNumber', GaussNumber_jH), [4,2,3,1,5,6]);
    end

	R_2l_rv1 = hJ_2l_rv1 + jH_2l_rv1;
	R_2m_rv1 = hJ_2m_rv1 + jH_2m_rv1;
	switch ip.Equation
		case 'Westervelt'
			;
		case 'Kuznetsov'
            if ip.isFocus
                potential_r.R_2l_rv1 = permute(...
                    cal_hDJ_norm_focus(2*l(:), k1*rv1, 0, k1*rv1, k1*focalLength) ...
                    + cal_jDH_norm_focus(2*l(:), k1*rv1, k1*rv1, k1*a, k1*focalLength), [3,2,1,4,5,6]);
                potential_r.R_2m_rv1 = permute(...
                    cal_hDJ_norm_focus(2*m(:), k2*rv1, 0, k2*rv1, k2*focalLength) ...
                    + cal_jDH_norm_focus(2*m(:), k2*rv1, k2*rv1, k2*a, k2*focalLength), [4,2,3,1,5,6]);
            else
                potential_r.R_2l_rv1 = permute(...
                    cal_hDJ_norm(2*l(:), k1*rv1, 0, k1*rv1) ...
                    + cal_jDH_norm(2*l(:), k1*rv1, k1*rv1, k1*a), [3,2,1,4,5,6]);
                potential_r.R_2m_rv1 = permute(...
                    cal_hDJ_norm(2*m(:), k2*rv1, 0, k2*rv1) ...
                    + cal_jDH_norm(2*m(:), k2*rv1, k2*rv1, k2*a), [4,2,3,1,5,6]);
            end
	end

	F1_buf = weight .* C1bar .* jbar_2n_karv1 .* hbar_2n_kar;
	F1.p = sum(F1_buf .* ka^3 .* rv1.^2 .* R_2l_rv1 .* conj(R_2m_rv1), 6);
		
	switch ip.Equation
		case 'Westervelt'
			;
		case 'Kuznetsov'
			F1.r = sum(F1_buf .* ka^3 .* rv1.^2 ...
				.* potential_r.R_2l_rv1 .* conj(potential_r.R_2m_rv1), 6);
			F1.theta = sum(F1_buf .* ka^3 ./ k1 ./ conj(k2) ...
				.* R_2l_rv1 .* conj(R_2m_rv1), 6);
	end
end
