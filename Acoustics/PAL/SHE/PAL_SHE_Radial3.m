% =========================================================================
% Calculate the third radial component
% -------------------------------------------------------------------------
% 输入
%	r				---	场点坐标；必须在第2维
%	l, m, n			--- 阶数；必须在第3、4、5维
%	ultra			---	超声
%	audio			---	音频声
%	a				---	活塞半径
% =========================================================================
function F3 = ...
    cal_radial_component3(transducer, ultra, audio, r, l, m, n, varargin)

	%% check the dimension
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
	addParameter(p, 'GaussNumber', 100);
	addParameter(p, 'PrintInfo', 'Default');
    addParameter(p, 'Method', 'Gauss');
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

	j_2n_kar = permute(SphBessel1(2*n(:), ka*r), [5,2,3,4,1]);
    if ip.isFocus
        J_2l_0_k1a = permute(cal_sphBessel_int_Gauss_0n_focus(l(end), ...
            0, k1*a, k1*focalLength, 60), [3,2,1]);
        J_2m_0_k2a = permute(cal_sphBessel_int_Gauss_0n_focus(m(end), ...
            0, k2*a, k2*focalLength, 60), [4,2,3,1]);
    else
        J_2l_0_k1a = permute(cal_sphBessel_int_Gauss_0n(l(end), ...
            0, k1*a, 60), [3,2,1]);
        J_2m_0_k2a = permute(cal_sphBessel_int_Gauss_0n(m(end), ...
            0, k2*a, 60), [4,2,3,1]);
    end
    J_2l_0_k1a = J_2l_0_k1a(:, :, l+1);
    J_2m_0_k2a = J_2m_0_k2a(:, :, :, m+1);
    
	gauss = GaussLegendreQuadParam(ip.GaussNumber);
	zero = permute(gauss.zero, [6,2,3,4,5,1]);
	weight = permute(gauss.weight, [6,2,3,4,5,1]);

	tic
	switch ip.Method
		case 'Complex'
			k_ave = 1/3 * abs(k1+conj(k2)+ka);
			k1 = k1/k_ave;
			k2 = k2/k_ave;
			ka = ka/k_ave;
	end
	switch ip.Region
		case 'inner'
			switch ip.Method
				case 'Complex'
					rv = a*k_ave + 1i* tan(pi/4*(zero+1));
				case 'Gauss'
					rv = a + tan(pi/4*(zero+1));
			end
		case 'outer'
			switch ip.Method
				case 'Complex'
					rv = r*k_ave + 1i* tan(pi/4*(zero+1));
				case 'Gauss'
					rv = r + tan(pi/4*(zero+1));
			end
	end
	exp_all = exp(1i*(k1-conj(k2)+ka) .* rv);
	switch ip.Method
		case 'Complex'
			C = ka^3 .* rv.^2 * 1i .* exp_all ...
				.* (sec(pi/4*(zero+1))).^2 * pi/4;
		case 'Gauss'
			C = ka^3 * rv.^2 .* exp_all ...
				.* (sec(pi/4*(zero+1))).^2 * pi/4 ;
	end

	hhat_2l_k1rv = ...
		permute(SphHankel1Scaled(2*l(:), k1*rv), ...
		[3,2,1,4,5,6]); 
	hhat_2m_k2rv = ...
		permute(SphHankel1Scaled(2*m(:), k2*conj(rv)),...
		[4,2,3,1,5,6]);
	hhat_2n_karv = ...
		permute(SphHankel1Scaled(2*n(:), ka*rv),...
		[5,2,3,4,1,6]);

	switch ip.Equation
		case 'Kuznetsov'
			buf = permute(SphHankel1Scaled((0:2*l(end)+1).', ...
				k1*rv), [3,2,1,4,5,6]);
			hhatD_2l_k1rv = buf(1,:,2*l+1,1,1,:) ./ (k1*rv).*2.*l ...
				- buf(1,:,2*l+2,1,1,:);
			buf = permute(SphHankel1Scaled((0:2*m(end)+1).',...
				k2*conj(rv)), [4,2,3,1,5,6]); 
			hhatD_2m_k2rv = buf(1,:,1,2*m+1,1,:) ./ (k2*conj(rv)) ...
				.* 2 .* m - buf(1,:,1,2*m+2,1,:);
	end
	buf = toc;
	if ~strcmp(ip.PrintInfo, 'Mute')
		fprintf('Elapsed time on preparing the third integral: %f s\n', buf);
	end

	% performan integral
	tic 
	F3_buf = C .* weight .* J_2l_0_k1a .* conj(J_2m_0_k2a) ...
		.* hhat_2n_karv .* j_2n_kar;
	F3.p = sum(F3_buf .* hhat_2l_k1rv .* conj(hhat_2m_k2rv), 6);

	switch ip.Equation
		case 'Westervelt'
			;
		case 'Kuznetsov'
			F3.r = sum(F3_buf .* hhatD_2l_k1rv .* conj(hhatD_2m_k2rv), 6);
			F3.theta = sum(F3_buf ./ (k1.*rv) ./ (conj(k2).*rv) ...
				.* hhat_2l_k1rv .* conj(hhat_2m_k2rv), 6);
	end
		
	buf = toc;
	if ~strcmp(ip.PrintInfo, 'Mute')
		fprintf('Elapsed time on preparing the integral: %f s\n', buf);
	end
end
