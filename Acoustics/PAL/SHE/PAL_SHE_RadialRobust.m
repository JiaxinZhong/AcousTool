% =========================================================================
% Calculate the radial components
% -------------------------------------------------------------------------
% 输入
%	r				---	the radial coordinate of the field point;
%						must be in the 2nd dimension
%	l, m, n			--- orders; must be in the 3rd, 4th, and 5th dimensions
%	ultra			---	ultrasounds
%	audio			---	audio sound
% =========================================================================
function [F, F1, F2, F3] = ...
    PAL_SHE_RadialRobust(transducer, ultra, audio, r, l, m, n, varargin)

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
	addParameter(p, 'PrintInfo', 'Default');
	% Focusing PAL
	addParameter(p, 'isFocus', 0);
	addParameter(p, 'AngularAperture', 40);
	parse(p, varargin{:});
	ip = p.Results;

	% 使用GaussLegendre进行积分
	GaussNumber1 = 80;
	GaussNumber2 = 80;
	GaussNumber3 = 200;
    
	%% calculate the first integral
	%	see Eqs. (23) and (24) of Zhong2020JASA
	tic
	F1 = PAL_SHE_Radial1(transducer, ...
		ultra, audio, r, l, m, n,...
		'Region', ip.Region, ...
		'Equation', ip.Equation, ...
        'GaussNumber', GaussNumber1, ...
        'isFocus', ip.isFocus,...
        'AngularAperture', ip.AngularAperture);
	buf = toc;
	if ~strcmp(ip.PrintInfo, 'Mute')
		fprintf('Elapsed time on preparing the first integral: %f s\n', buf);
	end

	tic
	F2 = PAL_SHE_Radial2(transducer, ...
		ultra, audio, r, l, m, n, ...
		'Region', ip.Region, ...
		'Equation', ip.Equation, ...
        'GaussNumber', GaussNumber2, ...
        'isFocus', ip.isFocus,...
        'AngularAperture', ip.AngularAperture);
	buf = toc;
	if ~strcmp(ip.PrintInfo, 'Mute')
		fprintf('Elapsed time on preparing the second integral: %f s\n', buf);
	end

	%% Calculate the 3rd integral
	% see Eqs. (23) and (24) of Zhong2020JASA
	% when l, m, and n is smaller than 20, it is better to use the complex 
	%	method to calculate the integral. see Eq. (32) of Zhong2020JASA
	timer3 = tic;
	F3 = PAL_SHE_Radial3Robust(transducer, ...
        ultra, audio, ...
		r, l, m, n, ...
		'Region', ip.Region, ...
		'PrintInfo', 'Mute', ...
        'Equation', ip.Equation,...
        'GaussNumber', GaussNumber3, ...
        'isFocus', ip.isFocus,...
        'AngularAperture', ip.AngularAperture);
    
	timer3 = toc(timer3);
	if ~strcmp(ip.PrintInfo, 'Mute')
		fprintf('Elapsed time on preparing the third integral: %f s\n', timer3);
	end

	%% performan integral
	timer3 = tic;

	F1.p(1, r==0, :, :, :) = 0;
	F2.p(1, r==0, :, :, 2:end) = 0;
	F3.p(1, r==0, :, :, 2:end) = 0;
	F.p = F1.p + F2.p + F3.p;
	F.r = [];
	F.theta = [];
	switch ip.Equation
		case 'Westervelt'
			;
		case 'Kuznetsov'
			if ip.isOrigin
				F1.r = 0;
				F1.theta = 0;
			end
			F1.r(1, r==0, :, :, :) = 0;
			F2.r(1, r==0, :, :, 2:end) = 0;
			F3.r(1, r==0, :, :, 2:end) = 0;
			F1.theta(1, r==0, :, :, :) = 0;
			F2.theta(1, r==0, :, :, 2:end) = 0;
			F3.theta(1, r==0, :, :, 2:end) = 0;
			F.r = F1.r + F2.r + F3.r;
			F.theta = F1.theta + F2.theta + F3.theta;
	end

	timer3_val = toc(timer3);
	if ~strcmp(ip.PrintInfo, 'Mute')
		fprintf('Elapsed time on preparing the integral: %f s\n', buf);
	end
end
