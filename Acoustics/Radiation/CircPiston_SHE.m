% =========================================================================
% Calculate the sound radiation from a circular piston using the 
%	spherical harmonics expansion (SHE) method
% -------------------------------------------------------------------------
% NOTE
%   - dim1 and dim2: fp.theta .* fp.r
%   - dim3: n
% -------------------------------------------------------------------------
% INPUT
%	fp.r, fp.theta, fp.phi		---	Spherical coordinates of field points
% -------------------------------------------------------------------------
% OUTPUT
%	prs, vel					--- sound pressure and velocity
%   lagrangian  --- Lagrangian density
% =========================================================================
function [prs, vel, Lagrangian] = CircPiston_SHE(wav, ...
    transducer, fp, varargin)

	% the amplitude of the vibration velocity on the transducer surface
    if ~isfield(transducer, 'vel0')
        transducer.vel0 = 1;
    end
	k = wav.wavnum;
	a = transducer.radius;

	p = inputParser;
	% Maximum truncated terms of the spherical harmonics
	addParameter(p, 'N_MAX', round(k*a));
	addParameter(p, 'is_cal_velocity', 0);
	addParameter(p, 'is_cal_lagrangian', 0);
	% for focusing PAL
	addParameter(p, 'isFocus', 0);
	addParameter(p, 'AngularAperture', 40);
	parse(p, varargin{:});
	ip = p.Results;
	if ip.is_cal_lagrangian
		ip.is_cal_velocity = 1;
	end
    
    if ip.isFocus
        FocalLength = a * cot(ip.AngularAperture/180*pi/2);
    end
	% 1: the elements in fp should be calculated one by one
	if ~isfield(fp, 'isOneByOne')
		fp.isOneByOne = 0;
	end

    N_MAX_outer = ip.N_MAX*1.5;
	N_MAX_inner = ip.N_MAX*100;

	% Separate the inner and outer points
	fp.r(fp.r==0) = a*1e-10;
	holder_fp = 0 * fp.theta .* fp.r;
	if fp.isOneByOne
		fp.r = fp.r + holder_fp;
		fp.theta = fp.theta + holder_fp;
	end
	fp.theta_inner = fp.theta;
	fp.theta_outer = fp.theta;
	idx_inner = fp.r<transducer.radius;
	idx_outer = fp.r>=transducer.radius;
	fp.r_inner = fp.r(idx_inner);
	fp.r_outer = fp.r(idx_outer);
	if fp.isOneByOne
		fp.theta_inner = fp.theta(idx_inner);
		fp.theta_outer = fp.theta(idx_outer);
	end

	% pressure
	prs = holder_fp;

	%% Pressure in the outer region
	if ~isempty(fp.r_outer)
        n_outer = permute((0:N_MAX_outer).', [3,2,1]);
        Cn_outer = (-1).^n_outer .* (4*n_outer+1) ...
            / sqrt(pi) .* exp(...
			gammaln(n_outer+1/2)-gammaln(n_outer+1));

		leg_2n_outer = permute(...
			LegendrePolynomial(2*max(n_outer), ...
			cos(permute(fp.theta_outer, [3,2,1]))), [3,2,1]);
        leg_2n_outer = leg_2n_outer(:,:,1:2:end);
        if ip.isFocus
            R_2n_outer = ...
                permute(...
                cal_hJ_norm_focus(2*n_outer(:), ...
                k*permute(fp.r_outer, [3,2,1]), 0, k*a, ...
                k*FocalLength), ...
                [3,2,1]);
        else
            R_2n_outer = ...
                permute(...
                cal_hJ_norm(2*n_outer(:), ...
                k*permute(fp.r_outer, [3,2,1]), 0, k*a), ...
                [3,2,1]);
        end
		prs_n_outer = 1.21*343*transducer.vel0 ...
			* Cn_outer .* R_2n_outer .* leg_2n_outer;
		prs(idx_outer | holder_fp) = sum(prs_n_outer,3);
    end

    %% Pressure in  inner region
	if ~isempty(fp.r_inner)
		n_inner = permute((0:N_MAX_inner).', [3,2,1]);
		C_n_inner = (-1).^n_inner .* (4*n_inner+1) ...
			/ sqrt(pi) .* exp(gammaln(n_inner+1/2)-gammaln(n_inner+1));

		% leg_2n_inner = ...
			% permute(...
			% LegendrePolynomial(2*n_inner(:), ...
			% cos(permute(fp.theta_inner,[3,2,1]))), [3,2,1]);
		leg_2n_inner = ...
			permute(...
			LegendrePolynomial(2*max(n_inner), ...
			cos(permute(fp.theta_inner,[3,2,1]))), [3,2,1]);
        leg_2n_inner = leg_2n_inner(:,:,1:2:end);
		if ip.isFocus
			buf1 = permute(...
				cal_hJ_norm_focus(2*n_inner(:), ...
				k*permute(fp.r_inner,[3,2,1]), ...
				0, k*permute(fp.r_inner,[3,2,1]), ...
				k*FocalLength), [3,2,1]);
			buf2 = permute(...
				cal_jH_norm_focus(2*n_inner(:), ...
				k*permute(fp.r_inner,[3,2,1]), ...
				k*permute(fp.r_inner,[3,2,1]), k*a, k*FocalLength), [3,2,1]);
		else
			buf1 = permute(...
				cal_hJ_norm(2*n_inner(:), ...
				k*permute(fp.r_inner,[3,2,1]), ...
				0, k*permute(fp.r_inner,[3,2,1])), [3,2,1]);
			buf2 = permute(...
				cal_jH_norm(2*n_inner(:), ...
				k*permute(fp.r_inner,[3,2,1]), ...
				k*permute(fp.r_inner,[3,2,1]), k*a), [3,2,1]);
		end
		% special case for fp.r==0
		% buf1((fp.r==0 | holder_fp) | 0*n_inner) = 0;
		% buf2((fp.r==0 | holder_fp) & n_inner>0) = 0;

		R_2n_inner = buf1 + buf2;
		prs_n_inner = C_n_inner.*R_2n_inner.*leg_2n_inner;
		prs(idx_inner | holder_fp) = 1.21*343*transducer.vel0 ...
			.* sum(prs_n_inner,3);
	end

	vel = [];
	if ~ip.is_cal_velocity
		return
    end
    
	%% calculate the velocity
	vel.r = holder_fp;
	vel.theta = holder_fp;
	if ~isempty(fp.r_outer)
		% legD_2n_outer = permute(...
			% LegendrePolynomialD(2*n_outer(:), ...
			% cos(permute(fp.theta_outer,[3,2,1]))), [3,2,1]) ...
			% .* (-sin(fp.theta_outer));
		legD_2n_outer = permute(...
			LegendrePolynomialD(2*max(n_outer), ...
			cos(permute(fp.theta_outer,[3,2,1]))), [3,2,1]) ...
			.* (-sin(fp.theta_outer));
        legD_2n_outer = legD_2n_outer(:,:,1:2:end);

		if ip.isFocus
			RD_2n_outer = permute(...
				cal_hDJ_norm_focus(2*n_outer(:), ...
				k*permute(fp.r_outer,[3,2,1]), 0, k*a, ...
				k*FocalLength), [3,2,1]);
		else
			RD_2n_outer = permute(...
				cal_hDJ_norm(2*n_outer(:), ...
				k*permute(fp.r_outer,[3,2,1]), 0, k*a), [3,2,1]);
		end
		vel.r(idx_outer | holder_fp) = -1i * transducer.vel0 ...
            .* sum(Cn_outer .* RD_2n_outer .* leg_2n_outer, 3);
		
		vel.theta(idx_outer | holder_fp) = -1i * transducer.vel0 ...
			.* sum(Cn_outer .* R_2n_outer ...
            .* legD_2n_outer, 3) ./ (k.*(fp.r_outer));
	end

	%% Velocity in the inner region
	if ~isempty(fp.r_inner)
		if ip.isFocus
			buf1 = permute(...
				cal_hDJ_norm_focus(n_inner(:)*2, ...
				k*permute(fp.r_inner,[3,2,1]), ...
				0, k*permute(fp.r_inner,[3,2,1]), ...
				k*FocalLength), [3,2,1]);
			buf2 = permute(...
				cal_jDH_norm_focus(n_inner(:)*2, ...
				k*permute(fp.r_inner,[3,2,1]), ...
				k*permute(fp.r_inner,[3,2,1]), k*a, ...
				k*FocalLength), [3,2,1]);
		else
			buf1 = permute(...
				cal_hDJ_norm(n_inner(:)*2, ...
				k*permute(fp.r_inner,[3,2,1]), ...
				0, k*permute(fp.r_inner,[3,2,1])), [3,2,1]);
			buf2 = permute(...
				cal_jDH_norm(n_inner(:)*2, ...
				k*permute(fp.r_inner,[3,2,1]), ...
				k*permute(fp.r_inner,[3,2,1]), k*a), [3,2,1]);
		end

		% buf1((fp.r==0 | holder_fp) | 0*n_inner) = 0;
        RD_2n_inner = buf1 + buf2;
		vel.r(idx_inner | holder_fp) = -1i * transducer.vel0 ...
            .* sum(C_n_inner .* RD_2n_inner .* leg_2n_inner, 3);

		legD_2n_inner = permute(...
			LegendrePolynomialD(2*max(n_inner), ...
			cos(permute(fp.theta_inner,[3,2,1]))), [3,2,1]) ...
			.* (-sin(fp.theta_inner));
        legD_2n_inner = legD_2n_inner(:,:,1:2:end);
		vel.theta(idx_inner | holder_fp) = -1i * transducer.vel0 ...
            .* sum(C_n_inner .* R_2n_inner .* legD_2n_inner, 3) ...
			./ k ./ fp.r_inner;
	end
	
	vel.phi = vel.r *0;

	vel.x = vel.r .* sin(fp.theta) .* cos(fp.phi) ...
		+ vel.theta .* cos(fp.theta) .* cos(fp.phi) ...
		+ vel.phi .* (-sin(fp.phi));
	vel.y = vel.r .* sin(fp.theta) .* sin(fp.phi) ...
		+ vel.theta .* cos(fp.theta) .* sin(fp.phi) ...
		+ vel.phi .* cos(fp.phi);
	vel.z = vel.r .* cos(fp.theta) ...
		+ vel.theta .* (-sin(fp.theta));

	if check_infnan(vel.r, 'Mode', 'Mute') || check_infnan(vel.theta, 'Mode', 'Mute')
		warning('There are inf or nan!\n');
	end

	%% calculate Lagaragian density
	Lagrangian = 1.21/2 .* (abs(vel.r).^2 + abs(vel.theta).^2) ...
		- abs(prs).^2 / 2 / 1.21/343^2;
end
