% Dimension order: fp.theta -> fp.r

[transducer, ultra, audio] = LoadPalParam('set_name', 'Zhong');
% 'Westervelt' or 'Kuznetsov'
eqn = 'Westervelt';
fp.set_name = '2D_few';

%% field points settings
%	Warning: the memory sizes increase if you assign too many field points
%	Hint: If you want to calculate the sound fields at many field points
%		It is better to pre-calcualte the radial components F by 
%		the routine cal_radial_component.m and store them for further
%		calculations.
switch fp.set_name
    % By choosing this set, only the SPL at a single point will be
    % calculated
    case 'SinglePoint'
        fp.r = 0.2;
        fp.theta = 0;
        
    % By choosing this set, the axial SPL will show at a few points
    case '1D_few'
        fp.r = [linspace(0, 0.0199,30), linspace(0.0205, 0.05, 30),...
			linspace(0.051, .1, 25), linspace(.101,0.3,25)]; 
        fp.theta = 0;
        
    % By choosing this set, the axial SPL will show at many points
    case '1D_many'
        fp.r = [linspace(0, 0.0199,78), linspace(0.0205, 0.05, 30),...
            linspace(0.051, .1, 30), linspace(.101,0.3,30)]; 
        fp.theta = 0;

	% By choosing this set, a two-dimensional sound fields will show
	case '2D'		
%         fp.r = [linspace(0, 0.0199,30), linspace(0.0205, 0.05, 30),...
% 			linspace(0.051, .1, 25), linspace(.101,0.3,25)]; 
        fp.r = linspace(0,0.099,5);

		fp.theta = linspace(0,pi/2,20).';
    case '2D_few'
        fp.r = [0.1,0.2];
        fp.theta = [0, pi/2, pi/3].';
	case '2D_210715B'		
        % fp.r = [linspace(0, 0.0199,30), linspace(0.0205, 0.05, 30),...
            % linspace(0.051, .1, 25), linspace(.101,0.3,25)]; 
        fp.r = linspace(0, 4, 200);
		fp.theta = linspace(0,pi/2,100).';
end
fp.phi = 0;
fp.x = fp.r .* sin(fp.theta);
fp.z = fp.r .* cos(fp.theta);

% The truncated terms, see Eq. xx of the submitted manuscript
N = 70;

% the main function
timer_total = tic;
prs = PAL_SHE(transducer, ultra, audio, fp, ...
	'Equation', eqn, 'TruncatedTerms', N);
spl = prs2spl(prs);
time_total = toc(timer_total);
fprintf('The total elapsed time: %f s.\n', time_total)

if numel(prs) == 1
    fprintf('===========Results============\n');
	fprintf('Pressure: %f + %fi; Abs: %f; SPL: %f\n',...
		real(prs), imag(prs), abs(prs), spl);
end

if size(prs,1)>1 && size(prs,2) > 1
    x_new = [flipud(-fp.x(2:end,:)); fp.x];
    z_new = [flipud(fp.z(2:end,:)); fp.z];
    spl_new = [flipud(spl(2:end,:)); spl];
    fig2D = Figure;
    pcolor(z_new, x_new, spl_new);
    fig2D.Init;
%     ylim([-0.15,0.15])
%     xlim([0, 0.3])
%     caxis([60,80])
end
