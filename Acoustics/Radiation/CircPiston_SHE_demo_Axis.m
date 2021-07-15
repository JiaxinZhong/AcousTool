% =========================================================================
% 维度：fp.theta -> fp.r -> n
% =========================================================================
clear all
isCalculateVelocity = 1;
isCalculateLagrangian = 1;

wav = SoundWave(40e3);
wav.CalAbsorpCoef;
transducer.radius = 5e-2;
transducer.z = 0;
transducer.vel0 = 1;
fp.theta = 0;
fp.r = logspace(-3, 0, 5e2);
fp.phi = 0;

transducer.RayleighDist = 1/2*real(wav.wavnum) ...
    .* transducer.radius.^2;

% fp.r = linspace(0, 1.2*transducer.RayleighDist, 60);

fp.rho = fp.r .* sin(fp.theta);
fp.x = fp.rho;
fp.y = 0;
fp.z = fp.r .* cos(fp.theta);


% 使用球谐展开法进行计算
tic
[prs, vel, Lagrangian] = CircPiston_SHE(wav, transducer,...
    fp, 'isCalculateVelocity', isCalculateVelocity, ...
	'isCalculateLagrangian', isCalculateLagrangian);
sys.time_SHE = toc;
fprintf('Elapsed time for SHE: %fs\n', sys.time_SHE);

% normalization
prs = prs/2/1.21/343;
Lagrangian = Lagrangian/2/1.21/343;

spl = 20*log10(abs(prs));
Lagrangian_level = 20*log10(abs(Lagrangian));

fig = Figure;
semilogx(fp.r, spl)
fig.Init;

if isCalculateLagrangian
	fig2 = Figure;
	semilogx(fp.r, Lagrangian_level+40*log10(40));
	fig2.init;
end
