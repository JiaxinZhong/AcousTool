% =========================================================================
% 维度：fp.theta -> fp.r -> n
% =========================================================================
clear all
is_cal_velocity = 1;
wav = SoundWave(40e3);
transducer.radius = 5e-2;
transducer.z = 0;
transducer.vel0 = 1;
fp.theta = linspace(0, pi/2, 3e1).';
fp.r = linspace(0, 4.5, 3e1);
fp.isOneByOne = 1;
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
[prs, vel] = CircPiston_SHE(wav, transducer,...
    fp, 'is_cal_velocity', is_cal_velocity);
sys.time_SHE = toc;
fprintf('Elapsed time for SHE: %fs\n', sys.time_SHE);
% 若该点位于轴线，则使用精确解来计算声压

spl = 20*log10(abs(prs)/sqrt(2)/20e-6);

z_full = [flipud(fp.z(2:end,:)); fp.z];
x_full = [-flipud(fp.x(2:end,:)); fp.x];
spl_full = [flipud(spl(2:end,:)); spl];

fig_spl = Figure;
pcolor(z_full, x_full, spl_full)
fig_spl.Init;
switch dataSet
    case 'Cervenka'
        xlim([0,4]);
        ylim([-1.5,1.5])
        caxis([80,160])
    case 'Zhong'
        xlim([0,4]);
        ylim([-1.5,1.5])
        caxis([40,130])
        fig_spl.hColorbar.Ticks = [40:20:120,130];
end
fig_spl.setSize('big')
xlabel('Axial distance (m)')
ylabel('Transverse distance (m)')
fig_spl.set_ColorbarTitle('SPL (dB)')

% fig_spl.print('Acoustics/Radiation/fig/cal_CircPiston_SHE_demo_2D_spl_');

%% Plot velocity field
vel.r_full = [flipud(vel.r(2:end,:)); vel.r];
fig_vel = Figure;
pcolor(z_full, x_full, abs(vel.r_full));
fig_vel.init;
caxis([0,2])
xlim([0,0.3]);
ylim([-0.2,0.2])
