wav = SoundWave(64e3);

radius = 0.2;
fp.z = linspace(0, 10, 3e3).';

[prs, vel] = CircPiston_OnAxis(radius, wav.wavnum, fp.z);

prs = prs/1.21/343;

fig = Figure;
subplot(2,1,1)
plot(fp.z, abs(prs));
fig.Init;
title('pressuer/rho0c0v0')

subplot(2,1,2)
semilogx(fp.z, abs(vel.z));
fig.Init;
title('velocity_z')

% test plane waves
err = abs(prs./vel.z-1);
fig2 = Figure;
plot(fp.z, err)
fig2.Init
