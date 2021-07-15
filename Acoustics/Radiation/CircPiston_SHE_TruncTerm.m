% =========================================================================
% 维度：fp.theta -> fp.r -> n
% =========================================================================
clear all

wav = SoundWave(60e3);
transducer.radius = 0.1;

fp.r = linspace(0,10,50);
fp.theta = 0;
fp.phi = 0;

tic
prs = CircPiston_SHE(wav, transducer,fp);
sys.time_SHE = toc;
fprintf('Elapsed time for SHE: %fs\n', sys.time_SHE);
