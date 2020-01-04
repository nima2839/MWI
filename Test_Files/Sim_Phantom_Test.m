addpath(genpath('~/MWI'))

Phantom_MWF = phantom3d(100);
FA_Map = Create_FA_Map(size(Phantom_MWF));
FA_Map = Normalize_FA_Map(FA_Map);
Phantom_MWI = Create_MWI_Phantom4D(Phantom_MWF, FA_Map, 45);
tic;
[maps,distributions,~] = T2map_Nima(Phantom_MWI, 'Threshold', 0, 'MinRefAngle', 60, 'nAngles', 12, 'T2Range', [0.015, 2]);
runtime = toc;

Final_MWI = squeeze(sum(distributions(:,:,:,1:40),4)./sum(distributions,4));

Deccription = ',Threshold = 200; MinRefAngle= 60, nAngles = 10, T2Range = 1e-2 to 2, and T1 = 1, Cutoff = index 40';

cd ~/Simulation/Phantom

save('Sim_Phantom_Results')
