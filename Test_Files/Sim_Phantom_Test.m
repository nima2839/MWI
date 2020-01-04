addpath(genpath('~/MWI'))

MyPhantom = phantom3d(100);
FA_Map = Create_FA_Map(size(MyPhantom));
FA_Map = Normalize_FA_Map(FA_Map);

tic;
[maps,distributions,~] = T2map_Nima(MyPhantom, 'Threshold', 200, 'MinRefAngle', 60, 'nAngles', 12, 'T2Range', [0.015, 2]);
runtime = toc;

Final_MWI = squeeze(sum(distributions(:,:,:,1:40),4)./sum(distributions,4));

Deccription = ',Threshold = 200; MinRefAngle= 60, nAngles = 10, T2Range = 1e-2 to 2, and T1 = 1, Cutoff = index 40';

cd ~/Simulation/Phantom

save('Sim_Phantom_Results')
