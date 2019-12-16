addpath(genpath('~/MWI'))
load '~/GRASE/GRASE_Results/GRASE_Results_Special_Input_Can-05-Rrm-064-M0R.mat'

FA_Map = maps.alpha;
FA_Map(isnan(FA_Map)) = 0;
FA_Map = medfilt3(FA_Map, [5 5 3]);

tic;
[maps,distributions,~] = T2map_Nima(tf_mgrase,'T1',T1, 'Threshold', 200, 'MinRefAngle', 60, 'nAngles', 12, 'T2Range', [0.015, 2],'FlipAngleMap', FA_Map);
runtime = toc;
Deccription = 'Input FA_Map from previous results,Threshold = 200; MinRefAngle= 60, nAngles = 10, T2Range = 1e-2 to 2, and T1 = 0.2,0.5,1, Cutoff = index 40';

cd ~/GRASE/GRASE_Results

save('GRASE_Results_FA_Map_MedFilt_Can-05-Rrm-064-M0R.mat')
