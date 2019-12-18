addpath(genpath('~/MWI'))
load('~/GRASE/GRASE_Results/GRASE_Results_Special_Input_Can-05-Rrm-064-M0R.mat','tf_mgrase')

tic;
[maps,distributions,~] = T2map_Nima(tf_mgrase, 'Threshold', 200, 'MinRefAngle', 60, 'nAngles', 12, 'T2Range', [0.015, 2]);
runtime = toc;

Final_MWI = squeeze(sum(distributions(:,:,:,1:40),4)./sum(distributions,4));

Deccription = ',Threshold = 200; MinRefAngle= 60, nAngles = 10, T2Range = 1e-2 to 2, and T1 = 1, Cutoff = index 40';

cd ~/GRASE/GRASE_Results

save('GRASE_Results_B1Effect_Can-05-Rrm-064-M0R.mat')
