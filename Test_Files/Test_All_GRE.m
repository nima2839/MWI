addpath(genpath('~/MWI'))

ReadPath = '~/GRE/GRE_To_Do/';
SavePath = '~/GRE/GRE_Results/';

names{1} = 'GRE_ROFlip_Brain.mat';
names{2} = 'GRE_ROFlip_Brain_PhaseAP.mat';
names{3} = 'GRE_ROFlip_Brain_RFSpoil.mat';

for i = 1:length(names)
  ProcessGRE_Tukey(names{i}, ReadPath, SavePath);
end
