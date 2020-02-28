addpath(genpath('~/MWI'))

ReadPath = '~/GRE/GRE_To_Do/';
SavePath = '~/GRE/GRE_Results/';

names{1} = 'Feb20_GRE2D_25Cont_Bipolar.mat';
names{2} = 'Feb20_GRE_32cont.mat';

for i = 1:length(names)
  ProcessGRE(names{i}, ReadPath, SavePath);
end
