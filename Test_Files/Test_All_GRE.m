addpath(genpath('~/MWI'))

ReadPath = '~/GRE/GRE_To_Do/';
SavePath = '~/GRE/GRE_Results/';

names{1} = 'GRE_HighResJ.mat';

for i = 1:length(names)
  ProcessGRE_Tukey(names{i}, ReadPath, SavePath);
end
