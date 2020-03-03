addpath(genpath('~/MWI'))

ReadPath = '~/GRE/GRE_To_Do/';
SavePath = '~/GRE/GRE_Results/';

names{1} = 'Feb20_GRE2D_25Cont_Bipolar';
names{2} = 'Feb20_GRE_20cont_Monopolar';
names{3} = 'Feb20_GRE_32cont';

for i = 1:length(names)
  ProcessGRE_Tukey(names{i}, ReadPath, SavePath);
end
