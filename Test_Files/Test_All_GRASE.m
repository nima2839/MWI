addpath(genpath('~/GRASE/Postprocessing'))
addpath(genpath('~/MWI'))

cd ~/GRASE/GRASE_Results/GRASE_UBC_Results
files = dir('*.mat');

%name{1} = 'Can-05-Ppm-067-M0R';


for i = 1:length(files)
  %ProcessGRASE(files(i).name);
  ReProcessGRASE(files(i).name);
end
