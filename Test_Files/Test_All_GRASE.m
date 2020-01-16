addpath(genpath('~/GRASE/Postprocessing'))
addpath(genpath('~/MWI'))

cd ~/GRASE/GRASE_To_Do/
files = dir('GRASE_Results_*.mat');

for i = 1:length(files)
  ProcessGRASE(files(i).name);
end
