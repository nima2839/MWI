addpath(genpath('~/MWI'))
cd ~/GRASE/GRASE_Results/GRASE_ExtB1_Results

files = dir('*.mat');

for i = 1:length(files)
	Redo_B1_Ex(files(i).name);
end	