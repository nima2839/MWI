addpath(genpath('~/MWI'))
cd ~/GRASE/GRASE_To_Do/

files = dir('*.mat');

for i = 1:length(files)
	Redo_B1_Ex(files(i).name);
end
