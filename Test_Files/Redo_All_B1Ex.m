cd ~/GRASE/GRASE_Results/B1Ex

files = dir('*.mat');

for i = 1:length(files)
	Redo_B1_Ex(files(i).name);
end	