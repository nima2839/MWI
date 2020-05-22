% this is a parametr sweep code for B1 effect research
addpath(genpath('~/MWI'))
addpath(genpath('~/GRASE/Postprocessing'))

Chi2Factors = [1.005, 1.010, 1.015, 1.02];
MWFs = (1:2:21) / 100;

Models = cell(length(Chi2Factors), length(MWFs));

tic

for i = 1:length(Chi2Factors)
	for j = 1:length(MWFs)
		Models{i,j} = Template_Model(MWFs(j), Chi2Factors(i));
		disp(strcat(string((i+j)/ (length(Chi2Factors) + length(MWFs))) + '% ...'))
	end
end

runtime = toc;

save('B1_Sim_Result_PS','-v7.3')
disp('All Done!')