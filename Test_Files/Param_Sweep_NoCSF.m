% this is a parametr sweep code for B1 effect research
addpath(genpath('~/MWI'))
addpath(genpath('~/GRASE/Postprocessing'))
clear
clc

Chi2Factors = 1.02;
MWFs = [0,.15,.3];



MyInfo.NumWaterComp = 2;
MyInfo.Times = (1:32)*1e-2;
MyInfo.TimeConstRange{1} = [15 15]*1e-3;
MyInfo.TimeConstRange{1} = [75 75]*1e-3;
MyInfo.T1Val = [1];

MyInfo.FlipAngle = 180;
MyInfo.NumData = 5e3;
MyInfo.TrueFAFlag = false;
MyInfo.SNR = 0;



Models = cell(length(Chi2Factors), length(MWFs));

tic

for i = 1:length(Chi2Factors)
	MyInfo.Chi2Factor = Chi2Factors(i);
	for j = 1:length(MWFs)
		IE = 1 - MWFs(j);
		%MyInfo.FractionRange{1}= [MWFs(j), MWFs(j)];
		MyInfo.FractionRange{1}= [IE, IE];
		Models{i,j} = Template_Model(MyInfo);
		disp(strcat(string(100*(i+j)/ (length(Chi2Factors) + length(MWFs))) + '% ...'))
	end
end

runtime = toc;
cd ~/Simulation/B1_Research/
save('B1_Sim_Result_PS_FixedT2','-v7.3')
disp('All Done!')
clear
Summerize_Sim_Results('B1_Sim_Result_PS_FixedT2')