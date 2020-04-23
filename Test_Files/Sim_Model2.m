addpath(genpath('~/MWI'))
addpath(genpath('~/GRASE/Postprocessing'))

cd ~/Simulation/B1_Research/
load('Models') % predefined models


MyInfo.Times = (1:32)*1e-2;

MyInfo.T2Dist.T2Values = time;
MyInfo.T2Dist.Weights = M2;

MyInfo.T1Val = ones(1,length(M2));
MyInfo.T1Val(1:15) = 0.5;
MyInfo.T1Val(33:end) = 2.5;
MyInfo.FlipAngle = 180;
MyInfo.NumData = 500;
MyInfo.TrueFAFlag = false;
MyInfo.SNR = 0;

FA = 110:180;
SNR = [40:2:60, 80];
nFA = length(FA);
nSNR = length(SNR);
Dist =  cell(nSNR,nFA);
Maps = Dist;
parfor j = 1:nSNR
	temp = MyInfo;
	temp.SNR = SNR(j);
	temp_dist = cell(nFA);
	temp_maps = cell(nFA);
	temp_Tdist = cell(nFA);
	temp_Tmaps = cell(nFA);
	for i  = 1:nFA
		temp.FlipAngle = FA(i)
		temp.TrueFAFlag = false;
		a = SimClass(temp);
		[temp_dist{i},temp_maps{i}] = UBC_Fitting(a);
		a.MyInfo.TrueFAFlag = true;
		[temp_Tdist{i},temp_Tmaps{i}] = UBC_Fitting(a);
	end
	Dist{j}{:} = temp_dist;
	Maps{j}{:} = temp_maps;
	TrueFA_Dist{j}{:} = temp_Tdist;
	TrueFA_Maps{j}{:}= temp_Tmaps;
end
clear a temp
Description = 'second dim is FlipAngle, and first dim is SNR';
cd ~/Simulation/B1_Research/
save('B1_Sim_Result_M2','-v7.3')
