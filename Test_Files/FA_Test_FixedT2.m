addpath(pwd)
addpath ~/Simulation
addpath ~/MWI
addpath(genpath('~/GRASE/Postprocessing'))

MyInfo.NumWaterComp = 3;
MyInfo.Times = (1:32)*1e-2;

MyInfo.TimeConstRange{1} = [20 20]*1e-3;
MyInfo.TimeConstRange{2} = [80 80]*1e-3;
MyInfo.TimeConstRange{3} = [300 300]*1e-3;
MyInfo.T1Val = [.8 1 1.5];
MyInfo.FractionRange{1}= [0.15,0.15];
MyInfo.FractionRange{2}= [0.75,0.75];
MyInfo.FractionRange{3}= [0.1,0.1];
MyInfo.FlipAngle = 180;
MyInfo.NumData = 500;
MyInfo.TrueFAFlag = false;
MyInfo.SNR = 0;

FA = 80:2:180;
SNR = 30:10:100;
nFA = length(FA);
nSNR = length(SNR);

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
cd ~/Simulation/Flip_Angle_Test
save('FA_Result_FixedT2','-v7.3')
