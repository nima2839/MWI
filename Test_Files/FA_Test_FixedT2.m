addpath(pwd)
addpath ~/Simulation
addpath ~/MWI
addpath(genpath('~/GRASE/Postprocessing'))

MyInfo.NumWaterComp = 2;
MyInfo.Times = (1:32)*1e-2;

MyInfo.TimeConstRange{1} = [15 15]*1e-3;
MyInfo.TimeConstRange{2} = [75 75]*1e-3;
%MyInfo.TimeConstRange{3} = [300 300]*1e-3;
MyInfo.T1Val = [.6 1];
MyInfo.FractionRange{1}= [0.15,0.15];
MyInfo.FractionRange{2}= [0.85,0.85];
%MyInfo.FractionRange{3}= [0.1,0.1];
MyInfo.FlipAngle = 180;
MyInfo.NumData = 500;
MyInfo.TrueFAFlag = false;
MyInfo.SNR = 0;

FA = [150, 170];
SNR = [100,200,500];
nFA = length(FA);
nSNR = length(SNR);

for j = 1:nSNR
	temp = MyInfo;
	temp.SNR = SNR(j);
	temp_dist = cell(nFA);
	temp_maps = cell(nFA);
	temp_Tdist = cell(nFA);
	temp_Tmaps = cell(nFA);
	for i  = 1:nFA
		temp.FlipAngle = FA(i);
		temp.TrueFAFlag = false;
		a = SimClass(temp);
		[temp_dist{i},temp_maps{i}] = SimClass.UBC_Nima_Fitting(a.SimulatedData, a.MyInfo);
		a.MyInfo.TrueFAFlag = true;
		[temp_Tdist{i},temp_Tmaps{i}] = SimClass.UBC_Nima_Fitting(a.SimulatedData, a.MyInfo);
	end
	Dist{j}{:} = temp_dist;
	Maps{j}{:} = temp_maps;
	TrueFA_Dist{j}{:} = temp_Tdist;
	TrueFA_Maps{j}{:}= temp_Tmaps;
end
clear a temp
Description = 'second dim is FlipAngle, and first dim is SNR';
cd ~/Simulation/Flip_Angle_Test
save('FA_Result_150_170')%,'-v7.3')
