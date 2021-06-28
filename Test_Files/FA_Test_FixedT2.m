addpath(pwd)
addpath ~/Simulation
addpath ~/MWI
addpath(genpath('~/GRASE/Postprocessing'))

MyInfo.NumWaterComp = 2;
MyInfo.Times = (1:32)*1e-2;

IE = SimClass.Create_Guassian_Dist(75e-3); % intra/extra-cellular water 
MW = SimClass.Create_Guassian_Dist(15e-3); % myelin water

temp.T2Values = [MW.T2Values, IE.T2Values];
temp.Weights = [MWFs(i) * MW.Weights, (1 - MWFs(i)) * IE.Weights];
MyInfo.T2Dist = temp;
MyInfo.T1Val = [0.6 * ones(size(MW.Weights)), ones(size(1 * IE.Weights))]; % setting T1 parameter
clear temp

%MyInfo.TimeConstRange{1} = [15 15]*1e-3;
%MyInfo.TimeConstRange{2} = [75 75]*1e-3;
%MyInfo.TimeConstRange{3} = [300 300]*1e-3;
%MyInfo.T1Val = [.6 1];
%MyInfo.FractionRange{1}= [0.15,0.15];
%MyInfo.FractionRange{2}= [0.85,0.85];
%MyInfo.FractionRange{3}= [0.1,0.1];
MyInfo.FlipAngle = 180;
MyInfo.NumData = 1000;
MyInfo.TrueFAFlag = false;
MyInfo.SNR = 0;

FA = [0.7:0.05:1.3]*180;
SNR = 50:50:1e3;
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
	%Dist{j}{:} = temp_dist;
	Maps{j}{:} = temp_maps;
	%Dist_FA{j}{:} = temp_Tdist;
	Maps_FA{j}{:}= temp_Tmaps;
end
clear a temp temp_Tdist temp_dist temp_Tmaps temp_maps
Description = 'second dim is FlipAngle, and first dim is SNR';
cd ~/Simulation/Flip_Angle_Test
save('FA_Result_SNR_Sweep','-v7.3')
