addpath(genpath('~/MWI'))
addpath(genpath('~/GRASE/Postprocessing'))

cd ~/Simulation/B1_Research/

%%
tic

MyInfo.NumWaterComp = 3;
MyInfo.Times = (1:32)*1e-2;
MyInfo.TimeConstRange{1} = [3 25]*1e-3;
MyInfo.TimeConstRange{2} = [60 80]*1e-3;
MyInfo.TimeConstRange{3} = [500 2000]*1e-3;
MyInfo.T1Val = [.6 1 4.163];
MyInfo.FractionRange{1}= [0.07,0.07];
MyInfo.FractionRange{2}= [0.83,0.83];
MyInfo.FractionRange{3}= [0.1,0.1];



MyInfo.FlipAngle = 180;
MyInfo.NumData = 5e3;
MyInfo.TrueFAFlag = false;
MyInfo.SNR = 0;
%%
FA = 110:180;
SNR = [100:50:500, 750, 1e3, 1e4];
nFA = length(FA);
nSNR = length(SNR);
Dist =  cell(nSNR,nFA);
Maps = Dist;
TrueFA_Dist = Dist;
TrueFA_Maps = Dist;
parfor j = 1:nSNR
	temp = MyInfo;
	temp.SNR = SNR(j);
	temp_dist = cell(1,nFA);
	temp_maps = cell(1,nFA);
	temp_Tdist = cell(1,nFA);
	temp_Tmaps = cell(1,nFA);
	for i  = 1:nFA
		temp.FlipAngle = FA(i);
		temp.TrueFAFlag = false;
		a = SimClass(temp);
		[temp_dist{1,i}, temp_maps{1,i}] = UBC_Nima_Fitting(a);
		a.MyInfo.TrueFAFlag = true;
		[temp_Tdist{1,i}, temp_Tmaps{1,i}] = UBC_Nima_Fitting(a);
	end
	Dist(j,:) = temp_dist;
	Maps(j,:) = temp_maps;
	TrueFA_Dist(j,:) = temp_Tdist;
	TrueFA_Maps(j,:)= temp_Tmaps;
	disp(j)
end
clear a temp
runtime = toc;
Description = 'second dim is FlipAngle, and first dim is SNR';
%%
cd ~/Simulation/B1_Research/
save('B1_Sim_Result_M4','-v7.3')
disp('Done!')
