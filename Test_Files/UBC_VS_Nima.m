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
MyInfo.FlipAngle = 160;
MyInfo.NumData = 500;
MyInfo.TrueFAFlag = false;
MyInfo.SNR = 60;

FA = 80:2:180;
nFA = length(FA);

SNR = 40:5:100;

for j = 1:length(SNR)
	MyInfo.SNR = SNR(j);
	tic;
	parfor i  = 1:nFA
		temp = MyInfo;
		temp.FlipAngle = FA(i);
		a = SimClass(temp);
		[Ndist{i}{j},Nmaps{i}{j}] = UBC_Nima_Fitting(a);
		[dist{i}{j},maps{i}{j}] = UBC_Fitting(a);
	end
	time = toc;
	fprintf('ETA: %f minutes! \n', (length(SNR) - j)* time/60)
end

cd ~/Simulation/UBC_VS_Nima
save('UBCNIMA_Results')
