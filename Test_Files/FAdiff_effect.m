% this is to explore FA diff effect on high SNR decay curves

Chi2Factors = 1.02;

MWF = 1;
CSF = 0;
IE = 1 - MWF - CSF;

MyInfo.NumWaterComp = 1;
MyInfo.Times = (1:32)*1e-2;
MyInfo.TimeConstRange{1} = [5 20]*1e-3;
%MyInfo.TimeConstRange{2} = [70 80]*1e-3;
%MyInfo.TimeConstRange{3} = [500 2000]*1e-3;

MyInfo.FractionRange{1}= [MWF, MWF];
%MyInfo.FractionRange{2}= [IE, IE];
%MyInfo.FractionRange{3}= [CSF, CSF];

MyInfo.T1Val = [.6 1 4.163];
MyInfo.FlipAngle = 180;
MyInfo.NumData = 500;
MyInfo.TrueFAFlag = true;
MyInfo.SNR = 1e4;

FA = [110:1:160, 165:5:180];
FAdiff = -20:0.5:20;

Results = cell(length(FA), length(FAdiff));

MyInfo.NumData = MyInfo.NumData - mod(MyInfo.NumData, maxNumCompThreads);
Noise = Create_Noise(32, MyInfo.NumData, MyInfo.SNR);

tic
for i = 1:length(FA)
	MyInfo.FlipAngle = FA(i);
	TempSim = SimClass(MyInfo);
	T_E2_E1 = mean(TempSim.SimulatedData(:,2) ./ TempSim.SimulatedData(:,1));
	SimulatedData = TempSim.SimulatedData + Noise;
	for j = 1:length(FAdiff)
		MyInfo.FlipAngle = FA(i) + FAdiff(j);
		[~, TempMaps] = SimClass.UBC_Nima_Fitting(SimulatedData, MyInfo);
		Maps.MWF = mean(TempMaps.MWF(:));
		Maps.E2_E1 = mean(TempMaps.E2_E1(:));
		Maps.T_E2_E1 = T_E2_E1;
		Results{i,j} = Maps;
	end
	disp(strcat(string(100*i/length(FA)),'%'));
end
toc

cd ~/Simulation/B1_Research/
save('FAdiff_Effect_Results_100MWF','Results','FA','FAdiff','MyInfo')






function out = Create_Noise(ETL, NumData, SNR)
	out = zeros(NumData, ETL);
	for i = 1:NumData
		out(i,:) = SimClass.ADD_Noise(zeros(1,32), SNR, 1);
	end
end