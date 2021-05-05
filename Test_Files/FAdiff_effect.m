% this is to explore FA diff effect on high SNR decay curves

Chi2Factors = 1.02;

MWF = .15;
CSF = 0;
IE = 1 - MWF - CSF;

MyInfo.NumWaterComp = 2;
MyInfo.Times = (1:32)*1e-2;
MyInfo.TimeConstRange{1} = [15 15]*1e-3;
MyInfo.TimeConstRange{2} = [75 75]*1e-3;
%MyInfo.TimeConstRange{3} = [500 2000]*1e-3;

MyInfo.FractionRange{1}= [MWF, MWF];
MyInfo.FractionRange{2}= [IE, IE];
%MyInfo.FractionRange{3}= [CSF, CSF];

MyInfo.T1Val = [.6 1 4.163];
MyInfo.FlipAngle = 180;
MyInfo.NumData = 100;
MyInfo.TrueFAFlag = true;
MyInfo.SNR = 1e4;

B1diff = -.3:0.01:0.3;
FA = [130:10:180];%, 165:5:180];
FAdiff = B1diff*180;

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
		[Dist, TempMaps] = SimClass.UBC_Nima_Fitting(SimulatedData, MyInfo);
		Maps.MWF = mean(TempMaps.MWF(:));
		Maps.E2_E1 = mean(TempMaps.E2_E1(:));
		Maps.T_E2_E1 = T_E2_E1;
		Maps.Dist = Dist;
		Results{i,j} = Maps;
	end
	disp(strcat(string(100*i/length(FA)),'%'));
end
toc

cd ~/Simulation/B1_Research/
save('FAB1diff_Effect_Results_NoCSF_Dist_fixedT2','Results','FA','FAdiff','MyInfo')






function out = Create_Noise(ETL, NumData, SNR)
	out = zeros(NumData, ETL);
	for i = 1:NumData
		out(i,:) = SimClass.ADD_Noise(zeros(1,32), SNR, 1);
	end
end