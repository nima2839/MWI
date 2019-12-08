MyInfo.NumData = 1e5;
MyInfo.Times = (1:32) * 1e-2;
MyInfo.NumWaterComp = 3;
MyInfo.TimeConstRange{1} = [5 40] * 1e-3;
MyInfo.TimeConstRange{2} = [42 120] * 1e-3;
MyInfo.TimeConstRange{3} = [80 500] * 1e-3;
MyInfo.T1Val{1} = .65;
MyInfo.T1Val{2} = 1;
MyInfo.T1Val{3} = 1.2;
MyInfo.FractionRange{1} = [0 0.5];
MyInfo.FractionRange{2} = [0.5 1];
MyInfo.FractionRange{3} = [0 0.5];
FA = linspace(60,180,50);
SNR = linspace(40,300,n);
Data = zeros(length(SNR),length(FA),MyInfo.NumData,length(MyInfo.Times));
Labels =  zeros(length(SNR),length(FA),MyInfo.NumData);

parfor i = 1:length(SNR)
    for j = 1:length(FA)
        tempInfo = MyInfo;
        tempInfo.SNR = SNR(i);
        tempInfo.FlipAngle = FA(j);
        sim = SimClass(tempInfo);
        Data(i,j,:,:) = sim.SimulatedData;
        Labels(i,j,:) = sim.Compartment_Fraction_Map{1}/(sim.Compartment_Fraction_Map{1}+sim.Compartment_Fraction_Map{2}+sim.Compartment_Fraction_Map{3});
    end
end

if n > 5
	save('TrainData.mat','-v7.3')
else
	save('TestData.mat','-v7.3')
end
