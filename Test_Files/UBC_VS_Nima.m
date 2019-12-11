addpath(genpath('~/MWI'))
addpath(genpath('~/GRASE/Postprocessing'))

MyInfo.NumWaterComp = 3;
MyInfo.Times = (1:32)*1e-2;

MyInfo.TimeConstRange{1} = [20 20]*1e-3;
MyInfo.TimeConstRange{2} = [80 80]*1e-3;
MyInfo.TimeConstRange{3} = [500 500]*1e-3;
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
nSNR = length(SNR)

time = 0;

for j = 1:nSNR
	if time > 0
		fprintf('ETA: %f minutes! \n', time)
	else
		disp('Process Started...')
	end
	MyInfo.SNR = SNR(j);
	tic;
	parfor i  = 1:nFA
		temp = MyInfo;
		temp.FlipAngle = FA(i);
		a = SimClass(temp);
		[Ndist{i}{j}, Nmaps{i}{j}] = UBC_Nima_Fitting(a);
		[Dist{i}{j}, Maps{i}{j}] = UBC_Fitting(a);
	end
	time = toc;
	time = (length(SNR) - j)* time/60;
end

clear time i j

% Re-structure data:
for i = 1:nSNR
    for j = 1:nFA
			temp1(i,j) = Maps{j}{i};
			temp2{i,j} = Dist{j}{i};
      temp3(i,j) = Nmaps{j}{i};
			temp4{i,j} = Ndist{j}{i};
    end
end
Maps = temp1;
Dist = temp2;
Nmaps = temp3;
Ndist = temp4;

clear temp1 temp2 temp3 temp4 i j

cd ~/Simulation/UBC_VS_Nima
save('UBCNIMA_Results')
