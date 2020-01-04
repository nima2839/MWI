function Output = Create_MWI_Phantom4D(Phantom, FA_Map, SNR)
% Set the options value for the loop
if nargin < 3
	SNR = 100;
	if nargin < 2
		FA_Map = 180 * ones(size(Phantom));
	end
end
tic;
disp('Initiating variables...')

MyInfo.NumWaterComp = 3;
MyInfo.Times = (1:32)*1e-2;

MyInfo.TimeConstRange{1} = [20 20]*1e-3;
MyInfo.TimeConstRange{2} = [80 80]*1e-3;
MyInfo.TimeConstRange{3} = [500 500]*1e-3;
MyInfo.T1Val = [.25 1 1.5];
MyInfo.FractionRange{1}= [0,0];
MyInfo.FractionRange{2}= [0.9,0.9];
MyInfo.FractionRange{3}= [0.1,0.1];
MyInfo.FlipAngle = 180;
MyInfo.NumData = 1;
MyInfo.SNR = SNR;

sd = size(Phantom);
nv = sd(1);
np = sd(2);
ns = sd(3);

Output = zeros(nv,np,ns,length(MyInfo.Times));

disp('Simulation started ...');

for i = 1:nv
	for j = 1:np
		for k = 1:ns
			if (Phantom(i,j,k) > 0) & (Phantom(i,j,k) < 1)
				MyInfo.FlipAngle = FA_Map(i,j,k);
				MyInfo.FractionRange{1} = [0, 0] + Phantom(i,j,k);
				MyInfo.FractionRange{2} = [0.9, 0.9] - Phantom(i,j,k);
				temp = SimClass(MyInfo);
				Output(i,j,k,:) = squeeze(temp.SimulatedData);
			end
		end
	end
end
disp('Finished simulating!');
toc;
end