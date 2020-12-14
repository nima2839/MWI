function Output = Create_MWI_Phantom4D(Phantom, FA_Map, MyInfo)
% Set the options value for the loop
if nargin < 3
	SNR = 100;
	if nargin < 2
		FA_Map = 180 * ones(size(Phantom));
	end
end
tic;
disp('Initiating variables...')

if ~isfield(MyInfo, 'NumWaterComp')
	MyInfo.NumWaterComp = 3;
	MyInfo.Times = (1:32)*1e-2;

	MyInfo.TimeConstRange{1} = [20 20]*1e-3;
	MyInfo.TimeConstRange{2} = [80 80]*1e-3;
	MyInfo.TimeConstRange{3} = [500 500]*1e-3;
	MyInfo.T1Val = [.6 1 4];
	MyInfo.FractionRange{1}= [0,0];
	MyInfo.FractionRange{2}= [0.9,0.9];
	MyInfo.FractionRange{3}= [0.1,0.1];
	MyInfo.FlipAngle = 180;
	MyInfo.NumData = 1;
	MyInfo.SNR = 300;
end

MyInfo.NumData = 1; % Optimized for this function

sd = size(Phantom);
nv = sd(1);
np = sd(2);
ns = sd(3);

IE = MyInfo.FractionRange{2};

Output = zeros(nv,np,ns,length(MyInfo.Times));

disp('Simulation started ...');

parfor i = 1:nv
	temp_out = zeros(np, ns, sd(4));
	for j = 1:np
		for k = 1:ns
			if ~isnan(Phantom) & (Phantom(i,j,k) > 0) & (Phantom(i,j,k) < 1)
				MyInfo.FlipAngle = FA_Map(i,j,k);
				MyInfo.FractionRange{1} = [0, 0] + Phantom(i,j,k);
				MyInfo.FractionRange{2} = IE - Phantom(i,j,k);
				temp = SimClass(MyInfo);
				temp_out(j,k,:) = squeeze(temp.SimulatedData);
			end
		end
	end
	Output(i,:,:,:) = temp_out;
end
disp('Finished simulating!');
toc;
end