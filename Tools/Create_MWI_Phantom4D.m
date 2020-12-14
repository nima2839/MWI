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



Output = zeros(nv*np*ns,length(MyInfo.Times));
Phantom = reshape(Phantom, size(Output));
Mask = zeros(nv*np*ns,1);
Mask(~isnan(Phantom(:,1))) = true;
Mask(Phantom < 0) = false;
Mask(Phantom > 1) = false;
Phantom = round(Phantom,2); % rounding to 2 decimal point
FA_Map = reshape(FA_Map, size(Mask));

disp('Simulation started ...');



for i = 1:size(Output,1)
	if Mask(i)
		idx = find(Phantom == Phantom(i));
		tempInfo = MyInfo;
		tempInfo.FlipAngle = FA_Map(i,j,k);
		tempInfo.FractionRange{1} = [0, 0] + Phantom(i,j,k);
		tempInfo.FractionRange{2} = IE - Phantom(i,j,k);
		tempInfo.NumData = length(idx);
		temp = SimClass(tempInfo);
		sd = size(temp.SimulatedData);
		Data = reshape(temp.SimulatedData, sd(1)*sd(2)*sd(3), sd(4))
		Output(idx,:) = Data(1:length(idx), :);
		Mask(idx) = false;
	end
end

Output = reshape(Output, sd);
disp('Finished simulating!');
toc;
end