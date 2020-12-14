function Output = Create_MWI_Phantom4D(Phantom, FA_Map, MyInfo)
% Set the options value for the loop
if nargin < 3
	SNR = 100;
	if nargin < 2
		FA_Map = 180 * ones(size(Phantom));
	end
end

disp('Initiating variables...')

if ~isfield(MyInfo, 'NumWaterComp')
	MyInfo.NumWaterComp = 3;
	MyInfo.Times = (1:32)*1e-2;
	MyInfo.IE = SimClass.Create_Guassian_Dist(75e-3); % intra/extra-cellular water 
	MyInfo.MW = SimClass.Create_Guassian_Dist(15e-3); % myelin water
	MyInfo.T2Dist.T2Values = [MyInfo.MW.T2Values, MyInfo.IE.T2Values];
	MyInfo.T1Val = [.6*ones(size(MyInfo.MW.Weights)), ones(size(MyInfo.IE.Weights))];
	MyInfo.FlipAngle = 180;
	MyInfo.NumData = 1;
	MyInfo.SNR = 300;
end

MyInfo.NumData = 1; % Optimized for this function

sd = size(Phantom);
nv = sd(1);
np = sd(2);
ns = sd(3);
ne = length(MyInfo.Times)



Output = zeros(nv*np*ns,ne);
Mask = zeros(nv*np*ns,1);
Phantom = reshape(Phantom, size(Mask));

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
		tempInfo.FlipAngle = FA_Map(i);
		MW = Phantom(i);
		IE = 1-MW;
		tempInfo.T2Dist.Weights = [MW * MyInfo.MW.Weights, IE * MyInfo.IE.Weights];
		tempInfo.NumData = length(idx);
		temp = SimClass(tempInfo);

		Output(idx,:) = temp.SimulatedData;
		Mask(idx) = false;
	end
end

Output = reshape(Output, [sd,ne]);
disp('Finished simulating!');
end