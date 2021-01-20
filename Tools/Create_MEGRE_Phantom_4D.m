function Phantom =  Create_MEGRE_Phantom_4D(MWF, Options)
% Gets MWF map and corresponding functions to simulate a 4D ME-GRE phantom
	if ~isfield(Options, 'FreqShift')
		% FreqShift is a struct which contains frequency shift of each water compartment
		Options.FreqShift.MW = zeros(size(MWF));
		Options.FreqShift.IE = zeros(size(MWF));
	end
	
	if ~isfield(Options, 'Times')
		% a struct containing the time steps for simulation
		error('Options argument requires "Times" field!');
	end
	
	if ~isfield(Options, 'SNR')
		Options.SNR = 0;
	end
	
	Mask = zeros(size(MWF));
	Mask(MWF > 0) = 1;
	
	DS = size(MWF); % data size
	
	% reshape data for parallel processing
	MWF = reshape(MWF, [1, numel(MWF)]);
	Mask = reshape(Mask, size(MWF));
	Options.FreqShift.MW = reshape(Options.FreqShift.MW, size(MWF));
	Options.FreqShift.IE = reshape(Options.FreqShift.IE, size(MWF));
	Phantom = zeros(numel(MWF), numel(Options.Times));
	
	MyInfo.NumData = 1;
	MyInfo.NumWaterComp = 3;
	MyInfo.SNR = Options.SNR;
	MyInfo.Times = Options.Times;
	MyInfo.T1Val = ones(1,3);
	MyInfo.TimeConstRange{1} = [10 10] * 1e-3;
	MyInfo.TimeConstRange{2} = [64 64] * 1e-3;
	MyInfo.TimeConstRange{3} = [48 48] * 1e-3;
	
	MyInfo.FreqRange{1} = [0 0];
	MyInfo.FreqRange{2} = [0 0];
	MyInfo.FreqRange{3} = [0 0];
	
	MyInfo.FractionRange{1} = [0 0];
	MyInfo.FractionRange{2} = [0 0];
	MyInfo.FractionRange{3} = [0 0];
	
	parfor i = 1:numel(MWF)
		tempInfo = MyInfo;
		if Mask(i) > 0
			tempInfo.FractionRange{1} = [1 1] * MWF(i);
			IE = 1 - MWF(i);
			tempInfo.FractionRange{2} = [1 1] * IE *2/3;
			tempInfo.FractionRange{3} = [1 1] * IE *1/3;
			tempInfo.FreqRange{1} = [1 1] * Options.FreqShift.MW(i);
			tempInfo.FreqRange{3} = [1 1] * Options.FreqShift.IE(i);
			temp = SimClass(tempInfo);
			Phantom(i,:) = temp.SimulatedData;
		end
	end
	
	Phantom = reshape(Phantom, [DS, numel(Options.Times)]);
end