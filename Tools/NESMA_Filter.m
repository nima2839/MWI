function Out = NESMA_Filter(input_Data, Options)
% Author: Nima
% 2020/08/22
% Ref for code:
% Bouhrara M, Maring MC, Spencer RG. A simple and fast adaptive nonlocal multispectral filtering algorithm for efficient noise reduction in magnetic resonance imaging. Magn Reson Imaging. 2019 Jan;55:133-139. doi: 10.1016/j.mri.2018.08.011. Epub 2018 Aug 24. PMID: 30149058; PMCID: PMC6242759.
% Bouhrara M, Reiter DA, Maring MC, Bonny JM, Spencer RG. Use of the NESMA Filter to Improve Myelin Water Fraction Mapping with Brain MRI. J Neuroimaging. 2018;28(6):640-649. doi:10.1111/jon.12537
%
%	Inputs:
%		Options:
%			Options.Mask
%			Options.Normalize_Flag
%			Options.SNR_Min -> SNR of last echo used in calculation to estimate Threshold parameter!
%			Options.Threshold
%			Options.Num_Channels
%			Options.Method : Distance methods -> "RED" (Relative Euclidean Distance), "RMD" (Relative Manhattan Distance), "SCD" (Relative Correlation Difference -> Testing)
%
    tic;
	sd = size(input_Data);
	
	disp('NESMA filtering started!');
	
    if ~isfield(Options, 'Mask')
		Options.Mask = input_Data(:,:,:,1) > 0;
    end
    if ~isfield(Options, 'Normalize_Flag')
         Options.Normalize_Flag = false;
    end
	if ~isfield(Options, 'Num_Channels')
		Options.Num_Channels = sd(4);
	end
	if ~isfield(Options, 'SNR_Min')
		Options.SNR_Min = 100;
	end
	if ~isfield(Options, 'Method')
		Options.Method = "RED";
	end
	if ~isfield(Options, 'Threshold')
		Options.Threshold = 2 * (Options.SNR_Min)^-2;
		if strcmp(Options.Method, "RMD")
			Options.Threshold = 0.675*Op;
		end
	end
	
	if Options.Normalize_Flag
		disp('Normalizing Data!');
		alpha = sum(input_Data.^2,4).^-0.5;
		input_Data = alpha.*input_Data;	
	end
	
	input_Data = reshape(input_Data, sd(1)*sd(2)*sd(3), sd(4));
	Data = input_Data(:,1:Options.Num_Channels);
	Mask = reshape(Options.Mask, sd(1)*sd(2)*sd(3),1);
	Out = zeros(size(input_Data));
	p = floor(length(Mask)*0.1);
	%Iterate through first three dimonsions
    disp('Iterating through voxels!')
	parfor i = 1:length(Mask)
		v = squeeze(Data(i,:));
		if Mask(i)
			if strcmp(Options.Method, "RCD")
				Distance = Data * v';
				Distance = Distance./sum(v.^2);
				Distance = abs(Distance - 1);
			else
				Diff = v - Data;
				if strcmp(Options.Method, "RMD") 
					Distance = sum(abs(Diff),2)./ sum(v);
				else % Method = "RED"
					Distance = sum(Diff.^2,2)./ sum(v)^2;
				end
			end
            
			idx = find((Distance < Options.Threshold));
			% Using Rician corrected average method!
			Out(i,:) = sqrt(max((sum(input_Data(idx,:).^2,1)./length(idx)) - var(input_Data(idx,:)), zeros(1, sd(4))));%repmat(sum(Data(idx,:),1)./length(idx), [length(idx),1]);
			%Mask(idx) = 0;	
		end
		if mod(i,p) == 0
			disp(strcat("."));
		end
	end
	Out = reshape(Out,sd);
	toc
	disp('NESMA finished!');
end

