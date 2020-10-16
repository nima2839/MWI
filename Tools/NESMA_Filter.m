function Out = NESMA_Filter(Data, Mask, Normalize_Flag, Threshold)
% Author: Nima
% 2020/08/22
% Ref for code:
% Bouhrara M, Reiter DA, Maring MC, Bonny JM, Spencer RG. Use of the NESMA Filter to Improve Myelin Water Fraction Mapping with Brain MRI. J Neuroimaging. 2018;28(6):640-649. doi:10.1111/jon.12537
    tic;
	disp('NESMA filtering started!');
	if nargin < 3
		Normalize_Flag = false;
    end
    if nargin < 4
         Threshold = 0.02;
    end

	if Normalize_Flag
		disp('Normalizing Data!');
		alpha = sum(Data.^2,4).^-0.5;
		Data = alpha.*Data;	
	end
	sd = size(Data);
	Data = reshape(Data, sd(1)*sd(2)*sd(3),sd(4));
	Mask = reshape(Mask, sd(1)*sd(2)*sd(3),1);
	Out = zeros(size(Data));
	p = floor(length(Mask)*0.1);
	%Iterate through first three dimonsions
    disp('Iterating through voxels!')
	parfor i = 1:length(Mask)
		if Mask(i)
			%Diff = bsxfun(@minus, Data, Data(i,:));
			%RMD = sum(abs(Diff),2)./ sum(Data(i,:));
			v = squeeze(Data(i,:));
			RMD = Data * v';
            if ~Normalize_Flag
                RMD = RMD./sum(v.^2);
            end
            RMD = sqrt(abs(RMD - 1));
			idx = find((RMD < Threshold));
			Out(i,:) = sum(Data(idx,:),1)./length(idx);%repmat(sum(Data(idx,:),1)./length(idx), [length(idx),1]);
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

