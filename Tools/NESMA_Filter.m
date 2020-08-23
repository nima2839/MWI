function Out = NESMA_Filter(Data, Mask, Normalize_Flag, Threshold)
% Author: Nima
% 2020/08/22
% Ref for code:
% Bouhrara M, Reiter DA, Maring MC, Bonny JM, Spencer RG. Use of the NESMA Filter to Improve Myelin Water Fraction Mapping with Brain MRI. J Neuroimaging. 2018;28(6):640-649. doi:10.1111/jon.12537
    if nargin < 3
		Normalize_Flag = false;
    end
    if nargin < 4
         Threshold = 0.05;
    end

	if Normalize_Flag
		Data = Data./repmat(sum(abs(Data),4), [1,1,1,size(Data,4)]);	
	end
	sd = size(Data);
	Data = reshape(Data, sd(1)*sd(2)*sd(3),sd(4));
	Mask = reshape(Mask, sd(1)*sd(2)*sd(3),1);
	Out = zeros(size(Data));
	
	%Iterate through first three dimonsions
	parfor i = 1:length(Mask)
		if Mask(i)
		Out(i,:) = Apply_Filt(Data(i,:));
		end
	end
	Out = reshape(Out,sd);
	function output = Apply_Filt(a)
		Diff = bsxfun(@minus, Data,a);
		RMD = sum(abs(Diff),2);
		idx = find((RMD < Threshold) & (Mask > 0));
		output = sum(Data(idx,:),1)./length(idx);
	end
end

