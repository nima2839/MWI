function [Output, Gp, Gv, Gs] = LFG_Correction(Mag, Phase, Info)
%   Applies LFG correction and removes points that are highly effected by LFG
%   Author: Nima 
%   Date: 02/2019
% Inputs:
%	Phase data must be in radian
%	Info (strutct):	
%	- Method: an integer indicating the method to use:
%		1. "Du" -> source: 'Multi-echo acquistion of MR Angiography..." by Du and colleagues
%		2. "Hwang" -> source: 'In vivo multi-silice mapping of myelin ..." by Hwang and colleagues
%	- EchoSpacing (s)
%	- FirstTE (s)
%	- Vox = [dx dy dz] (m)
%	- EchoIndexes: a vector indicating indexes of echos which are used in the correction
%	- Threshold: the threshold for the correction factor to apply to data


% Code definitions
DuMethod = 1;
HwangMethod = 2;
Nima = 3;
SizeData = size(Mag);
Gamma = 42.575e6;
% Laying code constraints on inputs
	if nargin < 3
        error('Function requires 2 inputs!');
    end
	
	if length(size(Mag)) ~= 4 || length(size(Phase)) ~= 4
		error('Input data must be 4D!');
	end
	
	if isfield(Info,'EchoSpacing') == 0
		error('Info.EchoSpacing can not be undefined!');
	end
	
	if isfield(Info,'FirstTE') == 0
		error('Info.FirstTE can not be undefined!');
	end
	
	if isfield(Info,'Vox') == 0
		error('Info.Vox can not be undefined!');
	end

	if isfield(Info,'EchoIndexes') == 0
        Info.EchoIndexes = [5 6];
    end
	
	if isfield(Info,'Method') == 0
        Info.Method = DuMethod;
    end
	
	if isfield(Info,'Threshold') == 0
        Info.Threshold = 1e-8;
    end
	
	
	disp('Starting LFG correction process...')	
	tic
	Info2 = Info;
	Info2.deltaTE = Info.EchoSpacing * (Info.EchoIndexes(2) - Info.EchoIndexes(1));
	Info2.Phase = Phase;
	Info2.Mag = Mag;
	ComplexData = Mag.*exp(i*Phase);
	disp('EchoIndexes')
	disp(Info.EchoIndexes)
	Echo1 = ComplexData(:,:,:,Info.EchoIndexes(1));
	Echo2 = ComplexData(:,:,:,Info.EchoIndexes(2));
	[Gp, Gv, Gs] = Find_LFG(Echo1, Echo2, Info2);
	Gp = Gp.* Info.Mask;
	Gv = Gv.* Info.Mask;
	Gs = Gs.* Info.Mask;
	
	TE = zeros(SizeData(4));
	for k = 1:SizeData(4)
		TE(k) = Info.FirstTE + (k-1) * Info.EchoSpacing;
	end
	
	CorDenomWeights = zeros(SizeData);
	for k = 1:SizeData(4)
		CorDenomWeights(:,:,:,k) = 	sinc(Gamma * Gp * Info.Vox(1)  * TE(k) / 2) .* ...
						sinc(Gamma * Gv * Info.Vox(2)  * TE(k) / 2) .* ...
						sinc(Gamma * Gs * Info.Vox(3)  * TE(k) / 2) .* double(Info.Mask);
	end
%	for p = 1:SizeData(1)
%		for j = 1:SizeData(2)
%			for k = 1:SizeData(3)
%				
%				if min(CorDenomWeights(p,j,k,:)) < Info.Threshold
%					CorDenomWeights(p,j,k,:) = inf;
%				end
%			end
%		end
%	end
	CorDenomWeights(CorDenomWeights < Info.Threshold) = inf;
	Output = Mag./CorDenomWeights;
	toc
	disp('LFG correction is finished!')
end
