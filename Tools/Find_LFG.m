function [Gp, Gv, Gs] = Find_LFG(Echo1, Echo2, Info)
% Creates LFG maps 
%   Author: Nima 
%   Date: 02/2019
% Inputs:
%	Echo1/2:  Complex signal of each echo time ("Du" method doesn't need phase unwrapping)
%	Info (strutct):	
%	- Method: an integer indicating the method to use:
%		1. "Du" -> source: 'Multi-echo acquistion of MR Angiography..." by Du and colleagues
%		2. "Hwang" -> source: 'In vivo multi-silice mapping of myelin ..." by Hwang and colleagues
%	- deltaTE (s)
%	- Vox = [dx dy dz] (m)
%	- 

% Code definitions
DuMethod = 1;
HwangMethod = 2;
tempS =  size(Echo1);
if length(tempS) < 3
    error('Data must be 3D!');
end
SizeData = tempS(1:3);
np = SizeData(1);
nv = SizeData(2);
ns = SizeData(3);
Gp =  zeros(SizeData);
Gv = Gp;
Gs = Gv;

% Laying code constraints on inputs
%	if nargin < 3
%        error('Function requires 3 inputs!');
%    end
	
	if isfield(Info,'deltaTE') == 0
		error('Info.deltaTE can not be undefined!');
	end
	
	if isfield(Info,'Vox') == 0
		error('Info.Vox can not be undefined!');
	end
	
	if isfield(Info,'Method') == 0
        Info.Method = DuMethod;
    end
disp('Starting LFG acquiring process...')	
tic
	
	if Info.Method == DuMethod
		Denom = Info.Vox * 2*pi*42.575e6 * Info.deltaTE * 2;
		Wp = zeros(SizeData);
		Xp = Wp;Yp = Wp;Zp = Wp;
		Wv = Wp;Xv = Wp;Yv = Wp;Zv = Wp;
		Ws = Wp;Xs = Wp;Ys = Wp;Zs = Wp;
		
		% P-direction
		Wp(1:np-1,:,:) = Echo2(2:np,:,:)./abs(Echo2(2:np,:,:));
		Xp(1:np-1,:,:) = conj(Echo1(2:np,:,:))./abs(Echo1(2:np,:,:));
		Yp(2:np,:,:) = conj(Echo2(1:np-1,:,:))./abs(Echo2(1:np-1,:,:));
		Zp(2:np,:,:) = Echo1(1:np-1,:,:)./abs(Echo1(1:np-1,:,:));
		LFGp = angle(Wp.*Xp.*Yp.*Zp)/(Denom(1));
		% V-direction
		Wv(:,1:nv-1,:) = Echo2(:,2:nv,:)./abs(Echo2(:,2:nv,:));
		Xv(:,1:nv-1,:) = conj(Echo1(:,2:nv,:))./abs(Echo1(:,2:nv,:));
		Yv(:,2:nv,:) = conj(Echo2(:,1:nv-1,:))./abs(Echo2(:,1:nv-1,:));
		Zv(:,2:nv,:) = Echo1(:,1:nv-1,:)./abs(Echo1(:,1:nv-1,:));
		LFGv = angle(Wv.*Xv.*Yv.*Zv)/(Denom(2));
		% S-direction-direction
		Ws(:,:,1:ns-1) = Echo2(:,:,2:ns)./abs(Echo2(:,:,2:ns));
		Xs(:,:,1:ns-1) = conj(Echo1(:,:,2:ns))./abs(Echo1(:,:,2:ns));
		Ys(:,:,2:ns) = conj(Echo2(:,:,1:ns-1))./abs(Echo2(:,:,1:ns-1));
		Zs(:,:,2:ns) = Echo1(:,:,1:ns-1)./abs(Echo1(:,:,1:ns-1));
		LFGs = angle(Ws.*Xs.*Ys.*Zs)/(Denom(3));
		
		% Applying 5x5x3 median filter to reduce the effect of rapid phase change that may occur in vessels and other fine structures
		Kernel = [5,5,3];
		Gp = medfilt3(LFGp, Kernel);
		Gv = medfilt3(LFGv, Kernel);
		Gs = medfilt3(LFGs, Kernel);
        
	elseif Info.Method == HwangMethod
        	% This needs to be checked!!!
        	Denom = Info.Vox * 2*pi*42.575e6 * Info.deltaTE;
		Wp = zeros(SizeData);
		Yp = Wp;Zp = Wp;
		Wv = Wp;Xv = Wp;Yv = Wp;Zv = Wp;
		Ws = Wp;Xs = Wp;Ys = Wp;Zs = Wp;
		
		% P-direction
		Wp = Echo1;
		Xp = conj(Echo2);
		Yp(2:np,:,:) = conj(Echo1(1:np-1,:,:));
		Zp(2:np,:,:) = Echo2(1:np-1,:,:);
		LFGp = angle(Wp.*Xp.*Yp.*Zp)/(Denom(1));
		% V-direction
		Wv(:,1:nv-1,:) = Echo1(:,1:nv-1,:);
		Xv(:,1:nv-1,:) = conj(Echo2(:,1:nv-1,:));
		Yv(:,2:nv,:) = conj(Echo1(:,1:nv-1,:));
		Zv(:,2:nv,:) = Echo2(:,1:nv-1,:);
		LFGv = angle(Wv.*Xv.*Yv.*Zv)/(Denom(2));
		% S-direction-direction
		Ws(:,:,1:ns-1) = Echo1(:,:,1:ns-1);
		Xs(:,:,1:ns-1) = conj(Echo2(:,:,1:ns-1));
		Ys(:,:,2:ns) = conj(Echo1(:,:,1:ns-1));
		Zs(:,:,2:ns) = Echo2(:,:,1:ns-1);
		LFGs = angle(Ws.*Xs.*Ys.*Zs)/(Denom(3));

        	% Applying 5x5x3 median filter to reduce the effect of rapid phase change that may occur in vessels and other fine structures
		Kernel = [5,5,3];
		Gp = medfilt3(LFGp, Kernel);
		Gv = medfilt3(LFGv, Kernel);
		Gs = medfilt3(LFGs, Kernel);	
	end
toc
disp('Calculation of LFG is finished!')

end
