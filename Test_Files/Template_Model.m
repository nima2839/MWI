function out = Template_Model(MyInfo)
% to use a template for simulation results of SimClass
% Fully costumized!
% this shoulde be renamed to Parameter sweep(?)

% costumized to its core!
L = 20; % this is from SimClass.Create_Guassian_Dist()
T1Val = [];
for i = 1:MyInfo.NumWaterComp
	T1Val = [T1Val, MyInfo.T1Val(1) * ones(1,L)];
end
%%
FA = 110:180;
SNR = [100,200, 300, 500, 750, 1e3];
nFA = length(FA);
nSNR = length(SNR);
Dist =  cell(nSNR,nFA);
Maps = Dist;
TrueFA_Dist = Dist;
TrueFA_Maps = Dist;

MyInfo.TrueFAFlag = false;

parfor j = 1:nSNR
	Noise = Create_Noise(length(MyInfo.Times), MyInfo.NumData, SNR(j));
	
	temp = MyInfo;
	%temp.SNR = SNR(j);
	%temp_dist = cell(1,nFA);
	temp_maps = cell(1,nFA);
	%temp_Tdist = cell(1,nFA);
	temp_Tmaps = cell(1,nFA);
	for i  = 1:nFA
		temp.FlipAngle = FA(i);
		temp.TrueFAFlag = false;
		temp.SNR = 0; % to add costumized noise
		a = SimClass(temp);
		SimulatedData = a.SimulatedData + Noise;
		%temp_dist{1,i}
		[~, temp_maps{1,i}] = SimClass.UBC_Nima_Fitting(SimulatedData, temp);
		temp.TrueFAFlag = true;
		%temp_Tdist{1,i}
		[~, temp_Tmaps{1,i}] = SimClass.UBC_Nima_Fitting(SimulatedData, temp);
	end
	%Dist(j,:) = temp_dist;
	Maps(j,:) = temp_maps;
	%TrueFA_Dist(j,:) = temp_Tdist;
	TrueFA_Maps(j,:)= temp_Tmaps;
end


out.MyInfo = MyInfo;
out.FA = FA;
out.SNR = SNR;

%out.Dist = Dist; % to save storage space!
out.Maps = Maps;
%out.TrueFA_Dist = TrueFA_Dist; 
out.TrueFA_Maps = TrueFA_Maps;
end

function out = Create_Noise(ETL, NumData, SNR)
	out = zeros(NumData, ETL);
	for i = 1:NumData
		out(i,:) = SimClass.ADD_Noise(zeros(1,32), SNR, 1);
	end
end

%function out = Create_Distributions(MyInfo)
%	% returns an array of structs -> T2Dist
%	for i = 1:NumData
%		T2Dist.T2Values = [];
%		for k = 1:MyInfo.NumWaterComp
%			temp = 
%		end
%		out(i) = T2Dist;
%	end
%end