function out = Template_Model(MWF, Chi2Factor, TimeConstRange)
% to use a template for simulation results of SimClass

CSF = 0.1; % Fraction of CSF

IE = 1 - MWF - CSF; % Fraction of intra-extra cellular water

MyInfo.NumWaterComp = 3;
MyInfo.Times = (1:32)*1e-2;
if nargin < 3
	MyInfo.TimeConstRange{1} = [3 10]*1e-3;
	MyInfo.TimeConstRange{2} = [70 80]*1e-3;
	MyInfo.TimeConstRange{3} = [500 2000]*1e-3;
else
	MyInfo.TimeConstRange = TimeConstRange
end
MyInfo.T1Val = [.6 1 4.163];
MyInfo.FractionRange{1}= [MWF, MWF];
MyInfo.FractionRange{2}= [IE, IE];
MyInfo.FractionRange{3}= [CSF, CSF];



MyInfo.FlipAngle = 180;
MyInfo.NumData = 500;
MyInfo.TrueFAFlag = false;
MyInfo.SNR = 0;
%%
FA = 110:180;
SNR = [50,100:50:300, 500, 750, 1e3];
nFA = length(FA);
nSNR = length(SNR);
Dist =  cell(nSNR,nFA);
Maps = Dist;
TrueFA_Dist = Dist;
TrueFA_Maps = Dist;
parfor j = 1:nSNR
	temp = MyInfo;
	temp.SNR = SNR(j);
	temp_dist = cell(1,nFA);
	temp_maps = cell(1,nFA);
	temp_Tdist = cell(1,nFA);
	temp_Tmaps = cell(1,nFA);
	for i  = 1:nFA
		temp.FlipAngle = FA(i);
		temp.TrueFAFlag = false;
		a = SimClass(temp);
		[temp_dist{1,i}, temp_maps{1,i}] = UBC_Nima_Fitting(a);
		a.MyInfo.TrueFAFlag = true;
		[temp_Tdist{1,i}, temp_Tmaps{1,i}] = UBC_Nima_Fitting(a, ones(1,60), Chi2Factor);
	end
	Dist(j,:) = temp_dist;
	Maps(j,:) = temp_maps;
	TrueFA_Dist(j,:) = temp_Tdist;
	TrueFA_Maps(j,:)= temp_Tmaps;
end

out.MyInfo = MyInfo;
out.FA = FA;
out.SNR = SNR;
out.Dist = Dist;
out.Maps = Maps;
out.TrueFA_Dist = TrueFA_Dist;
out.TrueFA_Maps = TrueFA_Maps;
end