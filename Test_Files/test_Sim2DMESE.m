% This is to employ Sim2DMESE to simulate decay curves and apply Multi-component analysis

disp('Initializing simulation parameters...');
% First we define Distribution for each water compartment
IE = SimClass.Create_Guassian_Dist(75e-3); % intra/extra-cellular water 
MW = SimClass.Create_Guassian_Dist(15e-3); % myelin water

% Now define other sim parameters
MWFs = (1:2:30) *1e-2;
SNRs = [100,200, 300, 500, 750, 1e3];

% Get the sequence parameters
cd ~/MESE/
load('2D_32echos_MESE_Nima_dic', 'param')
MyInfo.SeqParams = param;
MyInfo.NumData = 500;
MyInfo.B1Range = 0.5:1e-2:1.2;

for i = 1:length(MWFs)
	temp.T2Values = [MW.T2Values, IE.T2Values];
	temp.Weights = [MWFs(i) * MW.Weights, (1 - MWFs(i)) * IE.Weights];
	MyInfo.T2Dist(i) = temp;
end

disp('Generating simulation object...');
tic;
SimObj = Sim2DMESE(MyInfo);
SimTime = toc;

disp('Applying multi-component analysis while iterating through all the SNRs...');

tic;
Maps = cell(1,length(SNRs));
Maps_B1 = cell(1,length(SNRs));
for i = 1:length(SNRs)
	[Maps{i}, Maps_B1{i}] = MC_Analyzer(SimObj, SNRs(i));
	disp(strcat(string(100*i/length(SNRs)),'%');
end
AnlysisTime = toc;

cd ~/Simulation/MESE2D/

save('MESE_2D_B1_suppliedANDestimated')