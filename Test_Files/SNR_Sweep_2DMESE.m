% This is for an SNR sweep of 2D MESE simulations

disp('Initializing simulation parameters...');
% First we define Distribution for each water compartment
IE = SimClass.Create_Guassian_Dist(75e-3); % intra/extra-cellular water 
MW = SimClass.Create_Guassian_Dist(15e-3); % myelin water

% Now define other sim parameters
MWFs = 0.15;
SNRs = 50:50:1e3;

% Get the sequence parameters
cd ~/MESE/
load('LookUpTable_z257x129')
%%
orig_info = MyInfo;
%MyInfo.SeqParams = Ryan.MyInfo.SeqParams;
%MyInfo.SeqParams.dp = linspace(-1.5,1.5,2001);
MyInfo.SeqParams.etl = 32;
MyInfo.SeqParams.dx = linspace(-1.5,1.5, 129);
MyInfo.SeqParams.dz = linspace(-1.5,1.5, 257);
MyInfo.SeqParams.DSF = 10;
MyInfo.NumData = 1000;
MyInfo.B1Range = [0.7:0.05:1.3];

MyInfo = rmfield(MyInfo, "LookUpTable");
%B1_diff = -0.3:.01:.3;

for i = 1:length(MWFs)
	temp.T2Values = [MW.T2Values, IE.T2Values];
	temp.Weights = [MWFs(i) * MW.Weights, (1 - MWFs(i)) * IE.Weights];
	temp.T1Range = [1 * ones(size(MW.Weights)), ones(size(1 * IE.Weights))]; % setting T1 parameter
	MyInfo.T2Dist(i) = temp;
end

disp('Generating simulation object...');
tic;
SimObj = Sim2DMESE(MyInfo, orig_info.LookUpTable);
SimTime = toc;

disp('Applying multi-component analysis while iterating through all the SNRs...');

tic;
Maps = cell(1,length(SNRs));
Maps_B1 = cell(1,length(SNRs));
for i = 1:length(SNRs)
	[Maps{i}, Maps_B1{i}] = MC_Analyzer(SimObj, SNRs(i));%,B1_diff);
	disp(strcat(string(100*i/length(SNRs)),'%'));
end
AnlysisTime = toc;
%%
cd ~/Simulation/MESE2D/

%save('MESE_2D_B1_suppliedANDestimated','-v7.3')

for  i = 1:numel(SNRs)
	Maps{i}.Residuals = [];
	Maps{i}.Distribution = [];
	Maps_B1{i}.Residuals = [];
	Maps_B1{i}.Distribution = [];
end
save('Sim_2DMESE_SNR_Sweep_Alt_T1')