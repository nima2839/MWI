% This is to employ Sim2DMESE to simulate decay curves and apply Multi-component analysis

disp('Initializing simulation parameters...');
% First we define Distribution for each water compartment
IE = SimClass.Create_Guassian_Dist(75e-3); % intra/extra-cellular water 
MW = SimClass.Create_Guassian_Dist(15e-3); % myelin water

% Now define other sim parameters
MWFs = [0:0.25:30] *1e-2;
SNRs = [100,200,500, 1e3];

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
MyInfo.NumData = 500;
MyInfo.B1Range = 0.5:.05:1.5;

MyInfo = rmfield(MyInfo, "LookUpTable");
%B1_diff = -0.3:.01:.3;

for i = 1:length(MWFs)
	temp.T2Values = [MW.T2Values, IE.T2Values];
	temp.Weights = [MWFs(i) * MW.Weights, (1 - MWFs(i)) * IE.Weights];
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

%for  i = 1:numel(SNRs)
%	Maps{i}.Residuals = [];
%	Maps{i}.Distribution = [];
%	Maps_B1{i}.Residuals = [];
%	Maps_B1{i}.Distribution = [];
%end
save('Sim_2DMESE_B1diffEffect')