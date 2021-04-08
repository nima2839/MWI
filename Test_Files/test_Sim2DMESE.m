% This is to employ Sim2DMESE to simulate decay curves and apply Multi-component analysis

disp('Initializing simulation parameters...');
% First we define Distribution for each water compartment
IE = SimClass.Create_Guassian_Dist(60e-3); % intra/extra-cellular water 
MW = SimClass.Create_Guassian_Dist(10e-3); % myelin water

% Now define other sim parameters
MWFs = 15 *1e-2;
SNRs = [200,  1e4];

% Get the sequence parameters
cd ~/MESE/
load('LookUpTable_101p')
%%
orig_info = MyInfo;
%MyInfo.SeqParams = Ryan.MyInfo.SeqParams;
%MyInfo.SeqParams.dp = linspace(-1.5,1.5,2001);
MyInfo.SeqParams.etl = 32;
MyInfo.NumData = 500;
MyInfo.B1Range = 0.6:.1:1.4;

MyInfo = rmfield(MyInfo, "LookUpTable");
B1_diff = -0.5:.05:.5;

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
%Maps = cell(1,length(SNRs));
Maps_B1 = cell(1,length(SNRs));
for i = 1:length(SNRs)
	[~, Maps_B1{i}] = MC_Analyzer(SimObj, SNRs(i),B1_diff);
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
save('Sim_2DMESE_B1diff', 'Maps' , 'Maps_B1')