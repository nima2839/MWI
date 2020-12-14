addpath(genpath('~\MWI'))
tic;
disp('Creating 3D GRASE phantom...')
N = 120;
p = phantom3d('Modified Shepp-Logan', N);
p(p==0) =nan;
[X, Y, Z] = meshgrid(linspace(1,120, N), linspace(1,160,N), 1:N);
[Xq, Yq, Zq] = meshgrid(1:120, 1:160, 1:N);
p1 = interp3(X,Y,Z, p, Xq,Yq,Zq);
Phantom_3D = abs(p1(:,:, (floor(N/2) - 15): (floor(N/2) +16)));
%%
disp('Generating a flip angle map...');
opt.Size = size(Phantom_3D);
opt.Vox = [1.5, 1.5, 5];
opt.deltaFA = 0.8;
FA_Map = Create_FA_Map(opt);
FA_Map(isnan(Phantom_3D)) = nan;

%%
disp('Generating 4D phantom: T2 weighted...');
MyInfo.NumWaterComp = 2;
MyInfo.Times = (1:32)*1e-2;
MyInfo.IE = SimClass.Create_Guassian_Dist(75e-3); % intra/extra-cellular water 
MyInfo.MW = SimClass.Create_Guassian_Dist(15e-3); % myelin water
MyInfo.T2Dist.T2Values = [MyInfo.MW.T2Values, MyInfo.IE.T2Values];
MyInfo.T1Val = [.6*ones(size(MyInfo.MW.Weights)), ones(size(MyInfo.IE.Weights))];
MyInfo.FlipAngle = 180;
MyInfo.NumData = 1;
MyInfo.SNR = 300;

Phantom_4D_T2 =  Create_MWI_Phantom4D(Phantom_3D, FA_Map, MyInfo);
%%
disp('Generating 4D phantom: T2* weighted...');
MyInfo.NumWaterComp = 2;
MyInfo.Times = (1:32)*1e-2;
MyInfo.IE = SimClass.Create_Guassian_Dist(55e-3); % intra/extra-cellular water 
MyInfo.MW = SimClass.Create_Guassian_Dist(10e-3); % myelin water
MyInfo.T2Dist.T2Values = [MyInfo.MW.T2Values, MyInfo.IE.T2Values];
MyInfo.T1Val = [.6*ones(size(MyInfo.MW.Weights)), ones(size(MyInfo.IE.Weights))];
MyInfo.FlipAngle = 180;
MyInfo.NumData = 1;
MyInfo.SNR = 300;

Phantom_4D_T2star =  Create_MWI_Phantom4D(Phantom_3D, FA_Map, MyInfo);

%%
disp('Generating3D GRASE sim phantom...');
sd = size(Phantom_4D_T2);
GRASE_Phantom = zeros(sd);
for z = 1:sd(3)
	for e = 1:sd(4)
		temp = fftshift(fft2(squeeze(Phantom_4D_T2(:,:,z,e))));
		tempstar = fftshift(fft2(squeeze(Phantom_4D_T2star(:,:,z,e))));
		temp(:, 1:floor(sd(2)/3)) = tempstar(:, 1:floor(sd(2)/3));
		temp(:, end-floor(sd(2)/3):end) = tempstar(:, end-floor(sd(2)/3):end);
		GRASE_Phantom(:,:,z,e) = ifft2(ifftshift(temp));
	end
end

Phantom_Sim_time = toc

disp('Applying MWI Analysis of GRASE_Phantom...');
 [maps,distributions,~] = T2map_Nima(GRASE_Phantom, 'Threshold', 200, 'T2Range', [0.008, 2],'nT2', 60);


  MWI = squeeze(squeeze(sum(distributions(:,:,:,1:18),4))./squeeze(sum(distributions(:,:,:,:),4)));



cd ~/Simulation/Phantom/
save('Phantom_3DGRASE_T2starEffect','-v7.3')