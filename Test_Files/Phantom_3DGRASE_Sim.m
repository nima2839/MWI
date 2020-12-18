%addpath(genpath('~\MWI'))
tic;
disp('Creating 3D GRASE phantom...')
N = 120;
p = phantom3d('Modified Shepp-Logan', N);
p(p==0) =nan;
p(p == 0.3) = 0.15;
p(p == 0.2) = 0.05;
[X, Y, Z] = meshgrid(linspace(1,120, N), linspace(1,160,N), 1:N);
[Xq, Yq, Zq] = meshgrid(1:120, 1:160, 1:N);
p1 = interp3(X,Y,Z, p, Xq,Yq,Zq);
Phantom_3D = abs(p1(:,:, (floor(N/2) - 15): (floor(N/2) +16)));
%%
disp('Generating a flip angle map...');
opt.Size = size(Phantom_3D);
opt.Vox = [1.5, 1.5, 5];
opt.deltaFA = 0.5;
opt.CenterFA = 180;
FA_Map = Create_FA_Map(opt);
FA_Map(isnan(Phantom_3D)) = 0;

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
MyInfo.SNR = 200;

Phantom_4D_T2_1 =  Create_MWI_Phantom4D(Phantom_3D, FA_Map, MyInfo);
MyInfo.IE = SimClass.Create_Guassian_Dist(100e-3); % intra/extra-cellular water 
MyInfo.MW = SimClass.Create_Guassian_Dist(20e-3); % myelin water
Phantom_4D_T2_2 = Create_MWI_Phantom4D(Phantom_3D, FA_Map, MyInfo);
Phantom_4D_T2 = Phantom_4D_T2_2;
reshapedPhantom = reshape(Phantom_3D, [size(Phantom_3D,1)*size(Phantom_3D,2)*size(Phantom_3D,3),1]);
Phantom_4D_T2 =  reshape(Phantom_4D_T2,[size(reshapedPhantom,1), length(MyInfo.Times)]);
reshped_T2_4D_1 = reshape(Phantom_4D_T2_1, size(Phantom_4D_T2));
idx = find(reshapedPhantom > 0.1);
Phantom_4D_T2(idx,:) = reshped_T2_4D_1(idx,:); 
Phantom_4D_T2 = reshape(Phantom_4D_T2, size(Phantom_4D_T2_2));
clear Phantom_4D_T2_2 Phantom_4D_T2_1 reshped_T2_4D_1
%%
disp('Generating 4D phantom: T2* weighted...');
MyInfo.NumWaterComp = 2;
MyInfo.Times = (1:32)*1e-2;
MyInfo.IE = SimClass.Create_Guassian_Dist(40e-3); % intra/extra-cellular water 
MyInfo.MW = SimClass.Create_Guassian_Dist(5e-3); % myelin water
MyInfo.T2Dist.T2Values = [MyInfo.MW.T2Values, MyInfo.IE.T2Values];
MyInfo.T1Val = [.6*ones(size(MyInfo.MW.Weights)), ones(size(MyInfo.IE.Weights))];
MyInfo.FlipAngle = 180;
MyInfo.NumData = 1;
MyInfo.SNR = 200;

Phantom_4D_T2star_1 =  Create_MWI_Phantom4D(Phantom_3D, FA_Map, MyInfo);
MyInfo.IE = SimClass.Create_Guassian_Dist(30e-3); % intra/extra-cellular water 
MyInfo.MW = SimClass.Create_Guassian_Dist(10e-3); % myelin water
Phantom_4D_T2star_2 =  Create_MWI_Phantom4D(Phantom_3D, FA_Map, MyInfo);
Phantom_4D_T2star = Phantom_4D_T2star_2;
reshapedPhantom = reshape(Phantom_3D, [size(Phantom_3D,1)*size(Phantom_3D,2)*size(Phantom_3D,3),1]);
Phantom_4D_T2star =  reshape(Phantom_4D_T2star, [size(reshapedPhantom,1), length(MyInfo.Times)]);
reshped_T2star_4D_1 = reshape(Phantom_4D_T2star_1, size(Phantom_4D_T2star));
idx = find(reshapedPhantom > 0.1 | reshapedPhantom < 0.01);
Phantom_4D_T2star(idx,:) = reshped_T2star_4D_1(idx,:).*sinc(2*MyInfo.Times); 
Phantom_4D_T2star = reshape(Phantom_4D_T2star, size(Phantom_4D_T2star_2));
clear Phantom_4D_T2star_2 Phantom_4D_T2star_1 reshped_T2star_4D_1 reshapedPhantom
%%
disp('Generating3D GRASE sim phantom...');
sd = size(Phantom_4D_T2);
GRASE_Phantom = zeros(sd);

for e = 1:sd(4)
	Coil_effect = (sin(FA_Map*pi/360).^-12);
	Coil_effect(Coil_effect > 10) = 0;
	temp = fftshift(fftn(squeeze(Phantom_4D_T2(:,:,:,e).* Coil_effect)));
	tempstar = fftshift(fftn(squeeze(Phantom_4D_T2star(:,:,:,e).* Coil_effect)));
	temp(:, 1:floor(sd(2)/3), :) = tempstar(:, 1:floor(sd(2)/3), :);
	temp(:, end-floor(sd(2)/3):end, :) = tempstar(:, end-floor(sd(2)/3):end, :);
	GRASE_Phantom(:,:,:,e) = ifftn(ifftshift(temp));
end


Phantom_Sim_time = toc

disp('Applying MWI Analysis of GRASE_Phantom...');
 [maps,distributions,~] = T2map_Nima(abs(GRASE_Phantom), 'Threshold', 0.1, 'T2Range', [0.008, 2],'nT2', 60);


  MWI = squeeze(squeeze(sum(distributions(:,:,:,1:18),4))./squeeze(sum(distributions(:,:,:,:),4)));

cd ~/Simulation/Phantom/
save('Phantom_3DGRASE_T2starEffect','-v7.3')
cd ~/Simulation/Phantom/NIFTI_Files/
niftiwrite(180 - abs(180 - FA_Map), 'FA_Map')
niftiwrite(Phantom_3D, 'Phantom_3D')
niftiwrite(Phantom_4D_T2, 'P_4D_T2')
niftiwrite(Phantom_4D_T2star, 'P_4D_T2S')
niftiwrite(GRASE_Phantom, 'GRASE_Phantom')
niftiwrite(maps.alpha , 'alpha')
niftiwrite(MWI, 'MWF')