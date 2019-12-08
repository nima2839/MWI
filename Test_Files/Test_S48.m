%Mag =   ReadAllDCM('~/Nm_Myelin_Test_1/study/gre_96contrasts_new_10',48);
%Phase = ReadAllDCM('~/Nm_Myelin_Test_1/study/gre_96contrasts_new_11',48);

%Phase = double(Phase);
%Phase = 2*pi* Phase/4095 - pi;
%Mag = double(Mag) / 4095;

cd ~

%FEcho(:,:,:) = Mag(:,:,:,1);
%save_nii(FEcho,'Mag');
%unix('bet Mag.nii mask.nii -m');
%unix('rm mask.nii.gz');
%unix('gunzip -f mask_mask.nii.gz');
%unix('rm mask_mask.nii.gz')
%Mask = load_nii('mask_mask.nii');
%unix('rm mask_mask.nii');

%Mask = M;
%clear M

Info.Mask = Mask;
Info.FirstTE = 1.45e-3;
Info.EchoSpacing = 1.05e-3;
Info.Vox = [2 2 3] * 1e-3;
Info.EchoIndexes = [6 8];%testing
Info.Method = 1; %Du
[LFGC,Gp, Gv, Gs] = LFG_Correction(Mag, Phase, Info);


Opt.T_range = [1e-3 0.120];
Opt.Number_T = 120;
Opt.Vox = Info.Vox;

E.FirstEchoTime = 1.45e-3;
E.EchoSpacing = 1.05e-3;

%Fitted = NNLS_T_Batch_Fitting(Mag,E,Opt);
%LFGC_Fitted = NNLS_T_Batch_Fitting(Mag,E,Opt,
    X0 = [0.1 9e-3 10 0.85 6e-2 0.5 10];
    lb = [0 0.1e-3 -50 0 20e-3 8e-2 -50];
    ub = [1 20e-3 50 1 155e-3 5 50];

sd = size(Mag);
res = zeros(sd(1),sd(2),sd(3));
Params = zeros(sd(1),sd(2),sd(3),length(X0));

disp('3PM fitting started!...')
tic
np = sd(1);
nv = sd(2);
ns = sd(3);
Threshold = 0.01;
LastIndex = 60;
ResThreshold = 0.2;
Info.Algorithm = 'trust-region-reflective';
Info1 = Info;
Info1.Algorithm = 'levenberg-marquardt';
parfor i = 1:np
	for j = 1:nv
		for k = 1:ns
			if Mask(i,j,k) > 0
				tmp = squeeze(LFGC(i,j,k,:));
	                	tmp1 = tmp;
				tmp(tmp > Threshold*max(tmp(:))) = 0;
	                	a = find(tmp,1);
				if isempty(a)
					a = LastIndex;
				elseif a > LastIndex
					a = LastIndex;
				end
	                	if  a > length(X0)
	                   		tmpd = tmp1(1:a);
					
	                   		[Temp, r] = ThreePoolM_NLLS(tmpd,Info);
					% If the threshold is exceeded parameters will be calculated by Levenberg-Marquardt Algorithm
					if r > ResThreshold
						[Params(i,j,k,:) , res(i,j,k)] = ThreePoolM_NLLS(tmpd,Info1);
					else
						Params(i,j,k,:) = Temp;
						res(i,j,k) = r;
					end		
	                	end
			end
		end
	end
end
disp('Finished fitting 3PM!')
runtime = toc;

MWF = Find_MWF(Params, [1 4 7], 'NLLS_Test');
%LFGC_MWF = Find_MWF(LFGC_Fitted, 35,'NNLS');

disp('Saving results...')
Description = 'Calculating MWF from LFGC! 7Param 3PM! First 60Echos!';
clear Mag Phase
save('NM_Test1_3PM_BoB_60E')
disp('Done!')
