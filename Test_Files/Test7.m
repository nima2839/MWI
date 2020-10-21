cd ~/GRE/GRE_To_Do/
load NMT3_2DGRE_18Cont_Monopolar.mat



sd = size(Mag);
K = repmat(Tukey3D(sd(1),sd(2),1,0.35),[1,1,sd(3)]);
complex_data = Mag.*exp(1i*Phase);
clear Mag Phase

for i = 1:sd(4)
	filtered(:,:,:,i) =  Info.Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K))./Mag_Bias ;
end

tic

test = TestClass(abs(filtered),angle(filtered),Info);
%test = CalcLFGC(test);


idx = 18:19;

opt.Mask = Info.Mask(:,:,idx);
opt.Method = "RMD"
temp =   NESMA_Filter(test.Mag(:,:,idx,:),opt);
test = SetLFGC(test,temp(:,:,1,:));

test = Calc_3PM(test);

disp('Saving results...')
test.Description = 'Calculating 8Param 3PM! LFGC!Tukey alpha = 0.35';
RunTime = toc;

MWF = test.MWF_3PM;
LFGC = temp;
Res = test.Res_3PM;

cd ~/GRE/GRE_Results/
save('18Cont_2DMonopolar_RMD_NESMA','MWF','LFGC','Res');
disp('Done!')
