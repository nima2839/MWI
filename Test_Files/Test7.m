cd ~/GRE/GRE_To_Do/
load NMT3_2DGRE_18Cont_Monopolar.mat



sd = size(Mag);
K = repmat(Tukey3D(sd(1),sd(2),1,0.35),[1,1,sd(3)]);
complex_data = Mag.*exp(1i*Phase);
clear Mag Phase 

for i = 1:sd(4)
	filtered(:,:,:,i) =  Info.Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i)./Mag_Bias)).*K)) ;
end
clear complex_data Mag_Bias
tic

test = TestClass(abs(filtered),angle(filtered),Info);
test = CalcLFGC(test);
clear filtered

idx = 10:20;

opt.Mask = Info.Mask(:,:,idx);
opt.Num_Channels = 10;
opt.Method = "RMD"

test = SetLFGC(test, NESMA_Filter(test.LFGC(:,:,idx,:),opt));
test.MyInfo.Mask = Info.Mask(:,:,idx);
test = Calc_SC(test,2);
test = Calc_3PM(test);

disp('Saving results...')
test.Description = 'Calculating 8Param 3PM! LFGC!Tukey alpha = 0.35';
RunTime = toc;


cd ~/GRE/GRE_Results/
save('18Cont_2DMonopolar_RMD_NESMA');
disp('Done!')
