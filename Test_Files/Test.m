cd ~/GRE/GRE_To_Do/
load NMT3_2DGRE_18Cont_Monopolar.mat



sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.25);
complex_data = Mag.*exp(1i*Phase);
clear Mag Phase

for i = 1:sd(4)
	filtered(:,:,:,i) =  Info.Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end

tic

test = TestClass(abs(filtered),angle(filtered),Info);
test.Mag(:,:,17:22,:) = NESMA_Filter(test.Mag(:,:,17:22,:),Info.Mask(:,:,17:22),true, 0.02);
test = CalcLFGC(test);

%test1 = Calc_Multi_Seed(test);
test = Calc_3PM(test);

disp('Saving results...')
test.Description = 'Calculating 8Param 3PM! LFGC!Tukey alpha = 0.3';
RunTime = toc;
MWF = test.MWF_3PM;
cd ~/GRE/GRE_Results/
save('18Cont_2DMonopolar_NESMA','MWF');
disp('Done!')
