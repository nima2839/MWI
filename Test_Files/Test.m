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
test.Mag(:,:,17:19,:) = NESMA_Filter(test.Mag(:,:,17:19,:),Info.Mask(:,:,17:19),true, 0.01);
%test = CalcLFGC(test);

%test2 = Calc_Multi_Seed(test);
test1 = Calc_3PM(test);

disp('Saving results...')
test.Description = 'Calculating 8Param 3PM! LFGC!Tukey alpha = 0.3';
RunTime = toc;
MWF = test1.MWF_3PM;
%MWF_MS = test2.MWF_3PM;
cd ~/GRE/GRE_Results/
save('18Cont_2DMonopolar_NESMA','MWF');%, 'MWF_MS');
disp('Done!')
