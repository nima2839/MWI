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
test.Mag(:,:,16,:) = NESMA_Filter(reshape(test.Mag(:,:,16,:), [sd(1:2),1,sd(4)]),Info.Mask(:,:,16),true, 0.05);
%test = CalcLFGC(test);

test = Calc_Multi_Seed(test);


disp('Saving results...')
test.Description = 'Calculating 8Param 3PM! LFGC!Tukey alpha = 0.3';
RunTime = toc;

cd ~/GRE/GRE_Results/
save('18Cont_2DMonopolar_MultiSeedTest');
disp('Done!')
