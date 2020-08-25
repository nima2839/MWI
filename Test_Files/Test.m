cd ~/GRE/GRE_To_Do/
load NMT3_2DGRE_18Cont_Monopolar.mat



sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.3);
complex_data = Mag.*exp(1i*Phase);
clear Mag Phase

for i = 1:sd(4)
	filtered(:,:,:,i) =  Info.Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end

tic
test = TestClass(NESMA_Filter(abs(filtered),Info.Mask, true,0.03),angle(filtered),Info);
test = CalcLFGC(test);


X0 = [0.1,   60,	  5,	0.6,	15,	0.3,	25,	   0];
test = Calc_3PM(test);


disp('Saving results...')
test.Description = 'Calculating 8Param 3PM! LFGC!Tukey alpha = 0.3';
RunTime = toc;

cd ~/GRE/GRE_Results/
save('18Cont_2DMonopolar_NESMA');
disp('Done!')
