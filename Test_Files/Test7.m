cd ~/GRE/GRE_To_Do/
load NMT3_2DGRE_18Cont_Monopolar.mat



sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.35);
complex_data = Mag.*exp(1i*Phase);
clear Mag Phase

myinfo = Info;
myinfo.Mask = Info.Mask(:,:,14:16);

disp('Applying Tukey filter!')
for i = 1:sd(4)
	filtered(:,:,:,i) =  Info.Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end

disp('Process started!')
tic
test = TestClass(abs(filtered(:,:,14:16,:)),angle(filtered(:,:,14:16,:)),myinfo);
test = CalcLFGC(test);
test.LFGC = NESMA_Filter(test.LFGC,Info.Mask(:,:,14:16),true, 0.02)

X0 = [0.1,   100,	  5,	0.6,	15,	0.3,	25,	   0];
test = Calc_3PM(test,X0);


disp('Saving results...')
test.Description = 'Calculating 8Param 3PM! LFGC!Tukey alpha = 0.3! NESMA without Normalizing!';
RunTime = toc;

cd ~/GRE/GRE_Results/
save('18Cont_2DMonopolar_NESMA_nonNormolized');
disp('Done!')
