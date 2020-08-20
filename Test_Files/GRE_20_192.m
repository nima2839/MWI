cd ~/GRE/GRE_To_Do/
load Feb20_GRE_20cont_Monopolar.mat



sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.25);
complex_data = Mag.*exp(1i*Phase);
clear Mag Phase

for i = 1:sd(4)
	filtered(:,:,:,i) =  Info.Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end

tic
test = TestClass(abs(filtered),angle(filtered),Info);
test = CalcLFGC(test);
%for i = 1:sd(4)
%	adfiltered(:,:,:,i) = imdiffusefilt(test.LFGC(:,:,:,i),'NumberOfIterations',6);
%end
%test = SetLFGC(test,adfiltered);
test = Calc_SC(test,2); % LOG method
%test = Calc_2PM(test);
%test = Calc_3PM(test);
X0 = [0.1,   300,	  5,	0.6,	15,	0.3,	25,	   0];
test1 = Calc_3PM(test);
X0 = [0.1,   100,	  5,	0.6,	15,	0.3,	25,	   0];
test2 = Calc_3PM(test);
X0 = [0.1,   80,	  5,	0.6,	15,	0.3,	25,	   0];
test3 = Calc_3PM(test);
X0 = [0.1,   60,	  5,	0.6,	15,	0.3,	25,	   0];
test4 = Calc_3PM(test);
disp('Saving results...')
test.Description = 'Calculating 8Param C3PM! LFGC!Tukey alpha = 0.25!4seed test!';
test.RunTime = toc;
cd ~/GRE/GRE_Results/
save('20Cont_192_Monopolar_Tukey_MultiSeed');
disp('Done!')
