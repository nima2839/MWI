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
test = TestClass(NESMA_Filter(abs(filtered),Info.Mask, true,0.05),angle(filtered),Info);
test = CalcLFGC(test);
%for i = 1:sd(4)
%	adfiltered(:,:,:,i) = imdiffusefilt(test.LFGC(:,:,:,i),'NumberOfIterations',6);
%end
%test = SetLFGC(test,adfiltered);
%test = Calc_SC(test,2); % LOG method
%test = Calc_2PM(test);
%test = Calc_3PM(test);
X0 = [0.1,   300,	  5,	0.6,	15,	0.3,	25,	   0];
test1 = Calc_3PM(test);
MWF{1} = test1.MWF_3PM;
Res{1} = test1.Res_3PM;
clear test1
X0 = [0.1,   100,	  5,	0.6,	15,	0.3,	25,	   0];
test2 = Calc_3PM(test);
MWF{2} = test2.MWF_3PM;
Res{2} = test2.Res_3PM;
clear test2
X0 = [0.1,   80,	  5,	0.6,	15,	0.3,	25,	   0];
test3 = Calc_3PM(test);
MWF{3} = test3.MWF_3PM;
Res{3} = test3.Res_3PM;
clear test3
X0 = [0.1,   60,	  5,	0.6,	15,	0.3,	25,	   0];
test4 = Calc_3PM(test);
MWF{4} = test4.MWF_3PM;
Res{4} = test4.Res_3PM;

for i = 1:sd(1)
	for j = 1:sd(2)
		for k = 18%:sd(3)
			res = [Res{1}(i,j,k), Res{2}(i,j,k), Res{3}(i,j,k), Res{4}(i,j,k)];
			idx = find(res == min(res));
			Final_MWF(i,j,k) = MWF{idx(1)}(i,j,k);
		end
	end
end

disp('Saving results...')
test.Description = 'Calculating 8Param C3PM! LFGC!Tukey alpha = 0.25!4seed test!';
RunTime = toc;
clear test filtered
cd ~/GRE/GRE_Results/
save('20Cont_192_Monopolar_NESMA_Tukey_MultiSeed');
disp('Done!')
