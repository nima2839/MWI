cd ~/GRE/GRE_To_Do/
load Feb20_GRE_20cont_Monopolar.mat



sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.35);
complex_data = Mag.*exp(1i*Phase);
clear Mag Phase

for i = 1:sd(4)
	filtered(:,:,:,i) =  Info.Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end

tic
test = TestClass(abs(filtered),angle(filtered),Info);
test = CalcLFGC(test);
clear complex_data filtered
%for i = 1:sd(4)
%	adfiltered(:,:,:,i) = imdiffusefilt(test.LFGC(:,:,:,i),'NumberOfIterations',6);
%end
test = SetLFGC(test,adfiltered);
test = Calc_SC(test,2); % LOG method
test = Calc_2PM(test);
test = Calc_3PM(test);
%test.Description = 'Calculating 8Param C3PM! LFGC!Tukey alpha = 0.35!!';
data = GetAllData(test);
RunTime = toc;
clear test filtered
cd ~/GRE/GRE_Results/
save('20Cont_192_Monopolar_MEGRE', 'data');
disp('Done!')
