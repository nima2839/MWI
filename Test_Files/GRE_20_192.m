load NMT3_25cont_192_bipolar_12.mat

myinfo.Mask = Mask(:,:,20:22);
myinfo.Method = 1;	% Du method
myinfo.FirstTE = 1.84e-3;
myinfo.EchoSpacing = 2.62e-3;
myinfo.Vox = [1.3 1.3 2] * 1e-3;

sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.5);
complex_data = Mag(:,:,20:22,:).*exp(1i*Phase(:,:,20:22,:));
clear Mag Phase
filtered = complex_data;
%for i = 1:sd(4)
%	filtered(:,:,:,i) =  Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
%end

tic
test = TestClass(abs(filtered),angle(filtered),myinfo);
test = CalcLFGC(test);
for i = 1:sd(4)
	adfiltered(:,:,:,i) = imdiffusefilt(test.LFGC(:,:,:,i),'NumberOfIterations',6);
end
test = SetLFGC(test,adfiltered);
test = Calc_SC(test,2); % LOG method
%test = Calc_2PM(test);
%test = Calc_3PM(test);
test = Calc_Complex3PM(test);
disp('Saving results...')
test.Description = 'Calculating 8Param C3PM! LFGC!Tukey alpha = 0.5!ADF!';
test.RunTime = toc;
data = GetAllData(test);
save('20Cont_192_Monopolar_Tukey_ADF_C3PM_s21','data');
disp('Done!')
