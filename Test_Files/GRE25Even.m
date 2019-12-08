load NMT3_25cont_160_bipolar_10.mat

myinfo.Mask = Mask(:,:,20:22);
myinfo.Method = 1;	% Du method
myinfo.FirstTE = 2.02e-3+1.87e-3;
myinfo.EchoSpacing = 2*1.87e-3;
myinfo.Vox = [1.5 1.5 3] * 1e-3;

sd = size(Mag(:,:,20:22,2:2:end));
K = Tukey3D(sd(1),sd(2),sd(3),0.5);
complex_data = Mag(:,:,20:22,2:2:end).*exp(1i*Phase(:,:,20:22,2:2:end));
clear Mag Phase
filtered=complex_data;
%for i = 1:sd(4)
%	filtered(:,:,:,i) =  ifftn(fftshift(fftshift(fftn(myinfo.Mask.* complex_data(:,:,:,i))).*K)) ;
%end

tic
test = TestClass(abs(filtered),angle(filtered),myinfo);
test = CalcLFGC(test);
test = Calc_SC(test,2); % LOG method
%test = Calc_2PM(test);
%test = Calc_3PM(test);
test = Calc_Complex3PM(test);
disp('Saving results...')
test.Description = 'Calculating C3PM ! LFGC!Tukey alpha = 0.5!Echoes 2:end';
test.RunTime = toc;
data = GetAllData(test);
save('25Cont_160_bipolar_Tukey_C3PM_s21_Even','data');
disp('Done!')
