load NMT3_32cont_192_bipolar_8.mat

myinfo.Mask = Mask;
myinfo.Method = 1;	% Du method
myinfo.FirstTE = 1.88e-3;
myinfo.EchoSpacing = 2*1.32e-3;
myinfo.Vox = [1.3 1.3 2] * 1e-3;


complex_data = Mag(:,:,:,1:2:end).*exp(1i*Phase(:,:,:,1:2:end));
clear Mag Phase
sd = size(complex_data);
K = Tukey3D(sd(1),sd(2),sd(3),0.5);

for i = 1:sd(4)
	filtered(:,:,:,i) =  ifftn(fftshift(fftshift(fftn(Mask.* complex_data(:,:,:,i))).*K)) ;
end

tic
test = TestClass(abs(filtered),angle(filtered),myinfo);
test = CalcLFGC(test);
test = Calc_SC(test,2); % LOG method
%test = Calc_2PM(test);
%test = Calc_3PM(test);
test = Calc_Complex3PM(test);
disp('Saving results...')
test.Description = 'Calculating C3PM! LFGC!Tukey alpha = 0.5!Odd echoes`';
test.RunTime = toc;
data = GetAllData(test);
save('32Cont_192_bipolar_Tukey_C3PM_OddE','data');
disp('Done!')
