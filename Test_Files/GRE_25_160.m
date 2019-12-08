load NMT3_25cont_160_bipolar_10.mat

myinfo.Mask = Mask(:,:,:);
myinfo.Method = 1;	% Du method
myinfo.FirstTE = 2.02e-3;
myinfo.EchoSpacing = 2*1.87e-3;
myinfo.Vox = [1.5 1.5 3] * 1e-3;
myinfo.BipolarFlag = true;


complex_data = Mag(:,:,:,:).*exp(1i*Phase(:,:,:,:));
clear Mag Phase
%filtered =complex_data;

sd = size(complex_data);
K = Tukey3D(sd(1),sd(2),sd(3),0.4);
for i = 1:sd(4)
	filtered(:,:,:,i) =  ifftn(fftshift(fftshift(fftn(myinfo.Mask.* complex_data(:,:,:,i))).*K)) ;
end

tic
test = TestClass(abs(filtered),angle(filtered),myinfo);
%test = CalcLFGC(test);
test = Calc_SC(test,2); % LOG method
%test = Calc_2PM(test);
%test = Calc_3PM(test);
test = Calc_Complex3PM(test);
disp('Saving results...')
test.Description = 'Calculating C3PM ! LFGC!Tukey alpha = 0.5! bipolar correction + ES = 2*ES ';
test.RunTime = toc;
data = GetAllData(test);
save('25Cont_160_bipolar_C3PM_false','data');
disp('Done!')
