load Gre2D_18Con_Monopolar.mat

myinfo.Mask = Mask;
myinfo.Method = 1;	% Du method
myinfo.FirstTE = 2.01e-3;
myinfo.EchoSpacing = 2.87e-3;
myinfo.Vox = [1.5 1.5 3] * 1e-3;

sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.5);
complex_data = Mag(:,:,:,:).*exp(1i*Phase(:,:,:,:));
clear Mag Phase

for i = 1:sd(4)
	filtered(:,:,:,i) = Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end
%filtered = complex_data;
tic
test = TestClass(abs(filtered),angle(filtered),myinfo);
test = CalcLFGC(test);
test = Calc_SC(test,2); % LOG method
%test = Calc_2PM(test);
%test = Calc_3PM(test);
test = Calc_Complex3PM(test);
disp('Saving results...')
test.Description = 'Calculating 8Param complex 3PM ! LFGC!!';
test.RunTime = toc;
data = GetAllData(test);
save('18Cont_Monopolar_Tukey_C3PM','data');
disp('Done!')
