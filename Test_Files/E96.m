%Mag = double(ReadAllDCM('/media/work/aelkady/3T/Ahmed_96Echo_3Dtest_2/study/gre_96contrasts_3avg_7',48));
%Phase = double(ReadAllDCM('/media/work/aelkady/3T/Ahmed_96Echo_3Dtest_2/study/gre_96contrasts_3avg_8',48));
%Phase = 2*pi*Phase/4095;
%Phase = Phase - pi;

Threshold = 100;
MyInfo.FirstTE = 1.45e-3;
MyInfo.EchoSpacing =  2*1.05e-3;
MyInfo.Vox = [2 2 3]*1e-3;
MyInfo.Mask = abs(Mag(:,:,23:25,1));
MyInfo.Mask(MyInfo.Mask > Threshold) = 1;


sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.5);

complex_data = Mag.*exp(1i*Phase);


clear Mag Phase

for i = 1:sd(4)
	filtered(:,:,:,i) =  ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end

tic
test = TestClass(abs(filtered(:,:,23:25,1:2:60)),angle(filtered(:,:,23:25,1:2:60)),MyInfo);
test = CalcLFGC(test);
test = Calc_SC(test,2);
%test = Calc_2PM(test);
test = Calc_3PM(test);
disp('Saving results...')
test.Description = 'Odd echoes from the first average!LFGC!Tukey3D!3PM!';
test.RunTime = toc;
data = GetAllData(test);

save('GP_60_Tukey_OE','data');
disp('Done!')
