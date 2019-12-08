load NMT2_FS_data.mat

MyInfo.FirstTE = 1.45e-3;
MyInfo.EchoSpacing = 2* 1.05e-3;
MyInfo.Vox = [1 1 3]*1e-3;
MyInfo.Mask = Mask;

K = Tukey3D(6,6,1,0.5);

complex_data = Mag.*exp(1i*Phase);

sd = size(Mag);

Mask = double(Mask);

clear Mag Phase
for i = 1:sd(4)
	filtered(:,:,:,i) =  convn(Mask.* complex_data(:,:,:,i),K,'same') ;
end

tic
test = TestClass(abs(filtered(:,:,:,1:2:35)),angle(filtered(:,:,:,1:2:35)),MyInfo);
test = CalcLFGC(test);
test = Calc_SC(test,2);
test = Calc_2PM(test);
test = Calc_3PM(test);
disp('Saving results...')
test.Description = 'Calculating MWF from LFGC! 8Param 3PM! First Echoes 1:2:35! Tukey3D(6,6,1,0.5)!';
test.RunTime = toc;
data = GetAllData(test);
save('T2_35_Tukey','data');
disp('Done!')
