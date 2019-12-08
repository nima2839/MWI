load combined_img.mat


Threshold = 200;
MyInfo.FirstTE = 1.45e-3;
MyInfo.EchoSpacing =  2*1.05e-3;
MyInfo.Vox = [1 1 3]*1e-3;
MyInfo.Mask = abs(combined_img{1}(:,:,:,1));
MyInfo.Mask(MyInfo.Mask < Threshold) = 0;

MyInfo2.FirstTE = 1.45e-3;
MyInfo2.EchoSpacing =  2*1.05e-3;
MyInfo2.Vox = [1 1 3]*1e-3;
MyInfo2.Mask = abs(combined_img{2}(:,:,:,1));
MyInfo2.Mask(MyInfo2.Mask < Threshold) = 0;

sd = size(combined_img{1});
K = Tukey3D(sd(1),sd(2),sd(3),0.5);

complex_data = combined_img{1};
complex_dat2 = combined_img{2};

clear combined_img

for i = 1:sd(4)
	filtered(:,:,:,i) =  ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
	filtered2(:,:,:,i) =  ifftn(fftshift(fftshift(fftn(complex_data2(:,:,:,i))).*K)) ;
end

tic
test = TestClass(abs(filtered(:,:,:,1:2:end)),angle(filtered(:,:,:,1:2:end)),MyInfo);
test = CalcLFGC(test);
test = Calc_SC(test,2);
test = Calc_2PM(test);
test = Calc_Complex3PM(test);
disp('Saving results...')
test.Description = 'Odd echoes from the first average!LFGC!Tukey3D!';
test.RunTime = toc;
data = GetAllData(test);

tic
test2 = TestClass(abs(filtered2(:,:,:,1:2:end)),angle(filtered2(:,:,:,1:2:end)),MyInfo);
test2 = CalcLFGC(test2);
test2 = Calc_SC(test2,2);
test2 = Calc_2PM(test);
test2 = Calc_Complex3PM(test2);
disp('Saving results...')
test2.Description = 'Odd echoes from the second average!LFGC!Tukey3D!';
test2.RunTime = toc;
data2 = GetAllData(test2);
save('GP_40_Tukey_OE','data','data2');
disp('Done!')
