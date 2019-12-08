load NMT3_32cont_192_bipolar_8.mat

Mask = double(Mask);
myinfo.Mask = Mask;
myinfo.Method = 1;	% Du method
myinfo.FirstTE = 1.88e-3;
myinfo.EchoSpacing = 1.32e-3;
myinfo.Vox = [1.3 1.3 2] * 1e-3;
Kernel = ones(5,5,3);

disp('Applying Mask!...')
tic
sd = size(Mag);
K = Hann3D(5,5,3);
complex_data = Mag.*exp(1i*Phase);
for i = 1:sd(4)
	maskedMag(:,:,:,i) =  convn(Mask.* complex_data(:,:,:,i),K,'same') ;
end
toc
test = TestClass(abs(maskedMag),angle(maskedMag),myinfo);
test = Calc_NNLS(test);
disp('Saving results...')
test.Description = 'Calculating NNLS! LFGC!  Hanning windowed (5,5,3)';
test.RunTime = toc;
data = GetAllData(test);
save('32cont_192_bipolar_8_NNLS','data');
disp('Done!')
