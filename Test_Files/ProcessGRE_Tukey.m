function out  = ProcessGRE_Tukey(FileName,ReadPath,SavePath)

cd(ReadPath)
disp(['Reading ', FileName,]);
load(FileName)

%Threshold = 100;
MyInfo.FirstTE = 2.05e-3;
MyInfo.EchoSpacing =  1.01e-3;
MyInfo.Vox = [1 1 3]*1e-3;
MyInfo.Mask = Mask;
%MyInfo.Mask(MyInfo.Mask > Threshold) = 1;


sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.5);

complex_data = Mag.*exp(1i*Phase);


clear Mag Phase
disp('Applying Tukey3D filter!');
for i = 1:sd(4)
	filtered(:,:,:,i) =  ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end

tic
test = TestClass(abs(filtered),angle(filtered),MyInfo);
test = CalcLFGC(test);
test = Calc_SC(test,2);
%test = Calc_2PM(test);
test = Calc_Complex3PM(test);
disp('Saving results...')
test.Description = 'LFGC!Tukey3D!';
test.RunTime = toc;
data = GetAllData(test);

cd(SavePath)

save(['C3PM_Results' , FileName]), 'data');
end
