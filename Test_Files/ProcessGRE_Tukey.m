function ProcessGRE_Tukey(FileName,ReadPath,SavePath)

cd(ReadPath)
disp(['Reading ', FileName,]);
load(FileName)

sd = size(Mag);
K = Tukey3D(sd(1),sd(2),sd(3),0.35);

complex_data = Mag.*exp(1i*Phase);


clear Mag Phase
disp('Applying Tukey3D filter!');
filtered = zeros(sd);
for i = 1:sd(4)
	filtered(:,:,:,i) =  ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end

tic
test = TestClass(abs(filtered),angle(filtered),Info);
clear filtered complex_data
test = Normal_Procedure(test);
disp('Saving results...')
test.RunTime = toc;
GRE_MWF_data = GetAllData(test);
clear test
cd(SavePath)

save(['C3PM_Results' , FileName, '_Tukey3D'],GRE_MWF_data);
end
