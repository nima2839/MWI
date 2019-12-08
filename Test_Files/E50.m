load VSFC_Mag_Phase.mat

myinfo.Mask = Mask;
myinfo.Method = 1;	% Du method
myinfo.FirstTE = 1.45e-3;
myinfo.EchoSpacing = 2*1.05e-3;
myinfo.Vox = [1 1 3] * 1e-3;


tic
%Threshold = 0.01;
%ResThreshold = 0.2;
%complex_data = Mag(:,:,:,1:2:59).*exp(1i*Phase(:,:,:,1:2:59));
%clear Mag Phase
%vsfc = VSF(abs(complex_data),angle(complex_data),myinfo);
%ffunc = Calc_F_Func(vsfc);
%complex_data = complex_data.* (ffunc.^-1);
%Mag = abs(complex_data);
%Phase = angle(complex_data);
%filt = ADF(Mag,[1 1 3],15,3e-4);
%fMag = ApplyADF(filt);
%myinfo.Mask = Mask(:,:,3:45);
test = TestClass(Mag(:,:,:,1:15),Phase(:,:,:,1:15),myinfo);
test = Calc_SC(test,2); % LOG method
test = Calc_3PM(test);
disp('Saving results...')
test.Description = 'Calculating MWF from LFGC! 8Param 3PM! Silices 3:45, 1:2:29 Echos, ADF was used on the VSF corrected data!';
test.RunTime = toc;
data = GetAllData(test);
save('NM_Test2_FS_30ms_VSF','data');
disp('Done!')
