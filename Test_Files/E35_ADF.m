load NMT2_FS_data.mat

MyInfo.FirstTE = 1.45e-3;
MyInfo.EchoSpacing = 2* 1.05e-3;
MyInfo.Vox = [1 1 3]*1e-3;
MyInfo.Mask = Mask;

% mytest.EchoIndexes = [1, 96];

tic
test = TestClass(Mag(:,:,:,1:2:35),Phase(:,:,:,1:2:35),MyInfo);
test = CalcLFGC(test);
adf = ADF(test.LFGC,MyInfo.Vox,10,100);
filtered = ApplyADF(adf);
test = SetLFGC(test,filtered);
test = Calc_SC(test,2);
test = Calc_2PM(test);
test = Calc_3PM(test);
disp('Saving results...')
test.Description = 'Calculating MWF from LFGC! 8Param 3PM! First Echoes 1:2:35! ADF!';
test.RunTime = toc;
data = GetAllData(test);
save('T2_35_ADF','data');
disp('Done!')
