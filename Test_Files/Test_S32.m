load NM_ExSp_2Avg.mat

Mask = double(Mask);
myinfo.Mask = Mask;
myinfo.Method = 1;	% Du method
myinfo.FirstTE = 2e-3;
myinfo.EchoSpacing = 2*1.05e-3;
myinfo.Vox = [2 2 4] * 1e-3;
K = Tukey3D(10,10,3,0.5);

disp('Applying Mask!...')
tic
sd = size(Mag);

test = TestClass(Mag(:,:,:,1:2:35),Phase(:,:,:,1:2:35),myinfo);
test = CalcLFGC(test);
sd = size(test.LFGC);
for i = 1:sd(4)
	filtered(:,:,:,i) =  convn(test.LFGC(:,:,:,i),K,'same') ;
end
test = SetLFGC(test,filtered);
test = Calc_SC(test,2); % LOG method
test = Calc_2PM(test);
test = Calc_3PM(test);
disp('Saving results...')
test.Description = 'Calculating 8Param 3PM and 5Param 2PM! LFGC! Echoes 1:2:35! Tukey windowed(10,10,3,0.5)!';
test.RunTime = toc;
data = GetAllData(test);
save('ExSp_2Avg_MWI_OE35_tukey','data');
disp('Done!')
