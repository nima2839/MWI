load 7gre.mat

Mask = double(Mask);
myinfo.Mask = Mask;
myinfo.Method = 1;	% Du method
myinfo.FirstTE = 3.69e-3;
myinfo.EchoSpacing = 5.13e-3;
myinfo.Vox = [1 1 1.5] * 1e-3;

disp('Applying Mask!...')
tic
for i = 1:7
	maskedMag(:,:,:,i) = Mask.*Mag(:,:,:,i);
end
toc
test = TestClass(Mag(:,:,20:60,:),Mag(:,:,20:60,:),myinfo);
test = Calc_SC(test,2); % LOG method
test = Calc_2PM(test);
disp('Saving results...')
test.Description = 'Calculating MWF! 8Param 3PM! Silices 20:60, 7 Echos, no LFGC!';
test.RunTime = toc;
data = GetAllData(test);
save('7GRE_MWI','data');
disp('Done!')
