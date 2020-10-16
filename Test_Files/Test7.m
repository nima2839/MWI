cd ~/GRE/GRE_To_Do/
load NESMA_Data.mat



sd = size(Mag);
K = Tukey3D(sd(1),sd(2),1,0.35);
K = repmat(K, [1,1,sd(3)]);
complex_data = Mag.*exp(1i*Phase);
clear Mag Phase

myinfo = Info;
slices = 10:20;
myinfo.Mask = Info.Mask(:,:,slices);

disp('Applying Tukey filter!')
for i = 1:sd(4)
	filtered(:,:,:,i) =  Info.Mask.*ifftn(fftshift(fftshift(fftn(complex_data(:,:,:,i))).*K)) ;
end

disp('Process started!')
tic
test = TestClass(abs(filtered(:,:,slices,:)),angle(filtered(:,:,slices,:)),myinfo);
test = CalcLFGC(test);
test.LFGC = NESMA_Filter(test.LFGC,myinfo.Mask,false, 0.05)

test = Calc_3PM(test);


disp('Saving results...')
test.Description = 'Calculating 8Param 3PM! LFGC!Tukey alpha = 0.3! NESMA without Normalizing!';
RunTime = toc;

cd ~/GRE/GRE_Results/
save('18Cont_2DMonopolar_NESMA_nonNormolized');
disp('Done!')
