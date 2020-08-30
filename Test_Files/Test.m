cd ~/GRE/GRE_To_Do/
load NMT3_2DGRE_18Cont_Monopolar.mat



sd = size(Mag);
K = repmat(Tukey3D(sd(1),sd(2),1,0.25),[1,1,sd(3)]);
complex_data = Mag.*exp(1i*Phase);
clear Mag Phase

for i = 1:sd(4)
	filtered(:,:,:,i) =  Info.Mask.*ifftn(fftshift(fftshift(fftn((Mag_Bias.^-2).*complex_data(:,:,:,i))).*K)) ;
end

tic

test = TestClass(abs(filtered),angle(filtered),Info);
test = CalcLFGC(test);
idx = 15:17;
test.LFGC(:,:,idx,:) = NESMA_Filter(test.LFGC(:,:,idx,:),Info.Mask(:,:,idx), false, 0.05);


%test2 = Calc_Multi_Seed(test);
test1 = Calc_3PM(test);

disp('Saving results...')
test.Description = 'Calculating 8Param 3PM! LFGC!Tukey alpha = 0.3';
RunTime = toc;
MWF = test1.MWF_3PM;
%MWF_MS = test2.MWF_3PM;
cd ~/GRE/GRE_Results/
save('18Cont_2DMonopolar_NESMA','MWF');%, 'MWF_MS');
disp('Done!')
