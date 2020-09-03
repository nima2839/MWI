function ReProcessGRASE(name)
disp(['ReProcessGRASE on: ',name])
try
cd ~/GRASE/GRASE_To_Do/
load(name,'tf_mgrase')
cd ~/GRASE/GRASE_Results/
load(strcat('GRASE_Results_',name), 'maps')
alpha = maps.alpha;
clear maps
alpha(isnan(alpha)) = 0;
K = Tukey3D(9,9,5, 1);
alpha = convn(alpha,K, 'same');
tic

[maps,distributions,~] = T2map_Nima(tf_mgrase, 'Threshold', 200,'FlipAngleMap', alpha);


Final_MWI = squeeze(squeeze(sum(distributions(:,:,:,1:18),4))./squeeze(sum(distributions(:,:,:,:),4)));

Final_MWI(isnan(Final_MWI)) = 0 ;


runtime=toc;
Description = 'Alpha map is convolved with Tukey3D(9,9,5, 1) -> Hanning window!';
cd ~/GRASE/GRASE_Results/Filtered_B1_Map_Results
save(strcat('Filtered_B1_Map_Results_',name))
catch ME
	disp(ME.message)
	end
end
