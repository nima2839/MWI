function ReProcessGRASE(name)
disp(['ReProcessGRASE on: ',name])
try
cd ~/GRASE/GRASE_To_Do/
load(name,'tf_mgrase')
tic

[maps,distributions,~] = T2map_SEcorr(tf_mgrase, 'Threshold', 200,'nT2', 60,'T2Range', [0.008, 2], 'MinRefAngle', 100);


MWI = squeeze(squeeze(sum(distributions(:,:,:,1:18),4))./squeeze(sum(distributions(:,:,:,:),4)));

MWI_1 = ( ones( size( MWI )) - isnan( MWI ) ) ;
[ Xres , Yres , Zres ] = size( MWI );
Final_MWI = zeros( size( MWI ) ) ;
for c = 1 : Zres
    for a = 1 : Xres
        for b = 1 : Yres
            if MWI_1( a , b , c ) == 1
                Final_MWI( a , b , c ) = MWI( a , b , c ) ;
            end
        end
    end
end
clear MWI MWI_1 tf_mgrase
runtime=toc;
Description = 'Threshold = 200; nT2 = 60, T2Range = 8ms to 2s, MinRefAngle = 100 degrees ';
cd ~/GRASE/GRASE_Results/GRASE_UBC_Results
save(name)
catch ME
	disp(ME.message)
	end
end
