function  ProcessGRASE( FileName )
% function: basic process of GRASE MWI given by UBC
tic;


cd ~/GRASE/GRASE_To_Do
load(FileName, 'tf_mgrase')

[maps,distributions,~] = T2map_Nima(tf_mgrase, 'Threshold', 200, 'MinRefAngle', 60, 'nAngles', 12, 'T2Range', [0.015, 2]);


MWI = sqz(squeeze(sum(distributions(:,:,:,1:40),4))./squeeze(sum(distributions(:,:,:,:),4)));

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

runtime=toc;
cd ~/GRASE/GRASE_Results
Deccription = 'Threshold = 200; MinRefAngle= 60, nAngles = 10, T2Range = 1e-2 to 2, and T1 = 0.2,0.5,1, Cutoff = index 40';
save(['GRASE_Results_Extended_B1_', FileName])
end
