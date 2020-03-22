function ReProcessGRASE(name)
load(name,'tf_mgrase')
tic

[maps,distributions,~] = T2map_SEcorr(tf_mgrase, 'Threshold', 200,'nT2', 60);


MWI = squeeze(squeeze(sum(distributions(:,:,:,1:13),4))./squeeze(sum(distributions(:,:,:,:),4)));

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
Description = 'Threshold = 200; nT2 = 60 ';
save(name)
end
