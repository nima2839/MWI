function  ProcessGRASE( FileName )
% function: basic process of GRASE MWI given by UBC
tic;


cd ~/GRASE/GRASE_To_Do
load(FileName, 'mgrase')

fs = 1/3;
[ys xs zs es]=size(mgrase)
tf_mgrase = zeros(ys,xs,zs,es);
hfilt2=tukeywin(ys,fs)*tukeywin(xs,fs)';

for i = 1:es
    for j = 1:zs

        tf_mgrase(:,:,j,i) = abs(ifft2c(fft2c(mgrase(:,:,j,i)).*hfilt2));

    end
end

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

clear MWI MWI_1 mgrase hfilt2
runtime=toc;
cd ~/GRASE/GRASE_Results
Description = 'Threshold = 200; nT2 = 60, T2Range = 8ms to 2s, MinRefAngle = 100 degrees ';
save(['GRASE_Results_', FileName])

tic

[maps,distributions,~] = T2map_Nima(tf_mgrase, 'Threshold', 200,'nT2', 60,'T2Range', [0.008, 2], 'MinRefAngle', 100);


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

clear MWI MWI_1 
runtime=toc;
cd ~/GRASE/GRASE_Results
Description = 'Threshold = 200; nT2 = 60, T2Range = 8ms to 2s, MinRefAngle = 100 degrees ';
save(['GRASE_Results_Extended_B1_', FileName])
end
