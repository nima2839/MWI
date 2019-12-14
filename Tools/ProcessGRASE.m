function  ProcessGRASE( FileName )
% function: basic process of GRASE MWI given by UBC
tic;

% Nima : testing an arbitrary T1 map for GRASE data
T1 = ones(1,200);
T1(1:25) = .250; % 250 ms
T1(25:40) = .5; % 500 ms
% The rest can be left at 1 second since the difference is negligible!

cd ~/GRASE/GRASE_To_Do
load(FileName)

fs = 1/3;
[ys, xs, zs, es] = size(mgrase);

tf_mgrase = zeros(ys,xs,zs,es);
hfilt2=tukeywin(ys,fs)*tukeywin(xs,fs)';

for i = 1:es
    for j = 1:zs

        tf_mgrase(:,:,j,i) = abs(ifft2c(fft2c(mgrase(:,:,j,i)).*hfilt2));

    end
end

clear mgrase

[maps,distributions,~] = T2map_Nima(tf_mgrase,'T1',T1, 'Threshold', 200, 'MinRefAngle', 60, 'nAngles', 12, 'T2Range', [0.01, 2]);


MWI = sqz(squeeze(sum(distributions(:,:,:,1:45),4))./squeeze(sum(distributions(:,:,:,:),4)));

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
Deccription = 'Threshold = 200; MinRefAngle= 60, nAngles = 10, T2Range = 1e-2 to 2, and T1 = 0.25,0.5,1, Cutoff = index 45';
save(['GRASE_Results_Special_Input_', FileName])
end
