addpath(genpath('~/GRASE/Postprocessing'))
addpath(genpath('~/MWI'))

cd ~/GRASE/GRASE_To_Do

FileName = 'GRASE_Results_Can-05-Ppm-034-M0.mat';

load(FileName, 'tf_mgrase')

cd('~')
FlipAngleMap = 180 * double(niftiread('rB1_Phase.nii')) / (800*1.165);
FlipAngleMap = flip(permute(FlipAngleMap,[2 1 3]),1);
tic
[maps,distributions,~] = T2map_Nima(tf_mgrase, 'Threshold', 200, 'MinRefAngle', 60, 'nAngles', 12, 'T2Range', [0.015, 2], 'FlipAngleMap', FlipAngleMap);


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
Deccription = 'Threshold = 200; MinRefAngle= 60, nAngles = 12, FlipAngleMap from B1-map with factor 1.165';
save(['GRASE_Results_B1_map_', FileName])
