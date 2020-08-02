function UBC_Vs_BSB1(Subject)
disp(['Processing BSB1 on :',Subject])
try
  cd ~/GRASE/GRASE_To_Do/

  load(Subject, 'tf_mgrase','mgrase')
  
  if ~exist('tf_mgrase')
    tf_mgrase = mgrase;
  end

  cd(strcat('~/GRASE/B1_Maps/',Subject,'/'))

  FlipAngleMap = (156.2 * double(niftiread('rB1_Phase.nii'))) / (800);
  FlipAngleMap = flip(permute(FlipAngleMap,[2 1 3]),1);
  tic
  [maps,distributions,~] = T2map_Nima(tf_mgrase, 'Threshold', 200, 'T2Range', [0.008, 2], 'FlipAngleMap', FlipAngleMap,'nT2', 60);


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

  runtime=toc;

  clear tf_mgrase MWI MWI_1 FlipAngleMap

  cd ~/GRASE/GRASE_Results/GRASE_B1_Map_Results/
  Description = 'Threshold = 200; FlipAngleMap from B1-map with nominal angle of 156.2, nT2 = 60, T2Range = 8ms to 2 S';
  save(strcat('GRASE_Results_B1_map_', Subject))
catch ME
	disp(ME.message)
	end
end
