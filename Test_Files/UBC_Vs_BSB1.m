function UBC_Vs_BSB1(Subject)
  cd ~/GRASE/GRASE_Results

  FileName = ['GRASE_Results_B1_map_',Subject,'.mat'];

  load(FileName, 'tf_mgrase')

  cd(['~/GRASE/B1_Maps/',Subject,'/'])

  FlipAngleMap = (159.5 * double(niftiread('rB1_Phase.nii'))) / (800);
  FlipAngleMap = flip(permute(FlipAngleMap,[2 1 3]),1);
  tic
  [maps,distributions,~] = T2map_Nima(tf_mgrase, 'Threshold', 200, 'T2Range', [0.015, 2], 'FlipAngleMap', FlipAngleMap);


  MWI = squeeze(squeeze(sum(distributions(:,:,:,1:40),4))./squeeze(sum(distributions(:,:,:,:),4)));

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
  Description = 'Threshold = 200;FlipAngleMap from B1-map with nominal angle of 159.5';
  save(['GRASE_Results_B1_map_', Subject])

end
