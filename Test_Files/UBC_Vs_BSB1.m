function UBC_Vs_BSB1(Subject, Nominal_Angle)
addpath(genpath('~/GRASE'))
if nargin < 2
	Nominal_Angle = 156.2;
end
disp(['Processing BSB1 on :',Subject])
try
  cd ~/GRASE/GRASE_To_Do/

  load(Subject, 'tf_mgrase','mgrase')
  
  if ~exist('tf_mgrase')
	  if ~exist('mgrase')
		disp(strcat("Processing **",FileName,"** failed: could not load data!"));
		return;
	  end
	  fs = 1/3;
	  [ys xs zs es]=size(mgrase);
	  tf_mgrase = zeros(ys,xs,zs,es);
	  hfilt2=tukeywin(ys,fs)*tukeywin(xs,fs)';

	  for i = 1:es
		  for j = 1:zs

			  tf_mgrase(:,:,j,i) = abs(ifft2c(fft2c(mgrase(:,:,j,i)).*hfilt2));

		  end
	  end
end

  cd(strcat('~/GRASE/B1_Maps/',Subject,'/'))

  FlipAngleMap = (Nominal_Angle * double(niftiread('rB1_Phase.nii'))) / (800);
  FlipAngleMap = flip(permute(FlipAngleMap,[2 1 3]),1);
  tic
  [maps,distributions,~] = T2map_Nima(tf_mgrase, 'Threshold', 200, 'T2Range', [0.008, 2], 'FlipAngleMap', FlipAngleMap,'nT2', 60);


  MWI = squeeze(squeeze(sum(distributions(:,:,:,1:18),4))./squeeze(sum(distributions(:,:,:,:),4)));

  MWI(isnan(MWI)) = 0;
  runtime=toc;

  clear tf_mgrase MWI_1 FlipAngleMap tf_mgrase mgrase

  cd ~/GRASE/GRASE_Results/GRASE_B1_Map_Results/
  Description = 'Threshold = 200; FlipAngleMap from B1-map with nominal angle provided by user, nT2 = 60, T2Range = 8ms to 2 S';
  save(strcat('GRASE_Results_B1_map_', Subject))
catch ME
	disp(ME.message)
	end
end
