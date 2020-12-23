function ReProcessGRASE_SingleComponent(name)
addpath(genpath('~/GRASE'))
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
	  clear hfilt2 xs zs es fs ys
  end

  cd ~/GRASE/GRASE_Results/GRASE_B1_Map_Results/
  
  load(Subject, 'maps')
  
  [ys xs zs es] = size(tf_mgrase);
  tf_mgrase = reshape(tf_mgrase,  [ys*xs*zs, es]);
  Alpha = reshape(maps.alpha, [size(tf_mgrase,1),1]);
  clear maps ys xs zs es
  Maps.T2 = zeros(size(tf_mgrase,1),1);
  Maps.Alpha = zeros(size(tf_mgrase,1),1);
  Maps.Res = zeros(size(tf_mgrase,1),1);
  Maps_B1 = Maps;
  Threshold = 200;
  
  tic
  parfor i = 1:size(tf_mgrase,1)
	if tf_mgrase(i,1) > Threshold
		[Maps.T2, Maps.B1, Maps.Res, ~] = Single_Component_T2_B1(tf_mgrase(i,:));
		opt.B1 = Alpha(i);
		[Maps_B1.T2, Maps_B1.B1, Maps_B1.Res, ~] = Single_Component_T2(tf_mgrase(i,:), opt);
	end
  end
  
  runtime=toc;

  clear tf_mgrase maps tf_mgrase mgrase Alpha

  cd ~/GRASE/GRASE_Results/SingleComponent/
  Description = 'Using the flip angle map from the GRASE B1 Map Results';
  save(strcat('GRASE_Results_SC_', Subject))
catch ME
	disp(ME.message)
	end
end
