function ReProcessGRASE_SingleComponent(Subject)
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
  
  load(strcat("GRASE_Results_B1_map_",Subject), 'maps')
  
  [ys xs zs es] = size(tf_mgrase);
  tf_mgrase = reshape(tf_mgrase,  [ys*xs*zs, es]);
  Alpha = reshape(maps.alpha, [size(tf_mgrase,1),1]);
  clear maps ys xs zs es
  T2 = zeros(size(tf_mgrase,1),1);
  B1 = zeros(size(tf_mgrase,1),1);
  Res = zeros(size(tf_mgrase,1),1);
  T2_B1 = zeros(size(tf_mgrase,1),1);
  Res_B1 = zeros(size(tf_mgrase,1),1);
  Threshold = 200;
  opt.B1 = 180;
  tic
  parfor i = 1:size(tf_mgrase,1)
	temp = opt;
	if tf_mgrase(i,1) > Threshold
		[T2(i), B1(i), Res(i), ~] = Single_Component_T2_B1(tf_mgrase(i,:));
		temp.B1 = Alpha(i);
		[T2_B1(i), ~, Res_B1(i), ~] = Single_Component_T2(tf_mgrase(i,:), temp);
	end
  end
  
  Maps.T2 = T2;
  Maps.B1 = B1;
  Maps.Res = Res;
  Maps_B1.T2 = T2_B1;
  Maps_B1.B1 = Alpha;
  Maps_B1.Res = Res_B1;
  runtime=toc;

  clear tf_mgrase maps tf_mgrase mgrase Alpha T2 T2_B1 B1 Res Res_B1 Alpha

  cd ~/GRASE/GRASE_Results/SingleComponent/
  Description = 'Using the flip angle map from the GRASE B1 Map Results';
  save(strcat('GRASE_Results_SC_', Subject))
catch ME
	disp(ME.message)
	end
end
