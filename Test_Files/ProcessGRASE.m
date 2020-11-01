function  ProcessGRASE( FileName, Options )
% function: basic process of GRASE MWI using the code given by UBC
% Postprocessing directory must be added to path 
%
% Inputs:
%	-FileName: .mat file containting a 4D matrix of data name either 'mgrase' or 'tf_mgrase'
%		'mgrase' is the raw data and 'tf_mgrase' is the tukey filtered version
%	-Options: a struture contaiong Read_Dir and Save_Dir which are the directory to load and save

if nargin < 2
	Options.Empty = true;
end

if ~isfield(Options, 'Read_Dir')
	Options.Read_Dir = '~/GRASE/GRASE_To_Do';
end
if ~isfield(Options, 'Save_Dir')
	Options.Save_Dir = '~/GRASE/GRASE_Results';
end
tic;


cd(Options.Read_Dir)
load(FileName, 'tf_mgrase','mgrase')

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

[maps,distributions,~] = T2map_SEcorr(tf_mgrase, 'Threshold', 200,'nT2', 60,'T2Range', [0.008, 2], 'MinRefAngle', 100);



MWI = squeeze(squeeze(sum(distributions(:,:,:,1:18),4))./squeeze(sum(distributions(:,:,:,:),4)));

MWI(isnan(MWI)) = 0;

clear mgrase tf_mgrase

runtime=toc;
cd(Options.Save_Dir)
Description = 'Threshold = 200; nT2 = 60, T2Range = 8ms to 2s, MinRefAngle = 100 degrees ';
save(strcat('GRASE_Results_', FileName))
end
