function Mask = Create_Mask(Mag, Fraction)
  % Creates mask using fsl package
  % Mag must be a 3D matrix
  if nargin < 2
    Fraction = 0.5;
  end
  OriginalPath = pwd;
  disp('Creating temp environment!');
  unix('mkdir Create_Mask_temp_File');
  cd('Create_Mask_temp_File')

  disp('Writing Mag.nii file')
  niftiwrite(Mag, 'Mag');

  disp('Creating Mask ...');
  unix(char(['bet Mag.nii mask.nii -m -f ', string(Fraction)]));
  unix('gunzip -f mask_mask.nii.gz');

  disp('Mask Created!');
  Mask = double(niftiread('mask_mask.nii'));
  disp('Cleaning up...')
  cd(OriginalPath)
  unix('rm -r Create_Mask_temp_File');

  disp('Done!')
end
