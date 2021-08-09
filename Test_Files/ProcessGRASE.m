function  ProcessGRASE(Save_Name, MyInfo)
% Processes GRASE MWI directly from DICOM files;
% Requires '\MWI' and '\spm12' to be added to the working path;
% Estimated flip angle method will always be performed;
% Supplied flip angle method is only performed when dirctory to DICOM files of the B1 phase maps are specified; 
% A "Description" is saved along with the result in the case of using a flip angle map
%
% "Save_Name": name of the file that will be saved;
% "MyInfo" is a structure containing the following fields:
%   - 'GRASE_Dir': path to GRASE DICOM files;
%   - 'Save_Dir': path to where the data should be saved; (optional)
%   - 'B1Phase_Dir': path to B1 phase DICOM files (used to aquire flip angle map); (optional)
%   - 'B1Mag_Dir': path to B1 magnitude DICOM files (used for coregisteration purposes); (optional)
%   - 'B1_Scale': the scaling factor required to noramalize B1 Phase data (800); (optional) 
%   - 'Nominal_Angle': The scaling factor to aquire flip angle map from the normalized B1 map; (optional but crucial)
%
  % Checking and intializing the input
  if ~isfield(MyInfo, 'GRASE_Dir')
    error("GRASE_Dir is required!")
  end

  if isfield(MyInfo, 'Save_Dir')
    cd(MyInfo.Save_Dir);
  end

  if ~isfield(MyInfo, 'B1_Scale')
    MyInfo.B1_Scale = 800;
  end

  if ~isfield(MyInfo, 'Nominal_Angle')
    MyInfo.Nominal_Angle = 156.5; % calculated empirically
  end

  FlipAngleMap = [];
  % Reading GRASE images along with proper information for NIFTI format
  [Image, nifti_info] = ReadAllDCM(MyInfo.GRASE_Dir, 'spm');
  
  if isfield(MyInfo, 'B1Phase_Dir')
    % reading dicom files for coregisteration
    opt.ReadPath = MyInfo.B1Phase_Dir;
    opt.Name = strcat('B1_Phase_temp', date)
    temp_files{1} = opt.Name;
    SPM_Handler.DICOM_TO_NIFTI(opt);
    if isfield(MyInfo, 'B1Mag_Dir')
      opt.ReadPath = MyInfo.B1Mag_Dir;
      opt.Name = strcat('B1_Mag_temp', date)
      temp_files{2} = opt.Name;
      SPM_Handler.DICOM_TO_NIFTI(opt);
    end
    opt.ReadPath = MyInfo.GRASE_Dir;
    opt.Name = strcat('GRASE_Mag_temp', date)
    temp_files{numel(temp_files) + 1} = opt.Name;
    opt.Num_Of_Files = size(Image, 3);
    SPM_Handler.DICOM_TO_NIFTI(opt);
    % coregisteration part
    opt.ref = strcat(temp_files{end}, '.nii');
    if isfield(MyInfo, 'B1Mag_Dir')
      temp_files{4} = strcat('r', temp_files{1});
      temp_files{5} = strcat('r', temp_files{2});
      opt.source = strcat(temp_files{2}, '.nii');
      opt.other = strcat(temp_files{1}, '.nii');
    else
      temp_files{3} = strcat('r', temp_files{1});
      opt.source = temp_files{1};
    end
    SPM_Handler.Coregister(opt);
    FlipAngleMap = MyInfo.Nominal_Angle * double(niftiread(strcat('r',temp_files{1},'.nii'))) / MyInfo.B1_Scale;
    % clearing the temp files
    for i = 1:numel(temp_files)
      delete(strcat(temp_files{i}, '.nii'));
    end
  end
  
  % MWI analysis
  [Maps, Dist, ~] = T2map_Nima(Image);
  if isempty(FlipAngleMap)
    save(Save_Name, 'Maps', 'Dist', 'nifti_info');
  else
    [Maps_B1, Dist_B1, ~] = T2map_Nima(Image, 'FlipAngleMap', FlipAngleMap);
    Description = "The '_B1' extention specifies the supplied B1 method resutls";
    save(Save_Name, 'Maps', 'Dist', 'Maps_B1', 'Dist_B1', 'Description','nifti_info');
  end

end
