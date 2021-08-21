function  Results = ProcessGRASE(Options, Save_Name)
  % Processes GRASE MWI directly from DICOM files;
  % Requires '\MWI' and '\spm12' to be added to the working path;
  %
  % Supplied flip angle method is only performed when dirctory to DICOM files of the B1 phase maps are specified; 
  % A "Description" is retured along with the result in the case of using a flip angle map
  %
  % Output:
  % "Results" is a structure containing all results;
  %
  % Input arguments:
  % "Save_Name" if provided it would save results under the given name; (optional)
  % "Options" is a structure containing the following fields:
  %   - 'GRASE_Dir' is the path to GRASE DICOM files;
  %   - 'Save_Dir' is the path to where the data should be saved; (optional)
  %   - 'B1Phase_Dir' is the path to B1 phase DICOM files (used to aquire flip angle map); (optional)
  %   - 'B1Mag_Dir' is the path to B1 magnitude DICOM files (used for coregisteration purposes); (optional)
  %   - 'B1_Scale' is the scaling factor required to noramalize B1 Phase data (800); (optional) 
  %   - 'Nominal_Angle' is the scaling factor to aquire flip angle map from the normalized B1 map; (optional)
  %   - 'Find_Mask' is a function pointer that user specifies to find mask using only the 4D magnitude data
  %       + Ex: Options.Find_Mask = @(Mag) Create_Mask(Mag(:,:,:,1), 0.35); % "Create_Mask" can be found under "\MWI\Tools\"
  %       + (optional)
  %   - 'Find_NominalAngle' is a function pointer that accepts two input args: (B1, Maps)
  %       + Where 'B1' is the normalized coregistered B1 Phase map, and 'Maps' is the results of the estimated FA method; 
  %       + (handled by this funciton);
  %       + If the function fails with an exception it will skip supplied B1+ method;
  %       + (optional)
  %
    % Checking and intializing the input
    if ~isfield(Options, 'GRASE_Dir')
      error('GRASE_Dir is required!');
    end
  
    if isfield(Options, 'Save_Dir')
      cd(Options.Save_Dir);
    end
  
    if ~isfield(Options, 'B1_Scale')
      Options.B1_Scale = 800;
    end
  
    if ~isfield(Options, 'Nominal_Angle') 
      if ~isfield(Options, 'Find_NominalAngle')
        Options.Find_NominalAngle = @(B1, Maps) Estimate_NominalAngle(B1, Maps); % calculated empirically
      end
    end
  
    FlipAngleMap = [];
    % Reading GRASE images along with proper information for NIFTI format
    [Image, nifti_info] = ReadAllDCM(Options.GRASE_Dir, 'spm');
  
    if isfield(Options, 'Find_Mask')
      try
        Mask = Options.Find_Mask(Image);
        sd = size(Image);
        temp = reshape(Image, [numel(Mask), sd(4)]);
        idx = Mask > 0;
        temp(idx,:) = 0;
        Image = reshape(temp, sd);
      catch ME
        disp(ME);
        disp('Masking process failed!');
      end
    end
  
    % MWI analysis estimated FA method
    [Maps, Dist, ~] = T2map_Nima(Image);
    Results.Maps = Maps;
    Results.Dist = Dist;
    Results.nifti_info = nifti_info;
  
    
    if isfield(Options, 'B1Phase_Dir')
      % reading dicom files for coregisteration
      opt.ReadPath = Options.B1Phase_Dir;
      opt.Name = strcat('B1_Phase_temp', date);
      temp_files{1} = opt.Name;
      SPM_Handler.DICOM_TO_NIFTI(opt);
      if isfield(Options, 'B1Mag_Dir')
        opt.ReadPath = Options.B1Mag_Dir;
        opt.Name = strcat('B1_Mag_temp', date);
        temp_files{2} = opt.Name;
        SPM_Handler.DICOM_TO_NIFTI(opt);
      end
      opt.ReadPath = Options.GRASE_Dir;
      opt.Name = strcat('GRASE_Mag_temp', date);
      temp_files{numel(temp_files) + 1} = opt.Name;
      opt.Num_Of_Files = size(Image, 3);
      SPM_Handler.DICOM_TO_NIFTI(opt);
      % coregisteration part
      opt.ref = strcat(temp_files{end}, '.nii');
      if isfield(Options, 'B1Mag_Dir')
        temp_files{4} = strcat('r', temp_files{1});
        temp_files{5} = strcat('r', temp_files{2});
        opt.source = strcat(temp_files{2}, '.nii');
        opt.other = strcat(temp_files{1}, '.nii');
      else
        temp_files{3} = strcat('r', temp_files{1});
        opt.source = temp_files{1};
      end
      try
        SPM_Handler.Coregister(opt);
        % Reading the coregistered B1 Map
        B1_Normalized = double(niftiread(strcat('r',temp_files{1},'.nii'))) / Options.B1_Scale;
      catch ME
        disp("Failed in the coregistration process, returning without any results");
        rethrow(ME)
      end
      
      % If the estimation process fails Supplied B1 method will be skiped 
      Skip_Supplied_B1 = false;
      if ~isfield(Options, 'Find_NominalAngle')
        try 
          Options.Nominal_Angle = Options.Find_NominalAngle(B1_Normalized, Maps);
        catch ME
          disp(ME);
          disp('Nominal Angle estimation failed! Skipping B1 supplied method!');
          Skip_Supplied_B1 = true;
        end
      end
  
      if Skip_Supplied_B1
        FlipAngleMap = [];
      else
        FlipAngleMap = Options.Nominal_Angle * B1_Normalized;
      end
        % clearing the temp files
      for i = 1:numel(temp_files)
        delete(strcat(temp_files{i}, '.nii'));
      end
    end
    
  
    if ~isempty(FlipAngleMap)
      [Maps_B1, Dist_B1, ~] = T2map_Nima(Image, 'FlipAngleMap', FlipAngleMap);
      Results.Description = 'The %_B1% extention specifies the supplied B1 method resutls';
      Results.Maps_B1 = Maps_B1;
      Results.Dist_B1 = Dist_B1;
      Results.Nominal_Angle =  Options.Nominal_Angle;
    end
  
    if nargin > 1
      save(Save_Name, 'Results');
    end
  end
  