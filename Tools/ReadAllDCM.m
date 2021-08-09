function [Output, nifti_info] = ReadAllDCM(Path, Options)
	% Reads all the DICOM files (files ending with '.dcm') in the 'Path' and returns:
	%	'Output': 4D matrix;
	%	'nifti_info': when using 'spm' method returns the proper NIFTI information corresponding to the first echo data; 
	%				+ When using 'Normal' method returns an empty array;
	%
	% Input arguemnts:
	%	'Path': path to the DICOM images
	%	'Options': a structure containing the following fields:
	%		-'NumSlices': number of slices; if not provided code will decide to read the files based in the method;
	%		-'Method': the method to read the DICOM images->
	%			->'Normal': the default method; does not return any information and no dimension appropration is performed
	%			->'spm': uses spm12 toolbox to initialize NIFTI information and matches the dimention accordingly;
	%					+ This method requires the toolbox to be added to MATLAB path.
	%		+ 'Options' can be a sting specifying the 'Method' field as well!
	%

	if nargin < 2
		Options.Method = 'Normal';
	end

	if ~isstruct(Options)
		if ischar(char(Options))
			opt.Method = Options;
			[Output, nifti_info] = ReadAllDCM(Path, opt);
			return;
		else
			error('Invalid input argument!')
		end
	end

	if ~isfield(Options, 'Method')
		Options.Method = 'Normal';
	end

	nifti_info = [];

	if strcmp(Options.Method, 'Normal')
		orig_dir = pwd;
		cd(Path); 
		Files = dir('*.dcm'); % First two are the current and parent directories, the rest are assumend to be DICOM files
		disp('Reading DICOM Files...')
		disp(Path)
		s = 1;
		tic
		if nargin > 1 % this assumes files are in stacks of NumSlices
			for k = 1:length(Files)
				temp = mod(k, NumSlices);
				if temp == 0
					Output(:,:,NumSlices,s) = double(dicomread(Files(k).name));
					s = s + 1;
				else
					Output(:,:,temp,s) = double(dicomread(Files(k).name));
				end
			end
		else % otherwise it will return a cell sorted by EchoNumber and SliceLocation of DICOM headers
			for i = 1:numel(Files)
				temp_info =  dicominfo(Files(i).name);
				SliceLocations(i) = temp_info.SliceLocation;
				EchoNumbers(i) = temp_info.EchoNumber;
				D_images(:,:, i) = double(dicomread(Files(i).name));
			end
			ETL = max(EchoNumbers(:));
			nrows = size(D_images,1);
			ncols = size(D_images,2);
			Output = zeros(nrows, ncols, numel(Files)/ETL, ETL);
			for i = 1:ETL
				idx = find(EchoNumbers == i);
				[~, idx_sort] = sort(SliceLocations(idx));
				Output(:,:,:,i) = D_images(:,:, idx(idx_sort));
			end
		end
		% Matching the dimension to the format of the nifti output of spm12 toolbox
		Output = permute(flip(Output,1), [2, 1, 3, 4]);
		toc
		cd(orig_dir);
		disp('Finished reading!')
	elseif  strcmp(Options.Method, 'spm')
		temp_name = strcat('temp_File',date);
		if ~isfield(Options,'NumSlices')
			[Output, ~] = ReadAllDCM(Path);	
		else
			opt.NumSlices = Options.NumSlices;
			[Output, ~] = ReadAllDCM(Path, opt);
		end
		
		try
			opt.ReadPath = Path;
			opt.Num_Of_Files = size(Output,4);
			opt.Name = temp_name;
			SPM_Handler.DICOM_TO_NIFTI(opt);
			temp_name = strcat(temp_name, '.nii');	
			nifti_info = niftiinfo(temp_name);
			nifti_info.Datatype = 'double';
			delete(temp_name);
		catch ME
			disp('The following error occured when reading the nifti-info:');
			disp(ME);
		end

	else
		error(strcat("Method '",Options.Method,"' is not valid!"))
	end
end