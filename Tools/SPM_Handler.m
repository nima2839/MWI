classdef SPM_Handler
% Author: Nima
%%
% The objective is to use SPM12 functionality outside of GUI;
% Must add "spm12" directory to path of MATLAB before use;
% This class needs to write a MATLAB batch file in the output directory;
% For now, following the batch and scripts format of the GUI;

	methods(Static = true)
		function out = Coregister(Options)
			% Carries out the same function as GUI, to coregister (reslice and realign)
			% Options is a struct with the following fields:
			%	'ref' : Path of reference file (remains stationary)
			%	'source': Path of source file (jiggled to match the ref)
			%	'other': Path of the other file that should be coregistered into the same space as source file;
			%
			%
			%-------------------------------------------------------------------
			% Making sure filenames include the '.nii' extension
			if isempty(strfind(Options.ref,'.nii'))
				Options.ref = strcat(Options.ref, '.nii');
			end
			
			if isempty(strfind(Options.source,'.nii'))
				Options.source = strcat(Options.source, '.nii');
			end
			
			if isfield(Options, 'other')
				if isempty(strfind(Options.other,'.nii'))
					Options.other = strcat(Options.other, '.nii');
				end
			end

			if isempty(which('spm'))
				error(' Please make sure SPM12 is added to MATLAB path!');
			end

			matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(Options.ref);
			matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(Options.source);
			if isfield(Options, 'other')
				matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(Options.other);
			end
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = ...
							[0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
			try
				spm_jobman('run', matlabbatch);
			catch ME
				disp('SPM_Handler.Coregister failed!');
				rethrow(ME);
			end
		end

		function out = DICOM_TO_NIFTI(Options)
			% Reads DICOM files and writes them in a .nii format
			% Options is a struct with the following fields:
			%	'ReadPath': path of the dicom files;
			%	'OutDir': output directory, if not provided it would write into the current directory;
			%	'Num_Of_Files': If provided, will only read n first DICOM files, otherwise it will read all;
			%	'Name': name of output NIFTI file
			% On success returns a cell containing the name of the written files;
			%
			%-------------------------------------------------------------------
			if isempty(which('spm'))
				error(' Please make sure SPM12 is added to MATLAB path!');
			end

			if ~isfield(Options, 'Num_Of_Files')
				Options.Num_Of_Files = 0;
			end

			if ~isfield(Options, 'OutDir')
				Options.OutDir = pwd;
			end

			try
				Headers =  spm_dicom_headers(SPM_Handler.Read_File_Names(Options.ReadPath, Options.Num_Of_Files));
				spm_output = spm_dicom_convert(Headers,'all','flat','nii', Options.OutDir);
				% renaming the output files
				if isfield(Options, 'Name')
					% cd(Options.OutDir)
					n =  numel(spm_output.files);
					ext = "";
					out = cell(n, 1);
					for i = 1:n
						if n > 1
							ext = strcat("_Echo_", string(i));
						end
						out{i} = strcat(Options.OutDir,'\', Options.Name, ext, ".nii");
						movefile(spm_output.files{i}, out{i});
					end
				else
					out  = spm_output.files;
				end
			catch ME
				disp('SPM_Handler.DICOM_TO_NIFTI failed!');
				rethrow(ME);
			end
		end

		function FileNames = Read_File_Names(ReadDir, n, ext)
			% Reads file names to add to use in Job file
			% 'n' is the number of files in directory
			% 'ext' is the extention of the files
			%
			%-------------------------------------------------------------------

			if nargin < 3
				ext = '.dcm';
			end

			FileNames = "";
			files = dir(strcat(ReadDir, '/*', ext));
			
			if nargin < 2
				n = 0;
			end

			if n == 0 
				n = length(files);
			end
			FileNames = cell(1,n);
			for i = 1:n
				FileNames{i} = char(strcat(ReadDir, "/", files(i).name));
			end
		end

	end

end
