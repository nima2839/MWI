classdef SPM_Handler
% Author: Nima
%%
% The objective is to use SPM12 functionality outside of GUI;
% Must add "spm12" directory to path of MATLAB before use;
% This class needs to write a MATLAB batch file in the output directory;
% For now, following the batch and scripts format of the GUI;

	methods(Static = true)
		function out = Coregister(Options)
			% Carries out the same function as GUI, to coregister (reslice and realign) two images and other images
			% Options is a struct with the following fields:
			%	ref : Path of reference file (remains stationary)
			%	source: Path of source file (jiggled to match the ref)
			%	Other: Path of all other files, concatenated by [], not by strcat!
			%		and they must have single quatations arround each path!
			% On Failure out is the excption file, and on success it reurns 0;
			%
			%-------------------------------------------------------------------
			MyJob = ["matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {",...
					strcat("'",Options.ref, "'};"),...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.source = {",...
					strcat("'",Options.source, "'};"),...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.other = {", Options.Other,"};",...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';",...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];",...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];",...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];",...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;",...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];",...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;",...
					"matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';"
					];
			try
				SPM_Handler.Job_Handler(MyJob);
				out = 0;
			catch ME
				disp('SPM_Handler.Coregister failed!');
				out = ME;
			end
		end

		function out = DICOM_TO_NIFTI(Options)
			% Reads DICOM files and writes them in a .nii format
			% Options is a struct with the following fields:
			%	ReadPath: path of the dicom files;
			%	OutDir: output directory, if not provided it would write into the current directory;
			%	Num_Of_Files: If provided, will only read n first DICOM files, otherwise it will read all;
			%	Name: name of output NIFTI file
			% On Failure out is the excption file, and on success it reurns 0;
			%
			%-------------------------------------------------------------------
			if ~isfield(Options, 'Num_Of_Files')
				Options.Num_Of_Files = 0;
			end

			if ~isfield(Options, 'OutDir')
				Options.OutDir = pwd;
			end


			MyJob = ["matlabbatch{1}.spm.util.import.dicom.data = {",...
					SPM_Handler.Read_File_Names(Options.ReadPath, Options.Num_Of_Files),...
					"};",...
					"matlabbatch{1}.spm.util.import.dicom.root = 'flat';",...
					"matlabbatch{1}.spm.util.import.dicom.outdir = {", strcat("'",Options.OutDir,"'") ,...
					"};",...
					"matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';",...
					"matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';",...
					"matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;",...
					"matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;"
					];

			try
				SPM_Handler.Job_Handler(MyJob);
				out = 0;
				% renaming the output file
				if isfield(Options, 'Name')
					cd(Options.OutDir)
					file = dir('s*.nii');
					if length(file) > 1
						[~,idx] = sort([file.datenum]);
						file = file(idx(end));
					end
					movefile(file.name, strcat(Options.Name, ".nii"));
				end
			catch ME
				disp('SPM_Handler.DICOM_TO_NIFTI failed!');
				out = ME;
			end
		end

		function FileNames = Read_File_Names(ReadDir, n)
			% Reads file names to add to use in Job file
			% n input is the number of files in directory
			FileNames = "";
			files = dir(strcat(ReadDir,'/*.dcm'));
			if n == 0
				n = length(files);
			end

			for i = 1:n
				FileNames = [FileNames, strcat("'",ReadDir, "/", files(i).name, "'")];
			end
		end

		function Job_Handler(Job)
			% Writes down "Job" input into a file and then runs spm12 code
			% Job is a string array in the same format as spm12 scripts
			%
			%-------------------------------------------------------------------
			% Creating a new file (or clear file) to write the job into
			FileName = 'tempSPM12_Job.m';
			Header = "% This file is generated by SPM_Handler object!";
			SPM_Handler.Write_To_File(FileName, Header, 'wt');


			% Writing Job into file
			SPM_Handler.Write_To_File(FileName, Job, 'at');

			% Runnig the Job on SPM12
			spm('defaults', 'FMRI');
			spm_jobman('run', [pwd, '/', FileName]);

			% clearing Job file
			delete(FileName);
		end

		function Write_To_File(FileName, text, permission)
			% Writes text into file and throws excption on failure
			FID = fopen(FileName, permission);

			if FID == -1
				ME = MException('MyComponent:fopenFailed',['MATLAB could not open ', FileName, '!']);
				throw(ME);
			end

			if fprintf(FID,'%s\n', text)  == 0
				ME =  MException('MyComponent:fprintf_Failure',...
						'MATLAB failed to write in job file!');
				fclose(FID);
				throw(ME);
			end
			fclose(FID);
		end

	end

end
