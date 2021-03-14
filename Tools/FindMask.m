function Mask = FindMask(Mag, opt)
% Creates a brain mask for the input magnitude data using FSL bet (brain extraction func.)
% FSL path must be added to the current terminal for this function!
% /usr/local/fsl/bin/bet in out  -f 0.35 -g 0 -n -m
	if nargin < 2
		opt = "-f 0.3 -g 0 -n -m";
	end
	orig = pwd;
	% Create a temp file to save and read files
	if system('mkdir temp_FindMask_generated_dir')
		error('Could not create a directory to perform the process!')
	end
	cd('temp_FindMask_generated_dir')
	if numel(size(Mag)) > 4
		Mag = Mag(:,:,:,1);
	end
	niftiwrite(Mag,'Mag')
	% call FSL bet for brain extraction
	msg = '';
	if system(strcat("bet Mag.nii out ",opt));
		msg = 'There was an error issued by brain extraction process!';
		myfunc_Clean(msg);
		error(msg);
	else
		if system('gunzip out_mask.nii')
			msg = 'There was an issue in decompressing the fsl results!';
			myfunc_Clean(msg);
			error(msg);
		else
			Mask = double(niftiread('out_mask.nii'));
		end
	end
	myfunc_Clean('');
	function myfunc_Clean(msg)
		cd(orig)
		if system('rm -r temp_FindMask_generated_dir')
			error(strcat(msg,'+ FindMask() could not finish cleaning up files in the "temp_FindMask_generated_dir"!'));
		end
	end
end