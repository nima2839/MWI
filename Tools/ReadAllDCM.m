function Output = ReadAllDCM(directory,NumSlices)
	cd(directory) 
	Files = dir('*.*'); % First two are the current and parent directories, the rest are assumend to be DICOM files
	disp('Reading DICOM Files...')
	disp(directory)
	s = 1;
	tic
	if nargin > 1 % this assumes files are in stacks of NumSlices
		for k = 3:length(Files)
			temp = mod(k-2, NumSlices);
			if temp == 0
				Output(:,:,NumSlices,s) = double(dicomread(Files(k).name));
				s = s + 1;
			else
				Output(:,:,temp,s) = double(dicomread(Files(k).name));
			end
		end
	else % otherwise it will return a cell sorted by EchoNumber and SliceLocation of DICOM headers
		for i = 3:numel(Files)
			temp_info =  dicominfo(Files(i).name);
			SliceLocations(i-2) = temp_info.SliceLocation;
			EchoNumbers(i-2) = temp_info.EchoNumber;
			D_images(:,:, i-2) = double(dicomread(Files(i).name));
		end
		ETL = max(EchoNumbers(:));
		Output = cell(1,ETL);
		for i = 1:ETL
			idx = find(EchoNumbers == i);
			[~, idx_sort] = sort(SliceLocations(idx));
			Output{i} = D_images(:,:, idx(idx_sort));
		end
	end
	toc
	disp('Finished reading!')
end