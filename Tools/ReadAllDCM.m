function Output = ReadAllDCM(directory,NumSlices)
   cd(directory) 
   Files = dir('*.*');
   disp('Reading DICOM Files...')
   disp(directory)
   s = 1;
   tic
   for k= 3:length(Files)
       temp = mod(k-2, NumSlices);
       if temp == 0
           Output(:,:,NumSlices,s) = dicomread(Files(k).name);
           s = s + 1;
       else
           Output(:,:,temp,s) = dicomread(Files(k).name);
       end
   end
   toc
   disp('Finished reading!')
end
