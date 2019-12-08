function w = Hann3D(nx,ny,nz)
% creats a 3D Hanning window
	w = bsxfun(@times,bsxfun(@times,hann(nx),hann(ny).'),permute(hann(nz),[3 2 1])); 
end
