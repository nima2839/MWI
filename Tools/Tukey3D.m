function w = Tukey3D(nx,ny,nz,r)
% creats a 3D Tukey window
	w = bsxfun(@times,bsxfun(@times,tukeywin(nx,r),tukeywin(ny,r).'),permute(tukeywin(nz,r),[3 2 1])); 
end
