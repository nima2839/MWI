sd = size(img_cmb);
for e = 1:sd(4)
    temp = zeros(sd(1:3));
	temp(:,:,:) = img_cmb(:,:,:,e);
	temp = fftshift(fftn(temp));
	temp = padarray(temp,[sd(1)/2,sd(2)/2,0],'both');
	temp = ifftn(fftshift(temp));
	padded_img(:,:,:,e) = temp;
end
