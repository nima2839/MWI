function Mask = FindMask(Mag)
% Creates a brain mask for the input magnitude data
% Mag: Magnitude data (2D/3D)
% Threshold: Data less than this value will be removed

SM = size(Mag);
% Disk will erode about 10% of pixels
SE = strel('disk',round(min(SM(1),SM(2))*0.1));

if      length(SM) < 2
    error('Dimensions of Magnitude data can not be less than 2!')
elseif  length(SM) > 3
    error('Too many dimensions!')
elseif  length(SM) == 2
    Mask = imfill(Mag > mean(mag), 'holes');
    Mask = imerode(Mask, SE);
elseif  length(SM) == 3
    Mask = zeros(SM);
    for k = 1:SM(3)
        Mask(:,:,k) = imfill(Mag(:,:,k) > mean(Mag(:,:,k)), 'holes');
        Mask(:,:,k) = imfill(Mask(:,:,k), 'holes');
        Mask(:,:,k) = imerode(Mask(:,:,k), SE);
    end
end
Mask = logical(Mask);
end