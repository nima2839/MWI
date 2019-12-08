sliceNum= 2;
img(:,:) = Final_MWI(:,:,sliceNum);
Dist = distributions;
time=logspace(log10(0.015),log10(2),200);
figure
subplot(3,1,1)
imagesc(img)
while ~isempty(findobj(allchild(groot), 'flat', 'type', 'figure', 'number', 1))
    try
    [x,y] = ginput(1);
    subplot(3,1,2)
    semilogx(time,squeeze(Dist(floor(y),floor(x),sliceNum,:)))
    grid
    subplot(3,1,3)
    polt(squeeze(Mag(floor(y),floor(x),sliceNum,:)))
    grid
    catch
    end
end