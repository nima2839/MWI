sliceNum= 16;
img(:,:) = maps.ggm(:,:,sliceNum);
img(img <0.035 | img >0.08) = nan;
Dist = distributions;
time=logspace(log10(0.008),log10(2),size(Dist,4));
figure
subplot(1,2,1)
imagesc(img)
while ~isempty(findobj(allchild(groot), 'flat', 'type', 'figure', 'number', 1))
    try
    [x,y] = ginput(1);
    subplot(1,2,2)
    semilogx(time,squeeze(Dist(floor(y),floor(x),sliceNum,:)))
    grid
    %subplot(3,1,3)
    %plot(squeeze(Mag(floor(y),floor(x),sliceNum,:)))
    %grid
    catch ME
        disp(ME.message);
    end
end