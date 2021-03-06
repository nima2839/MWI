function [R , res] = SingleComponentFitting(Data,Info)
% Finds sigle component relaxation value
% Date: 6/3/2016
% Author:  Nima
% Outputs:
%   R: relation term
%   r: residue
if ~isfield(Info, 'ws')
    Info.ws = 0.4;
end
if max(Data(:)) ~= 0
    t = (Info.FirstTE):(Info.EchoSpacing):(Info.FirstTE + (length(Data)-1) * Info.EchoSpacing);
    if size(t) ~= size(Data)
        t= t';
    end
    % Applying a Low-Pass Filter
    fd = double(Data)/max(Data(:));%lowpass(Data/max(Data(:)),Info.ws);
    % removing the zero elements
    temp = fd;
    temp(fd==0) = 2;
    temp(fd~=2) = 0;
    index = find(temp,1,'first');
    if isempty(index)
        index = length(fd);
    end
    L = polyfit(t,log(fd(1:index)),1);
    R = abs(L(1));
    tmp = exp(-R*t);
    tmp = tmp/max(tmp(:)) -fd;
    tmp = tmp.^2;
    res = sum(tmp(:));
else
    R = inf;
    res = 0;
end
end
