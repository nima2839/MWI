function [Distribution,Times,Res] = NNLS_T_Batch_Fitting(Data,EchoTime_obj,Option_obj,LFG)
%
%   Author: Nima
%   Date: 07/2019
% Inputs:
%   Data: Contains the 4D data of all echeos
%   EchoTime_obj.FirstEchoTime: time of the first echo
%   EchoTime_obj.EchoSpacing: spacing time difference between echoes
%   Option_obj.T_range: Lower and upper bounds for T values
%   Option_obj.Number_T: Number of T values

tic
sizeData = size(Data);
if length(sizeData) ~= 4
	error('Data must be 4D!');
end


if ~isfield(Option_obj,'A') || isempty(Option_obj.A)
	T = linspace(Option_obj.T_range(1),Option_obj.T_range(2),...
	      Option_obj.Number_T);
	  % add a long T value
	  Option_obj.Number_T = Option_obj.Number_T + 1;
	  T(Option_obj.Number_T) = 1.5;
	sT = length(T);
	Option_obj.A = zeros(Number_echos,sT);
	for i = 1:sT
			Option_obj.A(:,i) = transpose(exp(-EchoTimes/T(i)));
	end
end
Distribution = zeros(sizeData(1), sizeData(2),...
                sizeData(3), Option_obj.Number_T + 1);

disp('Startiing fitting process...')
Opt_obj = Option_obj;
parfor k = 1:sizeData(3)
	temp_dist = zeros(sizeData(1),sizeData(2),1,Option_obj.Number_T + 1);
	temp_res = zeros(sizeData(1:2));
    for j = 1:sizeData(2)
        for i = 1:sizeData(1)
					if Data(i,j,k,1) > 0
            temp = Data(i,j,k,:);
			% The brain mask comes in handy here!
						if temp(1) > 0
							[~ ,temp_dist(i,j,1,:), temp_res(i,j)] = NNLS_T_Fitting(temp(:),EchoTime_obj,Opt_obj);
						end
					end
        end
    end
		Res(:,:,k) = temp_res(:,:);
		Distribution(:,:,k,:) = temp_dist(:,:,1,:);
end
Times = linspace(Option_obj.T_range(1),Option_obj.T_range(2),...
		Option_obj.Number_T);
% add a long T value
Times(Option_obj.Number_T + 1) = 1.5;
disp('Finished Calculating!')
toc
end
