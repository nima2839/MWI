function MWF_Map = Find_MWF(Data, Index , Method)
%
%   Author: Nima
%   Date: 02/2019
% Inputs:
%   Data: 4D distribtion / Params
%   Index: determines last data to be considered as myelin
    sd = size(Data);
    Data = double(Data);
    MWF_Map = zeros(sd(1),sd(2),sd(3));
    disp('MWF map calculation started...')
   if strcmp(Method,'NNLS')
    for k = 1:sd(3)
        for j = 1:sd(2)
            parfor i = 1:sd(1)
                temp = squeeze(Data(i,j,k,:));
		if sum(temp(:)) > 0
                	MWF_Map(i,j,k) = sum(temp(1:Index))/sum(temp(:));
            	end
	    end
        end
    end
   elseif strcmp(Method,'NLLS')
      for k = 1:sd(3)
        for j = 1:sd(2)
            parfor i = 1:sd(1)
                temp = squeeze(Data(i,j,k,:));
                MWF_Map(i,j,k) = temp(Index(1))/(temp(Index(1)) + temp(Index(2)) + temp(Index(3)));
            end
        end
      end
	elseif strcmp(Method,'NLLS_Test')
		np = sd(1);
		nv = sd(2);
		ns = sd(3);
% 		parfor k = 1:ns
%         		for j = 1:nv
%             			for i = 1:np
%                                 temp = squeeze(Data(i,j,k,:));
                                MWF_Map(:,:,:) = Data(:,:,:,Index(1)).*((Data(:,:,:,Index(1)) + Data(:,:,:,Index(2)) + Data(:,:,:,Index(3))).^-1);

%             			end
%         		end
%       	end
   else
       error('Specify method!')
   end
    disp('Finished Calculating!')
end
