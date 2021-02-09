function [Out, Residual] = Single_Component_T2_Dict(Signal, Dict, Range)
% Nima: This is currently too slow to use!
% Fits for a parameter using single component fitting using a dictionary matrix and interpolates the range
% Dict: the dictionary matrix which must be a 2D matrix of decay curves (ETLxN)
% Signal(ETLx1)
% Range: range of the parameter used to generate the dictionary matrix 
	if nargin < 3
		err('Not enough input arguments!');
	end
	
	ETL = length(Signal);
	% Reshaping matrices for fitting function
	Signal = reshape(Signal, [ETL,1]);
	res = zeros(1,length(Range));
	% setting options for lsqnonlin function
	options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','TolFun',1e-4,'MaxIter',400,'TolX',1e-4,'Display','off');
	options.MaxFunEvals = 200;
	for i = 1:length(Range)
		Model = @(X)(X* Dict(:,i)) -Signal;
		[~, res(i)] = lsqnonlin(Model,Signal(1),[],[],options);;
	end
	Range_spline = interp1(Range, 1:.1:numel(Range));
	res_spline = interp1(Range, res, Range_spline,'spline');
	[~,index] = min(res_spline);
	Out = Range_spline(index);
	Residual = min(res);
end