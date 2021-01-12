function Out = Single_Component_T2_Dict(Signal, Dict, Range)
% Fits for a parameter using single component fitting using a dictionary matrix and interpolates the range
% Dict: the dictionary matrix which must be a 2D matrix of decay curves (ETLxN)
% Signal(ETLx1)
% Range: range of the parameter used to generate the dictionary matrix 
	if nargin < 3
		err('Not enough input arguments!');
	end
	
	ETL = length(Signal);
	% Rshaping matrices for lsqnonneg function
	Signal = reshape(Signal, [ETL,1]);
	res = zeros(1,length(Range));
	for i = 1:length(Range)
		[~,res(i),~] = lsqnonneg(Dict, Signal);
	end
	Range_spline = logspace(log10(Range(1)), log10(Range(end), 1e3));
	res_spline = interp1(T2Range, res, Range_spline,'spline');
	[~,index] = min(res_spline);
	Out = Range_spline(index);
end