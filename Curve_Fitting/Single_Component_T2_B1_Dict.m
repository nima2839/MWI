function [T2, B1, Residual] = Single_Component_T2_B1_Dict(Signal, Dict, T2Range, B1Range)
% Fits for a parameter using single component fitting using a dictionary matrix and interpolates the range for T2
% Dict: the dictionary matrix which must be a 2D matrix of decay curves (MxNxETL)
% B1Range (Nx1)
% T2Range (Mx1)
% Signal(ETLx1)
% Range: range of the parameter used to generate the dictionary matrix 
	if nargin < 4
		err('Not enough input arguments!');
	end
	
	%Changing the order of dictionary matrix
	Dict = permute(Dict, [3,1,2]);
	res = zeros(1,length(B1Range));
	T2Vals = res;
	for i = 1:length(B1Range)
		temp_dict =  squeeze(Dict(:,:,i));
		[T2Vals(i),res(i)] = Single_Component_T2_Dict(Signal, temp_dict, T2Range , 2); % Using NNLS for fitting to speed up the process!
	end
	[~,index] = min(res);
	B1 = B1Range(index);
	T2 = T2Vals(index);
	Residual = min(res);
end