function [T2, B1, res, Options] = Single_Component_T2(Signal, Options)
% Fits T2  using single component fitting and EPG algorithm
% Here B1 is the refocusing angle
% 
	if nargin < 2
		Options.Empty =  true;
	end
	
	if ~isfield(Options, 'T2Range')
		Options.T2Range = [10e-3, 2]; % seconds
	end
	
	if ~isfield(Options, 'T1')
		Options.T1 = 1; %seconds
	end
	
	if ~isfield(Options, 'TE')
		Options.TE = 10e-3;
	end
	
	if ~isfield(Options, 'ETL')
		Options.ETL = length(Signal);
	end
	
	if ~isfield(Options, 'B1')
		[T2, B1, res, Options] = Single_Component_T2_B1(Signal, Options);
		return;
	else
		B1 = Options.B1;
	end
	
	opt = optimoptions('lsqnonlin','Algorithm', 'trust-region-reflective','TolFun',1e-12,'MaxIter',1e4,'TolX',1e-8,'Display','off');
	opt.MaxFunEvals = 1e4;
	[Param, res] = lsqnonlin(@CostFun ,...
		[60e-2 , Signal(1) * 2],... 			% Initial guess
		[Options.T2Range(1), Signal(1)/2],...	% Lower boundary
		[Options.T2Range(2), inf],...			% Upper boundary
		opt);
	T2 = Param(1);
	Options.Proton_Density = Param(2);
	function out = CostFun(X)
		T2 = X(1);
		Fit = EPGdecaycurve(Options.ETL,...
			B1,...
			Options.TE,...
			T2,...
			Options.T1,...
			180);
		out = Signal - X(2)*Fit;
	end
	
end