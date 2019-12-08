function [R, res] = SingleComponentNLLS(Data,Info,X0,lb,ub)

if ~isfield(Info, 'Algorithm')
	Info.Algorithm = 'trust-region-reflective';
end

t = (Info.FirstTE):(Info.EchoSpacing):(Info.FirstTE + (length(Data)-1) * Info.EchoSpacing);

if size(t) ~= size(Data)
    t= t';
end
try
	if max(Data(:)) ~= 0
		ND = Data/max(Data(:));
		Model = @(X)(X(1)*exp(-t.*X(2)) + X(3)) -ND;
		if nargin < 3              
		    X0 = [1	20	0];
		    lb = [0	0	0];
		    ub = [1.5	100	0.5];
		end
		Params = zeros(3);		
		options = optimoptions('lsqnonlin','Algorithm',Info.Algorithm,'TolFun',1e-12,'MaxIter',2e4,'TolX',1e-7,'Display','off');
		options.MaxFunEvals = 2e4;
	
		if strcmp(Info.Algorithm,'levenberg-marquardt')
			[Param, res] = lsqnonlin(Model,X0,[],[],options);
		elseif strcmp(Info.Algorithm,'trust-region-reflective')
			[Param, res] = lsqnonlin(Model,X0,lb,ub,options);
		else
			error('Invalid Algorithm input!');
		end
		R = Param(2);
	else
		res = 0;
		R = 0;
	end
catch ME
	disp('An exceptin was caught in the SingleComponentNLLS fitting!')
	res = nan;
	R = nan;
end
end
