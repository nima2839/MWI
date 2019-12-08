function [Param, res] = TwoPoolModel_NLLS(Data,Info,X0,lb,ub)

if ~isfield(Info, 'Algorithm')
	Info.Algorithm = 'trust-region-reflective';
end
if ~isfield(Info, 'MultiStart')
	Info.MultiStart = false;
elseif ~isfield(Info, 'NumMS') %Number of MultiStart points
	Info.NumMS = 25;
end

t = (Info.FirstTE):(Info.EchoSpacing):(Info.FirstTE + (length(Data)-1) * Info.EchoSpacing);

if size(t) ~= size(Data)
    t= t';
end
try
	if max(Data(:)) ~= 0
		ND = Data/max(Data(:));
		Model = @(X)abs(X(1)*	exp(-t.*(X(2)+ 1i*2*pi*X(3))) + ...
		                X(4)*	exp(-t.* X(5))) -ND; % R value should not be over
		if nargin < 3              
			X0 = [0.1,	100,	10,	0.9, 25];
	   		lb = [0,	50,     -40,	0.5, 0];
	   		ub = [0.5,	1000,	40,	1, 50];
		end
		options = optimoptions('lsqnonlin','Algorithm',Info.Algorithm,'TolFun',1e-12,'MaxIter',1e3,'TolX',1e-7,'Display','off');
		options.MaxFunEvals = 1e3;
	
		if strcmp(Info.Algorithm,'levenberg-marquardt')
			if Info.MultiStart
                		Problem = createOptimProblem('lsqnonlin','x0',X0,'objective',Model,'options',options);
				ms = MultiStart('Display','off');
				[Param, res] = run(ms, Problem, Info.NumMS);
			else
				[Param, res] = lsqnonlin(Model,X0,[],[],options);
			end
		elseif strcmp(Info.Algorithm,'trust-region-reflective')
			if Info.MultiStart
				Problem = createOptimProblem('lsqnonlin','x0',X0,'objective',Model,'options',options,...
					'lb',lb,'ub',ub);
				ms = MultiStart('Display','off');
				[Param, res] = run(ms, Problem, Info.NumMS);
			else
				[Param, res] = lsqnonlin(Model,X0,lb,ub,options);
			end
		else
			error('Invalid Algorithm input!');
		end
	else
		res = 0;
		Param = X0.*0;
	end
catch ME
	disp('An exception was caught in the TwoPoolModel_NLLS fitting:')
	disp(ME.message)
	res = nan;
	Param = X0.*nan;
end
end
