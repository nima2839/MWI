function [Param, res] = Sinc3PM_NLLS(Data,Info,X0,lb,ub)

if ~isfield(Info, 'Algorithm')
	Info.Algorithm = 'trust-region-reflective';
end
if ~isfield(Info, 'MultiStart')
	Info.MultiStart = false;
elseif ~isfield(Info, 'NumMS') %Number of MultiStart points
	Info.NumMS = 20;
end

t = (Info.FirstTE):(Info.EchoSpacing):(Info.FirstTE + (length(Data)-1) * Info.EchoSpacing);

if size(t) ~= size(Data)
    t= t';
end
try
	if max(Data(:)) ~= 0
		ND = Data/max(Data(:));
		Model = @(X)abs((X(1)*	exp(-t.*(X(2)+ 1i*2*pi*X(3))) + ...
		                 X(4)*	exp(-t.* X(5)) + ...
		                (X(6))*	exp(-t.*(X(7)+ 1i*2*pi*X(8))))).*sinc(X(9)  * t) -ND;
		if nargin < 3
			X0 = [0.1,   60,	  0,	0.6,	25,	0.3,	15,	   0,  0];
			lb = [0,     50,	 -70,	0,	  5,	0,	  0.1,	-30,  0];
			ub = [2,	  1000,	70,	2,	  50,	2,	  50,	   30,  20];
		end

    if X0(9) > 0
      ub(9) = 3*X0(9);
    end
		options = optimoptions('lsqnonlin','Algorithm',Info.Algorithm,'TolFun',1e-12,'MaxIter',1e3,'TolX',1e-8,'Display','off');
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
	disp('An exception was caught in the ThreePoolModel_NLLS fitting:')
    	disp(ME.message)
	res = nan;
	Param = X0.*nan;
end
end
