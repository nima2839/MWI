function [Param, res] = ThreePoolM_NLLS(Data,Info,X0,lb,ub)

if ~isfield(Info, 'Algorithm')
	Info.Algorithm = 'trust-region-reflective';
end
if ~isfield(Info, 'MultiStart')
	Info.MultiStart = false;
elseif ~isfield(Info, 'NumMS') %Number of MultiStart points
	Info.NumMS = 100;
end

t = (Info.FirstTE):(Info.EchoSpacing):(Info.FirstTE + (length(Data)-1) * Info.EchoSpacing);

if size(t) ~= size(Data)
    t= t';
end
try
	if max(Data(1)) ~= 0
		ND = Data/Data(1); % Normalize signal
		Model = @(X)abs((X(1)*	exp(-t.*((1/X(2))+ 1i*2*pi*X(3)))) + ...
		                 (X(4)*	exp(-t.* (1/X(5)))) + ...
		                (X(6)*	exp(-t.*((1/X(7))+ 1i*2*pi*X(8))))) -ND;
		if nargin < 3
			X0 = [0.1,   10e-3,	0,		0.6,	64e-3,	0.3,	48e-3,		0];
			lb = [0,     3e-3,	-25,	0,		25e-3,	0,		25e-3,		-10];
			ub = [2,	 25e-3,	25,		2,		150e-3,	2,		150e-3,		10];
		end
		options = optimoptions('lsqnonlin','Algorithm',Info.Algorithm,'TolFun',1e-7,'MaxIter',1e3,'TolX',1e-7,'Display','off'); % To include Jacobian Calculations must add: 'SpecifyObjectiveGradient',true,
		options.MaxFunEvals = 1e3;
		options.DiffMinChange = 1e-3;
		options.TypicalX = [0.1,   10e-3,	5,		0.6,	64e-3,	0.3,	48e-3,		1]; %for scaling finite differences for gradient estimation
		

		if strcmp(Info.Algorithm,'levenberg-marquardt')
			if Info.MultiStart
                		Problem = createOptimProblem('lsqnonlin','x0',X0,'objective',@MyFunc,'options',options);
				ms = MultiStart('Display','off');
				[Param, res] = run(ms, Problem, Info.NumMS);
			else
				[Param, res] = lsqnonlin(@MyFunc,X0,[],[],options);
			end
		elseif strcmp(Info.Algorithm,'trust-region-reflective')
			if Info.MultiStart
				Problem = createOptimProblem('lsqnonlin','x0',X0,'objective',@MyFunc,'options',options,...
					'lb',lb,'ub',ub);
				ms = MultiStart('Display','off');
				[Param, res] = run(ms, Problem, Info.NumMS);
			else
				[Param, res] = lsqnonlin(@MyFunc,X0,lb,ub,options);
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

	function F = MyFunc(X)
		F = Model(X);
		% Following is for Jacobian calculations:
		%J1 = exp(-t.*((1/X(2))+ 1i*2*pi*X(3)));
		%J2 = (1/X(2)^2) * t.*  J1;
		%J3 = -1i*2*pi* t.* J1;
		%J4 = exp(-t.* (1/X(5)));
		%J5 = (1/X(5)^2) * t.*  J4;
		%J6 = exp(-t.*((1/X(7))+ 1i*2*pi*X(8)));
		%J7 = (1/X(7)^2) * t.*  J6;
		%J8 = -1i*2*pi* t.* J6;
		%J = [J1', J2', J3', J4', J5', J6', J7', J8'];
	end
end
