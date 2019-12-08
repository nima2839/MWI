function [Param, res] = Complex3PM(Data,Info,X0,lb,ub)

Info.Algorithm = 'trust-region-reflective';

if ~isfield(Info,'Phi0')
	%error('Phi0 is missing!');
	Info.Phi0 = 0;
end



t = (Info.FirstTE):(Info.EchoSpacing):(Info.FirstTE + (length(Data)-1) * Info.EchoSpacing);

if size(t) ~= size(Data)
    t= t';
end
try
	if max(Data(:)) ~= 0
		ND = Data/max(abs(Data(:)));
		if nargin < 3
			X0 = [0.1,   60,	  5,	0.6,	25, 0,	0.3,	15,	   0];
			lb = [0,     50,	 -70,	0,	  5,-25,	0,	  0.1,	-30];
			ub = [2,	  1000,	70,	2,	  50,25,	2,	  50,	   30];
		end
		options = optimoptions('lsqcurvefit','Algorithm',Info.Algorithm,'TolFun',1e-12,'MaxIter',1e4,'TolX',1e-8,'Display','off');
		options.MaxFunEvals = 1e4;
        ydata(:,1) = real(ND);
		ydata(:,2) = imag(ND);
		[Param, res] = lsqcurvefit(@cplxreal,X0,t,ydata,lb,ub,options);
	else
		res = 0;
		Param = X0.*0;
	end
catch ME
	disp('An exception was caught in the Complex3PM fitting:')
    disp(ME.message)
	res = nan;
	Param = X0.*nan;
end

	function yout = cplxreal(X,t)
		Model = (X(1)*	exp(-t.*(X(2)- 1i*2*pi*X(3))) + ...
		            X(4)*	exp(-t.* (X(5) - 1i*2*pi*X(6))) + ...
		            (X(7))*	exp(-t.*(X(8) - 1i*2*pi*X(9)))).*exp(1i*(Info.Phi0));
		yout(:,1) = real(Model);
		yout(:,2) = imag(Model);
	end
end
