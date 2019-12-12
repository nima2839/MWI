function [Param, res] = VP_3PM(Data,Info,X0,lb,ub)
  if ~isfield(Info, 'Algorithm')
  	Info.Algorithm = 'trust-region-reflective';
  end

  t = (Info.FirstTE):(Info.EchoSpacing):(Info.FirstTE + (length(Data)-1) * Info.EchoSpacing);

  if size(t) ~= size(Data)
      t= t';
  end
%   try
  	if max(Data(:)) ~= 0
  		ND = Data/max(Data(:));
  		%Model = @(X)abs((X(1)*	exp(-t.*(X(2))) + ...
  		%                 X(3)*	exp(-t.* X(4)) + ...
  		%                (X(5))*	exp(-t.*(X(6))))) -ND;
  		if nargin < 3
  			X0 = [60,	20,	10];
  			lb = [35,	5,	0.1];
  			ub = [1000,	50,	50];
  		end
  		options = optimoptions('lsqnonlin','Algorithm',Info.Algorithm,'TolFun',1e-12,'MaxIter',1e3,'TolX',1e-8,'Display','off');
  		options.MaxFunEvals = 1e3;

  		if strcmp(Info.Algorithm,'levenberg-marquardt')
  			[Param, res] = lsqnonlin(@Model,X0,[],[],options);
  		elseif strcmp(Info.Algorithm,'trust-region-reflective')
  			[Param, res] = lsqnonlin(@Model,X0,lb,ub,options);
  		else
  			error('Invalid Algorithm input!');
  		end
  	else
  		res = 0;
  		Param = zeros(1,6);
  	end
%   catch ME
%   	disp('An exception was caught in the Mag_ThreePoolModel_NLLS fitting:')
%       	disp(ME.message)
%   	res = nan;
%   	Param = zeros(1,6).*nan;
%   end
    Basis_Curves(1,:) = exp(-Param(1)*t);
    Basis_Curves(2,:) = exp(-Param(2)*t);
    Basis_Curves(3,:) = exp(-Param(3)*t);
    [Param(4:6),res] = lsqnonneg(Basis_Curves', ND');
  function out = Model(X)
    Basis_Curves(1,:) = exp(-X(1)*t);
    Basis_Curves(2,:) = exp(-X(2)*t);
    Basis_Curves(3,:) = exp(-X(3)*t);
    [~,~,out] = lsqnonneg(Basis_Curves', ND');
  end
end
