function [Beta, res,pred] = MOIG_Fitting(Signal,Info,X0)
  % Non-linear least squares fitting
  % Info.ETL : echo train length
  % Info.TE : enter-echo EchoSpacing
  % Info.FA : Flip Angle
  % Info.T2 : T2 Values used in the fiiting process

  tic
  Signal = double(Signal)/max(Signal(:));

  if nargin < 3
    % [A1 T21 sigma1 A2 T22 sigma2 A3 T23 sigma3] sigma = standard deviation
    X0 = [0.1 25e-3 1e-3, 0.8 8e-2 1e-3, 0.1 0.5 1e-2];
    if nargin < 2
      Info.ETL = 32;
    end
  end

  if ~isfield(Info,'ETL')
  	Info.ETL = 32;
  end

  if ~isfield(Info,'TE')
  	Info.TE = 1e-2;
  end

  if ~isfield(Info,'FA')
  	Info.FA = 180;
  end

  if ~isfield(Info,'T2')
  	Info.T2 = [1:199,200:10:999, 1e3:1e2:2e3] * 1e-3;
  end

  lb = [0 1e-3 0, 0.5 4e-2 0, 0 0.1 0];
  ub = [1 35e-3 1e-2, 2 12e-2 2e-2, 2 2 1e-1];

  Times = (1:Info.ETL) * Info.TE;
  LookUpTable = zeros(length(Info.T2), Info.ETL);

  for i = 1:length(Info.T2)
    LookUpTable(i,:) = EPGdecaycurve(Info.ETL, Info.FA, Info.TE, Info.T2(i), 1, 180);
  end

  Weights = exp(Signal).^2;
  Model = @(b) Weights.*(ModelFunc(b,Times) - Signal);
  Info.Algorithm = 'trust-region-reflective';
  options = optimoptions('lsqnonlin','Algorithm',Info.Algorithm,'TolFun',1e-12,'MaxIter',1e4,'TolX',1e-12,'Display','off');
  options.MaxFunEvals = 1e4;
  [Beta, res] = lsqnonlin(Model,X0,lb,ub,options)
  pred = ModelFunc(Beta,Times);
  toc

  function out = ModelFunc(b,x)
    Params{1} = [b(2) b(3)];
    Amps(1) = b(1);
    Params{2} = [b(5) b(6)];
    Amps(2) = b(4);
    Params{3} = [b(8) b(9)];
    Amps(3) = b(7);
    Models{1} = @InverseGaussian;
    Models{2} = @Gaussian_Dist;
    Models{3} = @Gaussian_Dist;
    out = MixtureModel(Info.T2, Amps, Params, Models) * LookUpTable;
    out = out/max(out);
  end
end
