function [T, T_weights, resnorm] = NNLS_T_Fitting(Data, EchoTime_obj, Option_obj)
% Non-linear least squres TimeConstant fitting
%   Author: Nima
%   Date: 07/2019
% Inputs:
%   Data: Contains the signal value during different echoes
%   EchoTime_obj.FirstEchoTime: time of the first echo
%   EchoTime_obj.EchoSpacing: spacing time difference between echoes
%   Option_obj.T_range: Lower and upper bounds for T values
%   Option_obj.Number_T: Number of T values
%   Option_obj.A: Placeholder for solution space (enter empty array)
%   Option_obj.Mu: Regularization factor (Currently calulated in this func)
%
% Outouts:
%   T: Decay Times used in the fitting
%   T_weights: Weights associated with decay times calculated in fitting
%
Gamma = 42.575e6;
    if nargin < 2
        error('Function requires at least two inputs!');
    end

    if isempty(EchoTime_obj)
        error('Echotime_obj can not be empty!');
    end


    if nargin < 3 || isempty(Option_obj)
        % setting default values for Option_obj
        Option_obj.T_range = [5e-3 120e-3];
        Option_obj.Number_T = 100;
        Option_obj.A = [];
    end
    temp = size(Data);
    if temp(1) == 1
        Data = Data';
    end

    T = linspace(Option_obj.T_range(1),Option_obj.T_range(2),...
        Option_obj.Number_T);
    % add a long T value
    Option_obj.Number_T = Option_obj.Number_T + 1;
    T(Option_obj.Number_T) = 1;

    if max(Data(:)) == 0
      resnorm = 0;
      T_weights = T*0;
      return;
    else
      Data = Data/max(Data(:));
    end

    % Setting Echo Times
    Number_echos = length(Data);
    EchoTimes = zeros(Number_echos,1);
    EchoTimes(1) = EchoTime_obj.FirstEchoTime;
    for k = 2:Number_echos
       EchoTimes(k) = EchoTimes(1) + (k-1)*EchoTime_obj.EchoSpacing ;
    end

    if ~isfield(Option_obj,'A') || isempty(Option_obj.A)
        sT = size(T);
        temp = zeros(Number_echos,max(sT));
        for i = 1:(max(sT) -1)
            temp(:,i) = transpose(exp(-EchoTimes/T(i)));
        end
        temp(:,max(sT)) = transpose(exp(-EchoTimes/T(max(sT))));
        Option_obj.A = temp;
    end

    Reg =  eye(Option_obj.Number_T);
    optmp = optimset('fminbnd');
    optmp.TolX = eps;
    optmp.MaxFunEvals = 800;
    optmp.NaxIter = 800;
    Mu = fminbnd(@NNLSreg_obj,0,1,optmp);
%     Mu = Option_obj.Mu;
    [T_weights,resnorm] = lsqnonneg([Option_obj.A; Mu*Reg],...
        [Data; zeros(Option_obj.Number_T,1)]);

    function f = NNLSreg_obj(x)
        [~,f] = lsqnonneg([Option_obj.A;x*Reg], ...
            [Data;zeros(Option_obj.Number_T,1)]);
    end
end
