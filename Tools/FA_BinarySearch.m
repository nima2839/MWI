function output = FA_BinarySearch(decay_data, TE, T1, Center , radius)
  % Finds flip angle using binary search algorithm
  tic
  ETL = length(decay_data);
  nT2 = 60;
  T2_times=logspace(log10(0.01),log10(2),nT2);
  Num_Steps = 9;
  if nargin < 4
    radius = 32;
    Center = 180 - 2 * radius;
  end

  for i = 1:Num_Steps
    [~,r1,~] = lsqnonneg(Calc_basis_decay(ETL, nT2, Center + radius, TE, T2_times, T1, 180), decay_data');
    [~,r2,~] = lsqnonneg(Calc_basis_decay(ETL, nT2, Center - radius, TE, T2_times, T1, 180), decay_data');
    if r1 > r2
      Center = Center - radius;
    else
      Center = Center + radius;
    end
    radius = radius / 2;
  end
  output = Center;
  toc
end

function basis_decay = Calc_basis_decay(nechs, nT2, alpha, TE, T2_times, T1, RefCon)
  basis_decay=zeros(nechs,nT2);
  % Compute the NNLS basis over T2 space
  if alpha == 0 % this prevents the code from freezing
	alpha = 180;
  end
  for x=1:nT2
      echo_amp = EPGdecaycurve(nechs, alpha, TE, T2_times(x), T1, RefCon);
      basis_decay(:,x) = echo_amp';
  end
end
