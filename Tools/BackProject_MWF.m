function output = Backproject_MWF(decay_data, TE, FA)
  % Get MESE data and uses the tail to calculate MWF
  T1 = 1;
  ETL = length(decay_data);
  nT2 = 60;
  T2_times=logspace(log10(4e-2),log10(2),nT2);
  obs_weigts = ones(size(decay_data));
  Chi2Factor = 1.2;
  RefCon = 180;
  basis_decay = Calc_basis_decay(ETL, nT2, FA, TE, T2_times, T1, RefCon);
  index = (ETL/2):ETL;

  [T2_dis,mu,chi2] = Nima_UBC_NNLS(basis_decay(index,:), decay_data(index)', obs_weigts(index), 1.02);

  decay_pred = basis_decay * T2_dis;
  diff = (decay_data(:) - decay_pred(:));

  T2_times=logspace(log10(1e-2),log10(4e-2),nT2);
  [T2_dis2,mu,chi2] = Nima_UBC_NNLS(basis_decay, diff, obs_weigts, 1.02);
  output = sum(T2_dis2(:))/sum(T2_dis(:));
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
