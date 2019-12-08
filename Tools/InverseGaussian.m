function Dist = InverseGaussian(x, Params)
  % Create a pdf distribution
  Mu = Params(1);
  Sigma = Params(2);
  if Sigma == 0
    Sigma = 1e-5;
  end
  Lamda = (Mu^3)/Sigma;
  Dist = pdf(makedist('InverseGaussian','mu',Mu,'lambda',Lamda),x);
end
