function Dist = Gaussian_Dist(x, Params)
  % Create a pdf distribution
  Mu = Params(1);
  Sigma = Params(2);
  Dist = pdf(makedist('Normal','mu',Mu,'sigma',Sigma),x);
end
