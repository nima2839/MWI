function Dist = MixtureModel(x, Amps, Params, Models)
  % Amps: a vector containing relative amplitude of each pdf
  % Params is a cell containing parmeters for the model
  % Model: a cell of different Models
  Dist = 0;
  n = length(Params);
  for i = 1:n
    Dist = Dist + Amps(i) * Models{i}(x, Params{i});
  end
end
