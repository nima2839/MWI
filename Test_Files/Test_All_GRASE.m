addpath('~/GRASE/Postprocessing')
addpath(genpath('~/MWI'))
names{1} = 'Can-05-Rrm-064-M0R';

for i = 1:length(names)
  ProcessGRASE(names{i});
end
