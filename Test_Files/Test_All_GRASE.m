addpath(pwd)
cd ~/Postprocessing
addpath(genpath(pwd))
cd ~
addpath MWI/
names{1} = 'Can-05-Rrm-064-M0R';
names{2} = 'Can-05-Rrm-065-M0R';
names{3} = 'Can-05-Rrm-070-M0R';
names{4} = 'Can-05-Rrm-082-M0R';

for i = 1:length(names)
  ProcessGRASE(names{i});
end
