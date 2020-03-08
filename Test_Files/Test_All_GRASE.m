addpath(genpath('~/GRASE/Postprocessing'))
addpath(genpath('~/MWI'))

%cd ~/GRASE/GRASE_To_Do/
%files = dir('GRASE_Results_*.mat');

name{1} = 'Can-05-Ppm-067-M0R';
name{2} = 'Can-05-Ppm-076-M0R';
name{3} = 'Can-05-Ppm-078-M0R';
name{4} = 'Can-05-Rrm-087-M0R';
name{5} = 'Can-05-Ppm-088-M0R';
name{6} = 'Can-05-Rrm-089-M0R';
name{7} = 'Can-05-Ris-095-M0R';
name{8} = 'Can-05-Rrm-086-M0R';
name{9} = 'Can-05-Rrm-094-M0R';
name{10} = 'Can-05-Ris-105-M0R';
name{11} = 'Can-05-Ppm-107-M0R';
name{12} = 'Can-05-Rrm-113-M0R';

for i = 1:length(name)
  %ProcessGRASE(files(i).name);
  ProcessGRASE(name{i});
end
