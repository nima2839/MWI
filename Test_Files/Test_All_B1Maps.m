addpath(genpath('~/GRASE/Postprocessing'))
addpath(genpath('~/MWI'))

names{1} = 'Can-05-Rrm-005-M0';
names{2} = 'Can-05-Rrm-011-M0';
names{3} = 'Can-05-Rrm-013-M0';
names{4} = 'Can-05-Rrm-082-M0R';
names{5} = 'Can-05-Rrm-019-M0';
names{6} = 'Can-05-Rrm-023-M0';
names{7} = 'Can-05-Rrm-025-M0';
names{8} = 'Can-05-Rrm-031-M0';
names{9} = 'Can-05-Rrm-043-M0';
names{10} = 'Can-05-Rrm-053-M0';
names{11} = 'Can-05-Rrm-054-M0';
names{12} = 'Can-05-Rrm-057-M0';
names{13} = 'Can-05-Rrm-064-M0R';
names{14} = 'Can-05-Rrm-065-M0R';
names{15} = 'Can-05-Rrm-070-M0R';
%names{16} = 'Can-05-Rrm-018-M0';

for i = 1:length(names)
  UBC_Vs_BSB1(names{i});
end
