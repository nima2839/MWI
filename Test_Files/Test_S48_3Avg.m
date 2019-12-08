tic
Mag = ReadAllDCM('/media/work/aelkady/3T/Ahmed_96Echo_3Dtest_2/study/gre_96contrasts_3avg_7',48);
Phase = ReadAllDCM('/media/work/aelkady/3T/Ahmed_96Echo_3Dtest_2/study/gre_96contrasts_3avg_8',48);

Phase = double(Phase);
Phase = 2*pi* Phase/max(abs(Phase(:))) - pi;
Mag = double(Mag);

Info.FirstTE = 1.45e-3;
Info.EchoSpacing = 1.05e-3;
Info.Vox = [2 2 3] * 1e-3;
[Gp, Gv, Gs] = LFG_Correction(Mag, Phase, Info);

Opt.T_range = [1e-3 0.120];
Opt.Number_T = 120;
Opt.Vox = Info.Vox;
LFG.Gp = Gp;
LFG.Gs = Gs;
LFG.Gv = Gv;
E.FirstEchoTime = 1.45e-3;
E.EchoSpacing = 1.05e-3;

Fitted = NNLS_T_Batch_Fitting(Mag,E,Opt);
LFGC_Fitted = NNLS_T_Batch_Fitting(Mag,E,Opt,LFG);

MWF = Find_MWF(Fitted, 25);
LFGC_MWF = Find_MWF(LFGC_Fitted, 25);

RunT = toc;
disp('Saving results...')
% This is for linux matlab
cd ~
save('S48_3Avg')
disp('Done!')
