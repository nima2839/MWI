Mag =   ReadAllDCM('~/Nm_Myelin_Test_1/study/gre_96contrasts_new_10',48);
Phase = ReadAllDCM('~/Nm_Myelin_Test_1/study/gre_96contrasts_new_11',48);

Phase = double(Phase);
Phase = 2*pi* Phase/4095 - pi;
Mag = double(Mag) / 4095;
cd ~
save('NM_Mag_Phase')
disp('Done!')
