function ProcessGRE(FileName,ReadPath,SavePath)

cd(ReadPath)
disp(['Reading ', FileName,]);
load(FileName)

tic
test = TestClass(Mag,Phase,Info);
clear Mag Phase Info
test = Normal_Procedure(test);
disp('Saving results...')
test.RunTime = toc;
GRE_MWF_data = GetAllData(test);
clear test
cd(SavePath)
save(['C3PM_Results' , FileName], 'GRE_MWF_data');
end
