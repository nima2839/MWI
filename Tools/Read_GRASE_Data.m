function Read_GRASE_Data(Read_Path,name,Save_Path)
mgrase = ReadAllDCM(Read_Path,32);
cd(Save_Path)
save(name)
end