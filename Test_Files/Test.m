names = ["123","127","Jh","Js","Nm","Nn","Rws","Zz","Rh"];
nangles = [158.0052  147.3529  144.4640  152.8400  139.5369 149.6512  148.7811  145.5094  143.7493];



for i = 1:length(names)
	cd ~/GRASE/GRASE_To_Do/
	file = dir(strcat('*', names(i),'*'));
	UBC_Vs_BSB1(file.name, nangles(i));
	cd ~/GRASE/GRASE_Results/GRASE_B1_Map_Results/
	file = dir(strcat('*', names(i),'*'));
	load(file.name,'maps','MWI')
	cd InVivo_Results/
	mkdir(names(i))
	cd(names(i))
	save('Maps_B1', 'maps')
	niftiwrite(MWI, 'MWF_B1')
end