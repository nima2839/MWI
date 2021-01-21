% This simulates an MEGRE phantom
tic;
disp('Creating 3D phantom...')
N = 128;
p = phantom3d('Modified Shepp-Logan', N);
p(p==0) =nan;
p(p > 0.25 & p < 0.35) = 0.10;
p(p > 0.16 & p < 0.25) = 0.05;
MWF = abs(p(:,:, (floor(N/2) - 15): (floor(N/2) +16)));

temp = zeros(size(MWF));
temp(MWF == .1) = 15;
temp(MWF == .05) = -5;
opt.FreqShift.MW = temp;
temp(MWF == .1) = 5;
temp(MWF == .05) = 1;
opt.FreqShift.IE = temp;
opt.Times = (1:32) * 2e-3;

disp('Generating 4D MEGRE Phantom...')
SNR = [100 200 500 1e3];
Phantoms = cell(1,length(SNR));
data = cell(1,length(SNR));

MyInfo.Vox = [1 1 4] * 1e-3;
MyInfo.Mask = zeros(size(MWF));
MyInfo.Mask(MWF > 0) = 1;
MyInfo.FirstTE = opt.Times(1);
MyInfo.EchoSpacing = opt.Times(2) - opt.Times(1);

for i = 1:length(SNR)
	opt.SNR = SNR(i);
	Phantoms{i} = Create_MEGRE_Phantom_4D(MWF, opt);
	test = TestClass(abs(Phantoms{i}), angle(Phantoms{i}), MyInfo);
	test = Calc_SC(test,2); % LOG method
	test = Calc_3PM(test);
	data{i} = GetAllData(test);
end
RunTime = toc;

clear test p MyInfo N

cd ~/Simulation/Phantom/
save('Phantom_MEGRE_Sim','-v7.3')