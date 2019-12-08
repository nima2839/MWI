classdef ADF
% Source: Nonlinear Anisotropic Filtering of MRI Data by Greig et. al.
	properties
		K % This parameter is used to calculate the diffusion func.
		Mag
		SizeData
		EchoIndexes
		Kernel
		Vox
		IntegConst % Integration constant delta
		NumIter
	end
	methods
		function obj=ADF(MaskedMag,Vox,NumIter,K)
			obj.Mag = MaskedMag;
			obj.SizeData = size(MaskedMag);
			obj.Vox = Vox;
 			obj.NumIter = NumIter;
			obj.K = K;
			obj.EchoIndexes = 1:obj.SizeData(4);
			%obj.Kernel = Calc_Kernel(Vox);
			obj.IntegConst = (1/7);	% (1/7) was acquired from the article
		end

		%function k=Calc_Kernel(Vox)
			% It is assumed that Vox(1) = Vox(2)
		%	k = zeros(3,3,3);
			%a = 1/Vox(1);
			%b = 1/(sqrt(2)*Vox(1));
			%c = 1/sqrt(Vox(1)^2 + Vox(3)^2);
			%d = 1/sqrt(2*Vox(1)^2 + Vox(3)^2);
			%f = 1/Vox(3);
			%k(:,:,1) = -[d,c,d;c,f,c;d,c,d];
			%k(:,:,2) = [-b,-a,-b;-a,0,a;b,a,b];
			%k(:,:,3) = [d,c,d;c,f,c;d,c,d];
		%end


		function [fv,fp,fs] = CalcFlowFunc(obj,Inten)
			[gv,gp,gs] = obj.CalcGradient(Inten,obj.Vox);
			diffusionFunc = obj.CalcDiffusionFunc(gv,gp,gs,obj.K);
			fp =  gp.*diffusionFunc;
			fv =  gv.*diffusionFunc;
			fs =  gs.*diffusionFunc;
		end

		function Output=ApplyADF(obj)
			disp('ADF Started!')
			IntenN = obj.Mag;
			vox = obj.Vox;
			ei = obj.EchoIndexes;
			dt = obj.IntegConst;
			nv = obj.SizeData(1);
			np = obj.SizeData(2);
			ns = obj.SizeData(3);
			%Kernel = obj.Kernel;
			for i = 1:obj.NumIter
				fprintf('Iteration Number: %d ...',i)
				[FNv,FNp,FNs] = CalcFlowFunc(obj,IntenN); % Flow function in the nth iteration
					temp1 = zeros(size(FNv));
					temp2 = temp1;
					% V-direction
					temp1(2:nv,:,:,:) = FNv(1:nv-1,:,:,:);
					temp2(1:nv-1,:,:,:) = FNv(2:nv,:,:,:);
					temp = (temp2-temp1)/(2*vox(1));
					% P-direction
					temp1(:,2:np,:,:) = FNp(:,1:np-1,:,:);
					temp2(:,1:np-1,:,:) = FNp(:,2:np,:,:);
					temp = temp + (temp2-temp1)/(2*vox(2));
					% S-direction
					temp1(:,:,2:ns,:) = FNs(:,:,1:ns-1,:);
					temp2(:,:,1:ns-1,:) = FNs(:,:,2:ns,:);
					temp = temp + (temp2-temp1)/(2*vox(3));
					IntenN = IntenN + dt * temp;

				disp('done!')
			end
			Output = IntenN;
			disp('ADF is Finnished!')
		end
	end
	methods(Static = true)
		function [gv,gp,gs]=CalcGradient(Inten,Vox)
			sd = size(Inten);
			nv = sd(1);
			np = sd(2);
			ns = sd(3);
			temp1 = zeros(sd);
			temp2 = temp1;
			% V-direction
			temp1(2:nv,:,:,:) = Inten(1:nv-1,:,:,:);
			temp2(1:nv-1,:,:,:) = Inten(2:nv,:,:,:);
			gv = (temp2-temp1)/(2*Vox(1));
			% P-direction
			temp1(:,2:np,:,:) = Inten(:,1:np-1,:,:);
			temp2(:,1:np-1,:,:) = Inten(:,2:np,:,:);
			gp = (temp2-temp1)/(2*Vox(2));
			% S-direction
			temp1(:,:,2:ns,:) = Inten(:,:,1:ns-1,:);
			temp2(:,:,1:ns-1,:) = Inten(:,:,2:ns,:);
			gs = (temp2-temp1)/(2*Vox(3));

		end

		function gradcoeff=CalcDiffusionFunc(gv,gp,gs,K)
			% using the c1 and the multichannel alg.:
			sd = size(gv);
			temp = zeros(sd(1:3));
			for i = 1:sd(4)
				temp = temp + gv(:,:,:,i).^2 + gs(:,:,:,i).^2 + gp(:,:,:,i).^2;
			end
			gradcoeff = exp(-(temp)/(K^2));
		end
	end
end
