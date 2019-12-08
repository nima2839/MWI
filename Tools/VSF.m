classdef VSF
% Voxel Spread Function Method for Correction of Magnetic Field Inhomogeneity by Yablonsky et.al.
% using the same algorithm as Du et. al(2009) to deal with Phase unwrapping.
% Must check how much RAM is needed to peroform the calculations!!!!!
	properties(Constant)
		Gamma = 2*pi*42.575e6;
	end
	properties
		%F_Func		% F-function values
		Mag				% The assumtion is that first echo equals S(0)
		Phase
		MyInfo		% Contains: FirstTE, EchoSpacing, Vox, EchoIndexes(for finding gradient,etc.)
		bn				% Linear approximation of the distribution of inhomogeneous macroscopic magnetic field
		Flag_bn		% This flag shows weather bn values has been calculated or not.
		phi_zero 	%
		% RF field inhomogeneities
		phi_np
		phi_nv
		phi_ns
		% Field gradients
		gp
		gv
		gs
		% Phase dispersion across the 1D voxel
		qp
		qv
		qs
		% Dimesionless parameter 'q'obj.
		q_val_p
		q_val_v
		q_val_s
	end

	methods
		function obj = VSF(Mag,Phase,MyInfo)
			if nargin ~= 3
				error('Inadequate number of input for VSF(Mag,Phase,Inf)');
			end
			if ~isfield(MyInfo, 'Vox')
				error('MyInfo.Vox is missing!');
			end
			if ~isfield(MyInfo, 'FirstTE')
				error('MyInfo.FirstTE is missing!');
			end
			if ~isfield(MyInfo, 'EchoSpacing')
				error('MyInfo.EchoSpacing is missing!');
			end
			%sd = size(Mag);
			%obj.F_Func = zeros(sd);
			obj.MyInfo = MyInfo;
			% Hanning window helps with the Gibbs ringing artifact and effects the definition of 'm' by user
			complex_data = Mag.*exp(1i*Phase);
			sd = size(Mag);
			parfor i = 1:sd(4)
				complex_data(:,:,:,i) = obj.ApplyHanningWindow(complex_data(:,:,:,i), obj.Hann3D(5,5,5));
			end
			obj.Mag = abs(complex_data);
			obj.Phase = angle(complex_data);
			obj.Flag_bn = false;
			if ~isfield(MyInfo, 'EchoIndexes')
				obj.MyInfo.EchoIndexes = [1, 3];
			end

			obj.q_val_p = (-0.5+(1/sd(1))) : (1/sd(1)) : 0.5;
			obj.q_val_v = (-0.5+(1/sd(2))) : (1/sd(2)) : 0.5;
			obj.q_val_s = (-0.5+(1/sd(3))) : (1/sd(3)) : 0.5;
		end
		function obj = Calc_bn(obj)
			if obj.Flag_bn
				return;
			else
				obj.Flag_bn = true;
			end
			echo1(:,:,:) = obj.Phase(:,:,:,obj.MyInfo.EchoIndexes(1));
			echo2(:,:,:) = obj.Phase(:,:,:,obj.MyInfo.EchoIndexes(2));
			es = obj.MyInfo.EchoSpacing * (obj.MyInfo.EchoIndexes(2) - obj.MyInfo.EchoIndexes(1));
			p2 = exp(1i*echo2);
			p1 = exp(-1i*echo1);
			obj.bn = angle(p1.*p2)/(es*obj.Gamma);
		end
		function obj = Calc_phi_zero(obj)
			obj = Calc_bn(obj);
			pn(:,:,:) = exp(1i*obj.Phase(:,:,:,obj.MyInfo.EchoIndexes(1)));
			TE = obj.MyInfo.FirstTE + obj.MyInfo.EchoSpacing * (obj.MyInfo.EchoIndexes(1) - 1);
			temp = exp(-1i*TE*obj.Gamma*obj.bn);
			obj.phi_zero = angle(pn.*temp);
		end
		function obj = Calc_phi_nj(obj)
			disp('Calc_phi_nj Started!...');
			tic
			obj = Calc_phi_zero(obj);
			sd = size(obj.Phase);
			np = sd(1);
			nv = sd(2);
			ns = sd(3);
			% P-direction
			F = zeros(sd(1:3));
			S = F;
			F(1:np-1,:,:) = obj.phi_zero(2:np,:,:);
			S(2:np,:,:) = obj.phi_zero(1:np-1,:,:);
			obj.phi_np = angle(exp(1i*F).*exp(-1i*S))/(2*obj.MyInfo.Vox(1));
			% V-direction
			F = zeros(sd(1:3));
			S = F;
			F(:,1:nv-1,:) = obj.phi_zero(:,2:nv,:);
			S(:,2:nv,:) = obj.phi_zero(:,1:nv-1,:);
			obj.phi_nv = angle(exp(1i*F).*exp(-1i*S))/(2*obj.MyInfo.Vox(2));
			% S-direction
			F = zeros(sd(1:3));
			S = F;
			F(:,:,1:ns-1) = obj.phi_zero(:,:,2:ns);
			S(:,:,2:ns) = obj.phi_zero(:,:,1:ns-1);
			obj.phi_ns = angle(exp(1i*F).*exp(-1i*S))/(2*obj.MyInfo.Vox(3));
			toc
			disp('Calc_phi_nj Done!');
		end
		function obj = Calc_gradients(obj)
			disp('Calc_gradients Started!...');
			tic
			obj = Calc_bn(obj);
			sd = size(obj.Phase);
			np = sd(1);
			nv = sd(2);
			ns = sd(3);
			% P-direction
			F = zeros(sd(1:3));
			S = F;
			F(1:np-1,:,:) = obj.bn(2:np,:,:);
			S(2:np,:,:) = obj.bn(1:np-1,:,:);
			obj.gp = angle(exp(1i*F).*exp(-1i*S))/(2*obj.MyInfo.Vox(1));
			% V-direction
			F = zeros(sd(1:3));
			S = F;
			F(:,1:nv-1,:) = obj.bn(:,2:nv,:);
			S(:,2:nv,:) = obj.bn(:,1:nv-1,:);
			obj.gv = angle(exp(1i*F).*exp(-1i*S))/(2*obj.MyInfo.Vox(2));
			% S-direction
			F = zeros(sd(1:3));
			S = F;
			F(:,:,1:ns-1) = obj.bn(:,:,2:ns);
			S(:,:,2:ns) = obj.bn(:,:,1:ns-1);
			obj.gs = angle(exp(1i*F).*exp(-1i*S))/(2*obj.MyInfo.Vox(3));
			toc
			disp('Calc_gradients Done!');
		end
		function obj = Calc_qj(obj)
			disp('Calc_qj Started!...');
			tic
			obj = Calc_phi_nj(obj);
			obj = Calc_gradients(obj);
			sd = size(obj.Mag);
			temp_qp = zeros(sd);
			temp_qv = zeros(sd);
			temp_qs = zeros(sd);
			TE = obj.MyInfo.FirstTE:obj.MyInfo.EchoSpacing:((sd(4)-1)*obj.MyInfo.EchoSpacing + obj.MyInfo.FirstTE);
			gamma = obj.Gamma;
			Vp = obj.MyInfo.Vox(1);
			Vv = obj.MyInfo.Vox(2);
			Vs = obj.MyInfo.Vox(3);
			temp_phi_np = obj.phi_np;
			temp_phi_nv = obj.phi_nv;
			temp_phi_ns = obj.phi_ns;
			temp_gp = obj.gp;
			temp_gv = obj.gv;
			temp_gs = obj.gs;
			parfor e = 1:sd(4)
					temp_qp(:,:,:,e) =  (gamma * TE(e) * temp_gp + temp_phi_np) * Vp/(2*pi);
					temp_qv(:,:,:,e) =  (gamma * TE(e) * temp_gv + temp_phi_nv) * Vv/(2*pi);
					temp_qs(:,:,:,e) =  (gamma * TE(e) * temp_gs + temp_phi_ns) * Vs/(2*pi);
			end
			obj.qp = temp_qp;
			obj.qv = temp_qv;
			obj.qs = temp_qs;
			toc
			disp('Calc_qj Done!');
		end
		function F_Func = Calc_F_Func(obj)
			disp('Calc_F_Func Started!...');
			tic
			obj = Calc_qj(obj);
			sd = size(obj.Mag);
			numP = sd(1);
			numV = sd(2);
			numS = sd(3);
			numE = sd(4);
			F_Func = zeros(sd);
			s0(:,:,:) = obj.Mag(:,:,:,1);
			for e = 1:numE % testing without parfor
				for k = 3:numS-2
					for j = 3:numV-2
						for i = 3:numP-2
							if s0(i,j,k) == 0
								F_Func(i,j,k,e) = inf;
							else
								F_Func(i,j,k,e) = Calc_Inner_Sigma(obj,i,j,k,e) / s0(i,j,k);
							end
						end
					end
				end
				fprintf('%3.2f Percent...\n',100*double(e/numE));
			end
			toc
			disp('Calc_F_Func Done!');
		end
		function out = Calc_Inner_Sigma(obj,i,j,k,e)
			% Calculation of this method is dependent on the deffinition of 'm' (defined in the article) by the user
			out = 0;
			%sd = size(obj.Mag);
			%if (i < 3) || (i > sd(1)-2) || (j < 3) || (j > sd(2)-2) || (k < 3) || (k > sd(3)-2)
			%	out = inf;
			%	return;
			%end
			%n = [i, j, k];
			%s0(:,:,:) = obj.Mag(:,:,:,1);
			mp = i;%-2:i+2
			mv = j;%-2:j+2
					for ms = k-2:k+2
						%m = [mp, mv, ms];
						out = out + obj.Mag(mp,mv,ms,1) * exp(1i*obj.Phase(mp,mv,ms,e)) * obj.Calc_mu(k,ms,obj.q_val_s,obj.qs(mp,mv,ms));% *...
													%obj.Calc_mu(j,mv,obj.q_val_v,obj.qv(mp,mv,ms)) * obj.Calc_mu(i,mp,obj.q_val_p,obj.qp(mp,mv,ms));;%Calc_mu_product(obj,n,m);
					end
			%	end
			%end
		end
		function out = Calc_mu_product(obj,n,m)
			np = n(1);nv = n(2);ns = n(3);
			mp = m(1);mv = m(2);ms = m(3);
			%mu_p = obj.Calc_mu(np,mp,obj.q_val_p,obj.qp(mp,mv,ms));
		%	mu_v = obj.Calc_mu(nv,mv,obj.q_val_v,obj.qv(mp,mv,ms));
			mu_s = obj.Calc_mu(ns,ms,obj.q_val_s,obj.qs(mp,mv,ms));
			out = mu_s;% * mu_p * mu_v;
		end
	end
	methods(Static)
		function out = Calc_mu(n,m,q,qm)
				out = (sinc(q-qm).* (cos( pi* q).^2))* exp(1i*2*pi*q'*(n-m));
		end
		function out = ApplyHanningWindow(complex_data,HannWindow)
			out = convn(complex_data, HannWindow, 'same');
		end
		function w = Hann3D(nx,ny,nz)
			% creats a 3D Hanning window
			w = bsxfun(@times,bsxfun(@times,hann(nx),hann(ny).'),permute(hann(nz),[3 2 1]));
		end
	end
end
