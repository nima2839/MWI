classdef TestClass
   properties
       MyInfo   % Contains general info about the input data
       % MyInfo Fields:
		     % Mask
		     % Method: the method to acquire LFG=> 1=Du method and 2= Hwang method
		     % Vox: Voxel size
		     % FirstTE
		     % EchoSpacing
			 % EchoIndexes: used in the LFG correction code
		     % Algorithm: 'trust-region-reflective' or 'levenberg-marquardt' for NLLS fitting
			 % BipolarFlag: to apply bipolar phase correction
       Mag
       Phase
       SizeData
       LFGC
       Flag_UseLFGC % If true, LFGC data would be used for calculations
       Params_3PM   % Parameters calculated via 3PM will be stored in this variable
       Res_3PM
       Params_C3PM   % Parameters calculated via Complex-3PM will be stored in this variable
       Res_C3PM
       Params_S3PM
       Res_S3PM
       RSC          % Single Component Relaxation Parameter
       Res_SC
       Params_2PM
       Res_2PM
       Flag_UseSC   %If true, uses single component values in 3PM model
       Description
       RunTime
       NNLS
       MWF_3PM
       MWF_S3PM
       MWF_C3PM
       MWF_2PM
       % LFG values
       Gv
       Gp
       Gs

	   Freq_bg	% Background frequency
	   Phi0
   end
   methods
        function obj = TestClass(dataMag,dataPhase,myinfo)
         obj.Mag = dataMag;
         obj.Phase = dataPhase;
         obj.MyInfo = myinfo;
         obj.SizeData = size(dataMag);
         obj.Flag_UseLFGC = false;
         if ~isfield(myinfo,'BipolarFlag')
           obj.MyInfo.BipolarFlag = false;
         end
         disp('TestClass Initiated!')
       end

        function obj = Normal_Procedure(obj)
         obj = CalcLFGC(obj);
         obj = Calc_SC(obj,2);
         %obj = Calc_2PM(obj);
         obj = Calc_S3PM(obj);
         obj = Calc_3PM(obj);
         %obj = Calc_Complex3PM(obj);
       end

	    function obj = SetMag(obj,Mag)
         obj.Mag = Mag;
		 obj.SizeData = size(Mag);
       end

	    function obj = SetLFGC(obj,Mag)
         obj.LFGC = Mag;
		 obj.SizeData = size(Mag);
       end

        function obj = CalcLFGC(obj,method)
          if nargin > 1
              obj.Info.Method =  method;
          end
          [obj.LFGC,obj.Gp, obj.Gv, obj.Gs] = LFG_Correction(obj.Mag, obj.Phase, obj.MyInfo);
          obj.LFGC(isnan(obj.LFGC)) = 0;
          obj.Flag_UseLFGC = true;
       end

        function obj = Calc_3PM(obj,X0)
			np = obj.SizeData(1);
			nv = obj.SizeData(2);
			ns = obj.SizeData(3);
			mask = obj.MyInfo.Mask;
			if obj.Flag_UseLFGC
				mag = obj.LFGC;
			else
				mag = obj.Mag;
			end
			%ei = obj.EchoIndexes;
			% reshaping for parallel processing purposes
			mag = reshape(mag,  [np*nv*ns, obj.SizeData(4)]); 
			params = zeros(np*nv*ns , 8); % the number 8 is from the ThreePoolM_NLLS code 
			res = zeros(1, np*nv*ns);
			mask = reshape(mask, size(res));
			Info = obj.MyInfo;
			
			% initializing with values from https://doi.org/10.1016/j.neuroimage.2015.03.081
			% Frequency of MW has been initilized to 0 instead of 5
			% Upper boundary for MW T2s is set to 15 ms instead of 25 ms
			if nargin < 2
				X0 = [0.1,   10e-3,	0,		0.6,	64e-3,	0.3,	48e-3,		0];
			end
			lb = [0,     3e-3,	-75,	0,		25e-3,	0,		25e-3,		-25];
			ub = [2,	 15e-3,	75,		2,		150e-3,	2,		150e-3,		25];
			
			tic
			disp('3PM Started..!')
			parfor i = 1:numel(res)
				if mask(i) > 0
					[params(i,:), res(i)] = ThreePoolM_NLLS(mag(i,:),Info,X0,lb,ub);
				end
			end
			obj.Params_3PM = reshape(params, [np,nv,ns, size(params,2)]);
			obj.Res_3PM = reshape(res, [np,nv,ns]);
			obj.MWF_3PM = reshape(params(:,1)./((params(:,1) + params(:,4) + params(:,6))), [np,nv,ns]);
			disp('Finished fitting 3PM!')
			toc
		end

		function obj = Calc_Multi_Seed(obj)
			tic
			sd = obj.SizeData;
			Seed{1} = [0.1,   6e-3,		0,	0.6,	64e-3,	0.3,	48e-3,	   0];
			Seed{2} = [0.1,   10e-3,	0,	0.6,	64e-3,	0.3,	48e-3,	   0];
			Seed{3} = [0.1,   14e-3,	0,	0.6,	64e-3,	0.3,	48e-3,	   0];
			Seed{4} = [0.1,   18e-3,	0,	0.6,	64e-3,	0.3,	48e-3,	   0];
			Params = zeros([sd(1)*sd(2)*sd(3),8]);
			Res = inf*ones(sd(1)*sd(2)*sd(3),1);
			for i = 1:4
				disp(strcat("Processing Seed#",string(i)));
				temp = Calc_3PM(obj, Seed{i});
				temp_Params = reshape(temp.Params_3PM, [sd(1)*sd(2)*sd(3), size(temp.Params_3PM,4)]);
				temp_Res = reshape(temp.Res_3PM, [sd(1)*sd(2)*sd(3),1]);
				idx = find(temp_Res < Res);
				Res(idx) = temp_Res(idx);
				Params(idx,:) = temp_Params(idx,:);
			end
			obj.Params_3PM = reshape(Params, [sd(1:3),8]);
			obj.Res_3PM = reshape(Res, [sd(1:3)]);
			obj.MWF_3PM = obj.Params_3PM(:,:,:,1).*((obj.Params_3PM(:,:,:,1) + obj.Params_3PM(:,:,:,4) + obj.Params_3PM(:,:,:,6)).^-1);
			disp('Multi-seed process finished!');
			toc
		end
		
        function obj = Calc_S3PM(obj)
           np = obj.SizeData(1);
           nv = obj.SizeData(2);
           ns = obj.SizeData(3);
           mask = obj.MyInfo.Mask;
           mag = obj.Mag;
           gp = obj.Gp;
           gv = obj.Gv;
           gs = obj.Gs;
           Gradient = sqrt(gp.^2 + gv.^2 + gs.^2);
           Vox = obj.MyInfo.Vox;
           Distance = sqrt(sum(Vox(:).^2));
           %ei = obj.EchoIndexes;
           % reshaping for parallel processing purposes
		   mag = reshape(mag,  [np*nv*ns, obj.SizeData(4)]); 
           params = zeros(np*nv*ns , 9); % the number 9 is from the Sinc3PM_NLLS code 
           res = zeros(1, np*nv*ns);
		   mask = reshape(mask, size(res));
           Info = obj.MyInfo;
		   
           X0 = [0.1,   60,	  0,	0.7,	30,	0.2,	15,	   0,  0];
           lb = [0,     40,	 0,	0,	  10,	0,	  0.1,	0,  0];
           ub = [2,	  300,	25,	2,	  40,	2,	  40,	   25,  20];

           tic
           disp('S3PM Started..!')
           parfor i = 1:numel(res)
                if mask(i) > 0
                    [params(i,:), res(i)] = Sinc3PM_NLLS(mag(i,:),Info,X0,lb,ub);
                end
           end
           obj.Params_S3PM = reshape(params, [np,nv,ns, size(params,2)]);
           obj.Res_S3PM = reshape(res, [np,nv,ns]);
           obj.MWF_S3PM = reshape(params(:,1)./((params(:,1) + params(:,4) + params(:,6))), [np,nv,ns]);
           disp('Finished fitting S3PM!')
           toc
       end


	    function obj = Calc_Freq_bg(obj)
         tic
			   disp('Calculating Freq_bg!')
		     echo1(:,:,:) = obj.Phase(:,:,:,1);
			   echo2(:,:,:) = obj.Phase(:,:,:,3);
			   es = obj.MyInfo.EchoSpacing * (2);
			   p2 = exp(1i*echo2);
		     p1 = exp(-1i*echo1);
			   obj.Freq_bg = angle(p1.*p2)/(es*2*pi);
	       disp('Done!')
	       toc
       end

		function obj = Calc_phi_zero(obj)
         obj = Calc_Freq_bg(obj);
			   tic
         disp('Calculating Phi0!')
         pn(:,:,:) = exp(1i*obj.Phase(:,:,:,1));
         TE = obj.MyInfo.FirstTE;
         temp = exp(-1i*TE*obj.Freq_bg*2*pi);
         obj.Phi0 = angle(pn.*temp);
         disp('Done!')
         toc
       end

        function obj = Bipolar_Phase_Correction(obj)
         disp('Bipolar_Phase_Correction Started...!');
         tic
         TE = obj.MyInfo.EchoSpacing + obj.MyInfo.FirstTE;
         w = 2*pi*obj.Freq_bg;
         p(:,:,:) = exp(-1i*obj.Phase(:,:,:,2)).*exp(1i*(w*TE+obj.Phi0));
         for k = 2:2:obj.SizeData(4)
          % obj.Phase(:,:,:,k) = angle(exp(1i*(TE*w+obj.Phase(:,:,:,k-1))));
            obj.Phase(:,:,:,k) = angle(exp(1i*obj.Phase(:,:,:,k)).*p);
         end

         toc
         disp('Done!');
       end

	    function obj = Calc_Complex3PM(obj)
		    obj = Calc_phi_zero(obj);
            if obj.MyInfo.BipolarFlag
				obj = Bipolar_Phase_Correction(obj);
			end
			np = obj.SizeData(1);
			nv = obj.SizeData(2);
			ns = obj.SizeData(3);
			mask = obj.MyInfo.Mask;
			if obj.Flag_UseLFGC
				mag = obj.LFGC;
			else
				mag = obj.Mag;
			end
            %ei = obj.EchoIndexes;
			% reshaping for parallel processing purposes
			
			params = zeros(np*nv*ns , 9); % the number 9 is from the Complex3PM code 
			res = zeros(1, np*nv*ns);
			mask = reshape(mask, size(res));
			Info = obj.MyInfo;
			
			X0 = [0.1,   60,	  0,	0.7,	30,0,	0.2,	25,	   0];
	        lb = [0,     30,	 -25,	0,	  10,-25,	0,	  0.1,	-25];
	        ub = [2,	  300,	25,	2,	  40,25,	2,	  40,	   25];

		    fbg = obj.Freq_bg;
		    phi0 = obj.Phi0;
		    signal = mag.*exp(1i*obj.Phase);
			signal = reshape(signal,  [np*nv*ns, obj.SizeData(4)]); 
			
			tic;
			disp('Complex 3PM Started..!')
			parfor i = 1:numel(res)
                if mask(i) > 0
					Rx0 = X0;
					tlb = lb;
					tub = ub;
					tempInfo = Info;
					tempInfo.Phi0 = phi0(i,j,k);
					Rx0(3) = fbg(i,j,k);
					tlb(3) = Rx0(3) - 75;
					tub(3) = Rx0(3) + 75;
					Rx0(6) = fbg(i,j,k);
					tlb(6) = Rx0(3) - 25;
					tub(6) = Rx0(3) + 25;
					Rx0(9) = fbg(i,j,k);
					tlb(9) = Rx0(3) - 25;
					tub(9) = Rx0(3) + 25;
					[params(i,:), res(i)] = Complex3PM(signal(i,:),tempInfo,Rx0,tlb,tub);
                end
			end
			obj.Params_C3PM = reshape(params, [np,nv,ns, size(params,2)]);
			obj.Res_C3PM = reshape(res, [np,nv,ns]);
			obj.MWF_C3PM = reshape(params(:,1)./((params(:,1) + params(:,4) + params(:,7))), [np,nv,ns]);
			disp('Finished fitting C3PM!')
			toc
       end

        function obj = Calc_SC(obj,Method)
			if nargin < 2
				Method = 2;
			end
			
			if Method == 1
				disp('Using NLLS Method!')
			elseif Method == 2
				disp('Using Log Method!')
			else
				disp('Invalid Method for SC!')
				return;
			end
			
			
			disp('Excluding first 4 echoes for single component fitting!...')
			e1 = 5; % Based on Gelderen 2012 to exclude first 4 echoes of the decay
			np = obj.SizeData(1);
			nv = obj.SizeData(2);
			ns = obj.SizeData(3);
			
			if obj.Flag_UseLFGC
				mag = obj.LFGC;
			else
				mag = obj.Mag;
			end
			
			% reshaping for parallel processing purposes
			rc = zeros(1,np*nv*ns);
			res = rc;
			mag = reshape(mag,  [np*nv*ns, obj.SizeData(4)]);
			mask = reshape(obj.MyInfo.Mask, size(res));
			Info = obj.MyInfo;
			tic
			disp('Single Component Fitting Started!...');
			parfor i = 1:numel(res)
				if mask(i) > 0
					if Method == 1
						[rc(i) , res(i)] = SingleComponentNLLS(mag(i, e1:end),Info);
					else
						[rc(i) , res(i)] = SingleComponentFitting(mag(i, e1:end),Info);
					end
                end
			end
			obj.RSC = reshape(rc, size(obj.MyInfo.Mask));
			obj.Res_SC = reshape(res, size(obj.RSC));
			disp('SCF Completed!');
			obj.Flag_UseSC = true;
			toc
       end

        function obj = Calc_2PM(obj)
			np = obj.SizeData(1);
			nv = obj.SizeData(2);
			ns = obj.SizeData(3);
			mask = obj.MyInfo.Mask;
			if obj.Flag_UseLFGC
				mag = obj.LFGC;
			else
				mag = obj.Mag;
			end
			
			% reshaping for parallel processing purposes
			mag = reshape(mag,  [np*nv*ns, obj.SizeData(4)]); 
			params = zeros(np*nv*ns , 5); % the number 5 is from the TwoPoolModel_NLLS code 
			res = zeros(1, np*nv*ns);
			mask = reshape(mask, size(res));
			Info = obj.MyInfo;
			
			X0 = [0.1,	100,	5,		0.9, 	20];
			lb = [0,		40,		-25,	0.5, 	0.2];
			ub = [2,		300,	25,		2, 		40];
			
			tic
			disp('2PM Started..!')
			parfor i = 1:numel(res)
				if mask(i) > 0
					[params(i,:), res(i)] = TwoPoolModel_NLLS(mag(i,:),Info,X0,lb,ub);
				end
			end
			obj.Params_2PM = reshape(params, [np,nv,ns, size(params,2)]);
			obj.Res_2PM = reshape(res, [np,nv,ns]);
			obj.MWF_2PM = reshape(params(:,1)./((params(:,1) + params(:,4))), [np,nv,ns]);
			disp('Finished fitting 2PM!')
			toc
       end

        function obj = Calc_NNLS(obj)
		   if ~obj.Flag_UseLFGC
			obj = CalcLFGC(obj);
		   end
		   Chi2Factor = 1.02;
		   
		   % Calculating Basis Decay Curves
		   time =  obj.MyInfo.FirstTE:obj.MyInfo.EchoSpacing:(obj.MyInfo.FirstTE + (obj.SizeData(4) - 1)*obj.MyInfo.EchoSpacing);
		   T2Range = [1e-3, 2];
		   nT2 = 60;
		   T2_times=logspace(log10(T2Range(1)),log10(T2Range(2)),nT2);
		   basis_decay = zeros(obj.SizeData(4),nT2);
		   for x=1:nT2
			echo_amp = exp(-time/T2_times(x)); % Nima : T1 vector is used
			basis_decay(:,x) = echo_amp';
		   end
		   %initializing
		   % reshaping for parallel processing
		   SD = size(obj.LFGC);
		   
		   Mask = reshape(obj.MyInfo.Mask, [SD(1)*SD(2)*SD(3),1]);
		   img = reshape(obj.LFGC, [SD(1)*SD(2)*SD(3), SD(4)]);
		   distributions =zeros(size(img,1),nT2);
			gdnmap = zeros(size(Mask));
			ggmmap = zeros(size(Mask));
			gvamap = zeros(size(Mask));
			FNRmap = zeros(size(Mask));
			ResMap = zeros(size(img));
			obs_weights = ones(SD(4),1);
		   
		   tic;
		   disp('NNLS fitting started...');
		   % Iterating through voxels
			parfor i = 1:numel(Mask)
				if Mask(i) > 0
					decay_data = reshape(img(i,:), [SD(4),1]);
					if decay_data(1) > 0
						[T2_dis,~,~] = Nima_UBC_NNLS(basis_decay, decay_data, obs_weights, Chi2Factor);
						distributions(i,:) = T2_dis;
						% Compute parameters of distribution
						gdnmap(i) = sum(T2_dis);
						ggmmap(i) = exp(dot(T2_dis,log(T2_times))/sum(T2_dis));
						gvamap(i) = exp(sum((log(T2_times)-log(ggmmap(i))).^2.*T2_dis')./sum(T2_dis)) - 1;
						decay_calc = basis_decay*T2_dis;
						residuals = decay_calc-decay_data;
						ResMap(i,:) = residuals; % Nima
						FNRmap(i) = sum(T2_dis)/std(residuals); 
					end
				end
			end
			disp('NNLS Completed!');
			toc
			maps.gdn = reshape(gdnmap, SD(1:3));
			maps.ggm = reshape(ggmmap, SD(1:3));
			maps.gva = reshape(gvamap, SD(1:3));
			maps.FNR = reshape(FNRmap, SD(1:3)); 
			maps.Residuals = reshape(ResMap, SD); % Nima
			obj.NNLS.Distribution = reshape(distributions, [SD(1:3), nT2]);
			obj.NNLS.Maps = maps;
			disp('done!')
        end

        function data = GetAllData(obj)
          %if obj.Flag_UseLFGC
          %  data.LFGC = obj.LFGC;
          %else
          %  data.Mag = obj.Mag;
          %end
          data.Info = obj.MyInfo;
          data.RSC = obj.RSC;
          data.resRC = obj.Res_SC;
          data.Params_3PM = obj.Params_3PM;
          data.res3pm = obj.Res_3PM;
          data.Params_S3PM = obj.Params_S3PM;
          data.resS3pm = obj.Res_S3PM;
          data.Params_C3PM = obj.Params_C3PM;
          data.resC3pm = obj.Res_C3PM;
          data.Description = obj.Description;
          data.runtime = obj.RunTime;
          data.Params_2PM = obj.Params_2PM;
          data.res2pm = obj.Res_2PM;
          data.NNLS = obj.NNLS;
          data.MWF_3PM = obj.MWF_3PM;
          data.MWF_S3PM = obj.MWF_S3PM;
          data.MWF_C3PM = obj.MWF_C3PM;
          data.MWF_2PM = obj.MWF_2PM;
       end
   end
end
