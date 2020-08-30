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
       NNLS_W
       res_NNLS
       NNLS_Times
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
       end

	     function obj = SetLFGC(obj,Mag)
         obj.LFGC = Mag;
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
           params = zeros([obj.SizeData(1:3),8]);
           res = zeros(obj.SizeData(1:3));
           Info = obj.MyInfo;
		   if nargin < 2
				X0 = [0.1,   10e-3,	  0,	0.6,	64e-3,	0.3,	48e-3,	   0];
		   end
	       lb = [0,     1e-3,	 -25,	0,	  20e-3,	0,	  20e-3,	-25];
	       ub = [2,	  20e-3,	25,	2,	  1,	2,	  1,	   25];
           flag = obj.Flag_UseSC;
           if flag
               RC = obj.RSC;
           else
		       RC = X0(5)*ones(size(res));
		   end
           tic
           disp('3PM Started..!')
           parfor i = 1:np
               tempP = zeros(nv,ns,8);
               tempR = zeros(nv,ns);
                for j = 1:nv
                    for k = 16%1:ns
                        if mask(i,j,k) > 0
                            tmp = squeeze(mag(i,j,k,:));
                         	tmpd = tmp;	%tmp(ei);
                            Rx0 = X0;
                            if flag
                                Rx0(5) = RC(i,j,k) - 3;
                            end
                         	[tempP(j,k,:), tempR(j,k)] = ThreePoolM_NLLS(tmpd,Info,Rx0,lb,ub);
                        end
                    end
                end
                params(i,:,:,:) = tempP;
                res(i,:,:) = tempR;
           end
           obj.Params_3PM = params;
           obj.Res_3PM = res;
           obj.MWF_3PM = params(:,:,:,1).*((params(:,:,:,1) + params(:,:,:,4) + params(:,:,:,6)).^-1);
           disp('Finished fitting 3PM!')
           toc
       end

		function obj = Calc_Multi_Seed(obj)
			tic
			sd = obj.SizeData;
			Seed{1} = [0.1,   6e-3,	  0,	0.6,	64e-3,	0.3,	48e-3,	   0];
			Seed{2} = [0.1,   10e-3,	  0,	0.6,	64e-3,	0.3,	48e-3,	   0];
			Seed{3} = [0.1,   14e-3,	  0,	0.6,	64e-3,	0.3,	48e-3,	   0];
			Seed{4} = [0.1,   18e-3,	  0,	0.6,	64e-3,	0.3,	48e-3,	   0];
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
           params = zeros([obj.SizeData(1:3),9]);
           res = zeros(obj.SizeData(1:3));
           Info = obj.MyInfo;
           X0 = [0.1,   60,	  0,	0.7,	30,	0.2,	15,	   0,  0];
          lb = [0,     40,	 0,	0,	  10,	0,	  0.1,	0,  0];
          ub = [2,	  300,	25,	2,	  40,	2,	  40,	   25,  20];
           flag = obj.Flag_UseSC;
           if flag
               RC = obj.RSC;
           end
           tic
           disp('S3PM Started..!')
           parfor i = 1:np
               tempP = zeros(nv,ns,9);
               tempR = zeros(nv,ns);
                for j = 1:nv
                    for k = 1:ns
                        if mask(i,j,k) > 0
                            tmp = squeeze(mag(i,j,k,:));
                           tmpd = tmp;	%tmp(ei);
                            Rx0 = X0;
                            if flag
                                Rx0(5) = RC(i,j,k) - 3;
                            end
                            Rx0(9) = Gradient(i,j,k) * Distance * 42.575e6;
                           [tempP(j,k,:), tempR(j,k)] = Sinc3PM_NLLS(tmpd,Info,Rx0,lb,ub);
                        end
                    end
                end
                params(i,:,:,:) = tempP;
                res(i,:,:) = tempR;
           end
           obj.Params_S3PM = params;
           obj.Res_S3PM = res;
           obj.MWF_S3PM = params(:,:,:,1).*((params(:,:,:,1) + params(:,:,:,4) + params(:,:,:,6)).^-1);
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
           params = zeros([obj.SizeData(1:3),9]);
           res = zeros(obj.SizeData(1:3));
           Info = obj.MyInfo;
           X0 =   [0.1,   60,	  0,	0.7,	30,0,	0.2,	25,	   0];
	         lb = [0,     30,	 -25,	0,	  10,-25,	0,	  0.1,	-25];
	         ub = [2,	  300,	25,	2,	  40,25,	2,	  40,	   25];
           flag = obj.Flag_UseSC;
           if flag
               RC = obj.RSC;
           end
		       fbg = obj.Freq_bg;
		       phi0 = obj.Phi0;
		       signal = mag.*exp(1i*obj.Phase);
           tic
           disp('Complex 3PM Started..!')
           parfor i = 1:np
               tempP = zeros(nv,ns,9);
               tempR = zeros(nv,ns);
                for j = 1:nv
                    for k = 1:ns
                        if mask(i,j,k) > 0
                          tmp = squeeze(signal(i,j,k,:));
                         	tmpd = tmp;	%tmp(ei);
                          Rx0 = X0;
							            tlb = lb;
							            tub = ub;
                          if flag
                              Rx0(5) = RC(i,j,k) - 3;
                          end
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
                         	[tempP(j,k,:), tempR(j,k)] = Complex3PM(tmpd,tempInfo,Rx0,tlb,tub);
                        end
                    end
                end
                params(i,:,:,:) = tempP;
                res(i,:,:) = tempR;
           end
           obj.Params_C3PM = params;
           obj.Res_C3PM = res;
           obj.MWF_C3PM = params(:,:,:,1).*((params(:,:,:,1) + params(:,:,:,4) + params(:,:,:,7)).^-1);
           disp('Finished fitting C3PM!')
           toc
       end

       function obj = Calc_SC(obj,Method)
          if nargin < 2
              Method = 1;
          end
          if Method == 1
              disp('Using NLLS Method!')
          elseif Method == 2
              disp('Using Log Method!')
          else
              disp('Invalid Method for SC!')
              return;
          end
          rc = zeros(obj.SizeData(1:3));
          res = rc;
          e1 = 1; % testing!
          e2 = 40;
          if e2 > obj.SizeData(4)
              e2 = obj.SizeData(4);
          end
          np = obj.SizeData(1);
          nv = obj.SizeData(2);
          ns = obj.SizeData(3);
          if obj.Flag_UseLFGC
              mag = obj.LFGC;
          else
              mag = obj.Mag;
          end
          mask = obj.MyInfo.Mask;
          Info = obj.MyInfo;
          tic
          disp('Single Component Fitting Started!...');
          parfor i = 1:np
              temprc = zeros(nv,ns);
              tempres = temprc;
              for j = 1:nv
                for k = 1:ns
                    if mask(i,j,k) > 0
                        tmp = squeeze(mag(i,j,k,:));
                        tmpd = tmp(e1:e2);
                        if Method == 1
                            [temprc(j,k) , tempres(j,k)] = SingleComponentNLLS(tmpd,Info);
                        else
                            [temprc(j,k) , tempres(j,k)] = SingleComponentFitting(tmpd,Info);
                        end
                    end
                end
              end
              rc(i,:,:) = temprc;
              res(i,:,:) = tempres;
          end
          obj.RSC = rc;
          obj.Res_SC = res;
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
           %ei = obj.EchoIndexes;
           params = zeros([obj.SizeData(1:3),5]);
           res = zeros(obj.SizeData(1:3));
           Info = obj.MyInfo;
           X0 = [0.1, 60,	5,	 0.9, 30];
           lb = [0,   40,   0,	0.5, 1];
           ub = [2,	300,	25,	2, 40];
           flag = obj.Flag_UseSC;
           if flag
               RC = obj.RSC;
           end
           tic
           disp('2PM Started..!')
          parfor i = 1:np
               tempP = zeros(nv,ns,5);
               tempR = zeros(nv,ns);
                for j = 1:nv
                    for k = 1:ns
                        if mask(i,j,k) > 0
                            tmp = squeeze(mag(i,j,k,:));
                            tmpd = tmp;	%(ei);
                            Rx0 = X0;
                            if flag
                                Rx0(5) = RC(i,j,k) - 3;
                            end
                       	    [tempP(j,k,:), tempR(j,k)] = TwoPoolModel_NLLS(tmpd,Info,Rx0,lb,ub);
                        end
                    end
                end
                params(i,:,:,:) = tempP;
                res(i,:,:) = tempR;
           end
           obj.Params_2PM = params;
           obj.Res_2PM = res;
           obj.MWF_2PM = obj.Params_2PM(:,:,:,1).*(obj.Params_2PM(:,:,:,1)+obj.Params_2PM(:,:,:,4)).^-1;
           disp('Finished fitting 2PM!')
           toc
       end

       function obj = Calc_NNLS(obj)
           obj = CalcLFGC(obj);
           LFG.Gs = obj.Gs;
           LFG.Gv = obj.Gv;
           LFG.Gp = obj.Gp;
           et.FirstEchoTime = obj.MyInfo.FirstTE;
           et.EchoSpacing = obj.MyInfo.EchoSpacing;
           options.T_range = [1 120];
           options.Number_T = 120;
           options.Vox = obj.MyInfo.Vox;
           tic
           disp('Calculating Weights Using NNLS...')
           [W, T, r] = NNLS_T_Batch_Fitting(obj.Mag,et,options,LFG);
           obj.NNLS_W = W;
           obj.NNLS_Times =T;
           obj.res_NNLS = r;
           disp('done!')
           toc
         end

       function data = GetAllData(obj)
          if obj.Flag_UseLFGC
            data.LFGC = obj.LFGC;
          else
            data.Mag = obj.Mag;
          end
          data.Info = obj.MyInfo;
          data.RC = obj.RSC;
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
          data.NNLS_W = obj.NNLS_W;
          data.res_NNLS = obj.res_NNLS;
          data.NNLS_Times = obj.NNLS_Times;
          data.MWF_3PM = obj.MWF_3PM;
          data.MWF_S3PM = obj.MWF_S3PM;
          data.MWF_C3PM = obj.MWF_C3PM;
          data.MWF_2PM = obj.MWF_2PM;
       end
   end
end
