classdef SimClass
	properties
		MyInfo % Contains necessary info, this is user defined property
			% 	Note: all range values must be a vector of size 2, even if the corresponding is a constant value!
			%	Fields:
			%	NumData: Number of data point to be simulated
			%	Times: A vector containing all the echo times
			%	NumWaterComp: Number of water compartment in the simuation
			%	TimeConstRange:  a cell that contains range of each compartment time constant range (seconds)
			%		- edit1: if "T2Dist" is defined then this should not be defined at all!
			%	T1Val: an array containing all T1 values corresponding to each water compartment (seconds)
			%		- edit1: if "T2Dist" is defined then it should be a vector with the same length!
			%	FractionRange: a cell containing the fraction of each water compartment (only for T2* simuation)!
			%	FreqRange: a cell ... (Must not exist for T2 simulations!)
			%	FlipAngle: refocusing flip angle value (degrees)
			%	TrueFAFlag: If false code will estimate FA on its own
			%	SNR: uses adittive white Gaussian noise!
			%	T2Dist: This is a struct containing two vectors->
			%		- T2Values: contains T2 values in the distribution
			%		- Weights: contain the corresponding weights
			%
			
		Dist % Contains distribution weights corresponding to "T2Dist"  option!
		Compartment_T_Map % TimeConstant values of each water compartment
		Compartment_Freq_Map
		Compartment_Fraction_Map
		SimulatedData
	end
	methods
		function obj = SimClass(MyInfo)
			MyInfo.NumData = MyInfo.NumData + mod(MyInfo.NumData, maxNumCompThreads); % for acceleration purposes
			if ~isfield(MyInfo, 'TrueFAFlag')
				MyInfo.TrueFAFlag = false;
			end

			if ~isfield(MyInfo, 'FlipAngle')
				MyInfo.FlipAngle = 180;
			end
			
			if ~isfield(MyInfo, 'T1Val')
				error('MyInfo.T1Val needs to be defined!')
			end

			obj.MyInfo = MyInfo;
			if ~isfield(MyInfo, 'T2Dist')
				if length(MyInfo.TimeConstRange) ~= MyInfo.NumWaterComp || length(MyInfo.FractionRange) ~= MyInfo.NumWaterComp
					error('Number of water compartments and corresponding range values are inconsistent!');
				end
			end
			if ~isfield(MyInfo, 'SNR') 
				MyInfo.SNR = 100;
			elseif MyInfo.SNR == 0
				MyInfo.SNR = 100;
			end
			obj.MyInfo = MyInfo;
			SimulatedData = zeros(MyInfo.NumData,length(MyInfo.Times));
			
			if isfield(MyInfo,'FreqRange')
				Compartment_Freq_Map = cell(MyInfo.NumWaterComp,1);
			end
			tic;
			if ~isfield(MyInfo, 'T2Dist')
			Compartment_Fraction_Map = cell(MyInfo.NumWaterComp,1);
			Compartment_T_Map = zeros(MyInfo.NumWaterComp,MyInfo.NumData);
			
				for k = 1:MyInfo.NumWaterComp
					tempT = SimClass.GenerateRandomValues(MyInfo.TimeConstRange{k}, MyInfo.NumData);
					tempFraction = SimClass.GenerateRandomValues(MyInfo.FractionRange{k}, MyInfo.NumData);
					if isfield(MyInfo,'FreqRange')
						Compartment_Freq_Map{k} = SimClass.GenerateRandomValues(MyInfo.FreqRange{k}, MyInfo.NumData);
						temp = cell2mat(arrayfun(@(T,f,a) a * SimClass.CreateDecayCurve(T, MyInfo.Times, f,MyInfo.SNR), tempT, Compartment_Freq_Map{k},tempFraction,'UniformOutput',false));
					else

						temp = cell2mat(arrayfun(@(T,a) a * SimClass.GenerateT2DecayCurves_Delta(MyInfo.Times,T,MyInfo.T1Val(k),MyInfo.FlipAngle,MyInfo.SNR),...
								tempT,tempFraction,'UniformOutput',false));
					end
					SimulatedData = SimulatedData + reshape(temp, size(SimulatedData'))';
					Compartment_Fraction_Map{k} = tempFraction;
					Compartment_T_Map(k,:) = tempT;
				end
			else
				Compartment_Fraction_Map = {};
				Compartment_T_Map = {};
				
				SimResult = SimClass.GenerateT2DecayCurves(MyInfo.Times, MyInfo.T2Dist, MyInfo.T1Val, MyInfo.FlipAngle);
				for i = 1:MyInfo.NumData
					SimulatedData(i,:) = awgn(SimResult(:), MyInfo.SNR);
				end
			end
			SimTime = toc;
			obj.Compartment_Fraction_Map = Compartment_Fraction_Map;
			obj.Compartment_T_Map = Compartment_T_Map;
			obj.SimulatedData = SimulatedData;
			if isfield(MyInfo,'FreqRange')
				obj.Compartment_Freq_Map = Compartment_Freq_Map;
			end
			%disp(['Total Simulation Time: ', string(SimTime)])
		end

		function [Dist, Maps] = NNLS_Fitting(obj) % do not use this funtion, this is from the old code!
			TempInfo.FirstTE = obj.MyInfo.Times(1);
			TempInfo.EchoSpacing = obj.MyInfo.Times(2) - obj.MyInfo.Times(1);
			TempInfo.Range_T = [1 120] * 1e-3;
			TempInfo.Num_T =  100;

			temp(1,1,:,:) = obj.SimulatedData(:,:);

			[Dist(:,:),Maps.Times,Maps.res] = NNLS_T_Batch_Fitting(temp, TempInfo);
			threshold = obj.MyInfo.TimeConstRange{obj.MyInfo.NumWaterComp}(2);
			Maps.MWF = SimClass.Find_MWF(Dist, find(Maps.Times < threshold, 1, 'last' ), 'NNLS');
		end

		function [Dist, Maps] = UBC_Fitting(obj)
			temp(1,1,:,:) = abs(obj.SimulatedData(:,:));

			if obj.MyInfo.TrueFAFlag
				[Maps, Dist(:,:), ~] = T2map_SEcorr(temp, 'SetFlipAngle', obj.MyInfo.FlipAngle,'Threshold', 0,'nT2', 60,'T2Range', [0.008, 2], 'MinRefAngle', 100);
			else
				[Maps, Dist(:,:), ~] = T2map_SEcorr(temp,'Threshold', 0,'nT2', 60,'T2Range', [0.008, 2], 'MinRefAngle', 100);%, 'SetFlipAngle', obj.MyInfo.FlipAngle);
			end
			Maps.MWF = SimClass.Find_MWF(Dist, 18, 'NNLS');
		end

		function [Dist, Maps] = UBC_Nima_Fitting(obj, T1)
			if nargin < 2
				T1 = ones(1,200);
			end
			nT2 = 60;
			ns = obj.MyInfo.NumData / maxNumCompThreads;
			ne = length(obj.MyInfo.Times);
			temp = reshape(abs(obj.SimulatedData(:,:)), maxNumCompThreads,1,ns,ne); 

			if obj.MyInfo.TrueFAFlag
			disp('Yellow');
				[Maps, Dist, ~] = T2map_Nima(temp, 'FlipAngleMap', obj.MyInfo.FlipAngle*ones(maxNumCompThreads,1,ns), 'T1', T1,'Threshold', 0,'nT2', nT2,'T2Range', [0.008, 2], 'MinRefAngle', 100);
			else
				[Maps, Dist, ~] = T2map_Nima(temp, 'T1', T1,'Threshold', 0,'nT2', nT2,'T2Range', [0.008, 2], 'MinRefAngle', 100);%, 'SetFlipAngle', obj.MyInfo.FlipAngle);
			end
			Maps.MWF = SimClass.Find_MWF(Dist, 18, 'NNLS');
		end

		function Maps = NLLS_Fitting(obj, SNR)
			if nargin < 2
				temp(:,:) = (obj.SimulatedData(:,:));
			else
				temp = zeros(obj.MyInfo.NumData, length(obj.MyInfo.Times));
				for k = 1:obj.MyInfo.NumData
					temp(k,:) = awgn((obj.SimulatedData(k,:)), SNR);
				end
			end
			tic
			Maps.Params3PM = zeros(obj.MyInfo.NumData,8);
			Maps.ParamsC3PM = zeros(obj.MyInfo.NumData,9);
			Maps.ParamsMag3PM = zeros(obj.MyInfo.NumData,6);
			%Maps.Params2PM = zeros(obj.MyInfo.NumData,5);
			Maps.Res3PM = zeros(obj.MyInfo.NumData,1);
			Maps.ResC3PM = zeros(obj.MyInfo.NumData,1);
			Maps.ResMag3PM = zeros(obj.MyInfo.NumData,1);
			%Maps.Res2PM = zeros(obj.MyInfo.NumData,1);
			TempInfo.FirstTE = obj.MyInfo.Times(1);
			TempInfo.EchoSpacing = obj.MyInfo.Times(2) - obj.MyInfo.Times(1);
			TempInfo.MultiStart = true;
			for k = 1:obj.MyInfo.NumData
				[Maps.Params3PM(k,:), Maps.Res3PM(k)] = ThreePoolM_NLLS(abs(temp(k,:)), TempInfo);
				[Maps.ParamsMag3PM(k,:), Maps.ResMag3PM(k)] = Mag3PM(abs(temp(k,:)), TempInfo);
				[Maps.ParamsC3PM(k,:), Maps.ResC3PM(k)] = Complex3PM(temp(k,:), TempInfo);
				%[Maps.Params2PM(k,:), Maps.Res2PM(k)] = TwoPoolModel_NLLS(temp(k,:), TempInfo);
			end
			Maps.MWF3PM = Maps.Params3PM(:,1).* (Maps.Params3PM(:,1)+Maps.Params3PM(:,4)+Maps.Params3PM(:,6)).^-1;
			Maps.MWFC3PM = Maps.ParamsC3PM(:,1).* (Maps.ParamsC3PM(:,1)+Maps.ParamsC3PM(:,4)+Maps.ParamsC3PM(:,7)).^-1;
			Maps.MWFMag3PM = Maps.ParamsMag3PM(:,1).* (Maps.ParamsMag3PM(:,1)+Maps.ParamsMag3PM(:,3)+Maps.ParamsMag3PM(:,5)).^-1;
			%Maps.MWF2PM = Maps.Params2PM(:,1).* (Maps.Params2PM(:,1)+Maps.Params2PM(:,4)).^-1;
			toc
		end

	end
	methods(Static = true)
		function output = CreateDecayCurve(TimeConstant,t,Freq,SNR)
			if nargin < 3
				Freq = 0;
			elseif nargin < 4
				SNR = 0;
			end
			output = exp(-t*((1/TimeConstant) - 1i*2*pi*Freq));
			if SNR > 0
				output = awgn(output, SNR);
			end
		end

		function output = GenerateRandomValues(range,n)
			output = (range(2) - range(1)) * rand(1,n) + range(1);
		end

		function output = Find_MWF(Dist, index, Mode)
			if strcmp(Mode, 'NNLS')
				output = squeeze(squeeze(sum(Dist(:,1:index),2)./sum(Dist,2)));
			end
		end
		
		function output = GenerateT2DecayCurves_Delta(Times, T2 , T1, FA, SNR)
			% This function uses the "GenerateT2DecayCurves" for alternative input option
			T2Dist.T2Values = T2;
			T2Dist.Weights = 1;
			output = SimClass.GenerateT2DecayCurves(Times,T2Dist,T1,FA, SNR);
		end

		function output = GenerateT2DecayCurves(Times,T2Dist,T1Val,FA, SNR)
			if nargin < 5
				SNR = 0;
			end
			
			nT2 = length(T2Dist.T2Values);
			ETL = length(Times);
			TE = Times(2) - Times(1);
			RefCon = 180; % Predefined
	
			basis_decay = SimClass.Calc_basis_decay(ETL, nT2, FA, TE, T2Dist.T2Values, T1Val, RefCon);
			
			sbd = size(basis_decay);
			output = basis_decay * reshape(T2Dist.Weights, [sbd(2), 1]);
			% Normalize output
			output = output / sum(1);
			if SNR > 0
				output = awgn(output,SNR);
			end
		end

		function basis_decay = Calc_basis_decay(nechs, nT2, alpha, TE, T2_times, T1, RefCon)
			basis_decay=zeros(nechs,nT2);
			% Compute the NNLS basis over T2 space
			if alpha == 0 % this prevents the code from freezing
				alpha = 180;
			end
			for x=1:nT2
				echo_amp = EPGdecaycurve(nechs, alpha, TE, T2_times(x), T1(x), RefCon); % Nima : T1 vector is used
				basis_decay(:,x) = echo_amp';
			end
		end
	end
end
