classdef SimClass
% GRASE code must be included for T2 simulation
	properties
		MyInfo % Contains necessary info
			% 	Note: all range values must be a vector of size 2, even if the corresponding is a constant value!
			%	Fields:
			%	NumData: Number of data point to be simulated
			%	Times: A vector containing all the echo times
			%	NumWaterComp: Number of water compartment in the simuation
			%	TimeConstRange:  a cell that contains range of each compartment time constant range (seconds)
			%	T1Val: an array containing all T1 values corresponding to each water compartment (seconds)
			%	FractionRange: a cell containing the fraction of each water compartment
			%	FreqRange: a cell ... (for T2 simulations should not exist!)
			% FlipAngle: refocusing flip angle value (degrees)
			%	TrueFAFlag: If false code will estimate FA on its own
			% SNR

		Compartment_T_Map % TimeConstant values of each water compartment
		Compartment_Freq_Map
		Compartment_Fraction_Map
		SimulatedData
	end
	methods
		function obj = SimClass(MyInfo)

			if ~isfield(MyInfo, 'SNR')
				MyInfo.SNR = 0;
			end

			if ~isfield(MyInfo, 'TrueFAFlag')
				MyInfo.TrueFAFlag = false;
			end

			if ~isfield(MyInfo, 'FlipAngle')
				MyInfo.FlipAngle = 180;
			end

			obj.MyInfo = MyInfo;
			if length(MyInfo.TimeConstRange) ~= MyInfo.NumWaterComp || length(MyInfo.FractionRange) ~= MyInfo.NumWaterComp
				error('Number of water compartments and corresponding range values are inconsistent!');
			end

			SimulatedData = zeros(MyInfo.NumData,length(MyInfo.Times));
			Compartment_Fraction_Map = cell(MyInfo.NumWaterComp,1);
			Compartment_T_Map = cell(MyInfo.NumWaterComp,1);
			if isfield(MyInfo,'FreqRange')
				Compartment_Freq_Map = cell(MyInfo.NumWaterComp,1);
			end


			for k = 1:MyInfo.NumWaterComp
				tempT = SimClass.GenerateRandomValues(MyInfo.TimeConstRange{k}, MyInfo.NumData);
				tempFraction = SimClass.GenerateRandomValues(MyInfo.FractionRange{k}, MyInfo.NumData);

				if isfield(MyInfo,'FreqRange')
					Compartment_Freq_Map{k} = SimClass.GenerateRandomValues(MyInfo.FreqRange{k}, MyInfo.NumData);
					temp = cell2mat(arrayfun(@(T,f,a) a * SimClass.CreateDecayCurve(T, MyInfo.Times, f,MyInfo.SNR), tempT, Compartment_Freq_Map{k},tempFraction,'UniformOutput',false));
				else
					%ES = MyInfo.Times(2) - MyInfo.Times(1);
					temp = cell2mat(arrayfun(@(T,a) a * SimClass.GenerateT2DecayCurves(MyInfo.Times,T,MyInfo.T1Val(k),MyInfo.FlipAngle,MyInfo.SNR),...
							tempT,tempFraction,'UniformOutput',false));
				end
				SimulatedData = SimulatedData + reshape(temp, size(SimulatedData'))';
				Compartment_Fraction_Map{k} = tempFraction;
				Compartment_T_Map{k} = tempT;
			end
			obj.Compartment_Fraction_Map = Compartment_Fraction_Map;
			obj.Compartment_T_Map = Compartment_T_Map;
			obj.SimulatedData = SimulatedData;
			if isfield(MyInfo,'FreqRange')
				obj.Compartment_Freq_Map = Compartment_Freq_Map;
			end
		end

		function [Dist, Maps] = NNLS_Fitting(obj)
			MyInfo.FirstTE = obj.MyInfo.Times(1);
			MyInfo.EchoSpacing = obj.MyInfo.Times(2) - obj.MyInfo.Times(1);
			MyInfo.Range_T = [1 120] * 1e-3;
			MyInfo.Num_T =  100;

			temp(1,1,:,:) = obj.SimulatedData(:,:);

			[Dist(:,:),Maps.Times,Maps.res] = NNLS_T_Batch_Fitting(temp, MyInfo);
			threshold = obj.MyInfo.TimeConstRange{obj.MyInfo.NumWaterComp}(2);
			Maps.MWF = SimClass.Find_MWF(Dist, find(Maps.Times < threshold, 1, 'last' ), 'NNLS');
		end

		function [Dist, Maps] = UBC_Fitting(obj)
			temp(1,1,:,:) = abs(obj.SimulatedData(:,:)) *1e3; % Scale is for thresholding in the code

			if obj.MyInfo.TrueFAFlag
				[Maps, Dist(:,:), ~] = T2map_SEcorr(temp, 'SetFlipAngle', obj.MyInfo.FlipAngle);
			else
				[Maps, Dist(:,:), ~] = T2map_SEcorr(temp);%, 'SetFlipAngle', obj.MyInfo.FlipAngle);
			end
			Maps.MWF = SimClass.Find_MWF(Dist, 40, 'NNLS');
		end

		function [Dist, Maps] = UBC_Nima_Fitting(obj, T1)
			if nargin < 2
				T1 = ones(1,200);
			end
			temp(1,1,:,:) = abs(obj.SimulatedData(:,:)) *1e3; % Scale is for thresholding in the code

			if obj.MyInfo.TrueFAFlag
				[Maps, Dist(:,:), ~] = T2map_Nima(temp, 'FlipAngleMap', obj.MyInfo.FlipAngle*ones(size(temp)), 'T1', T1);
			else
				[Maps, Dist(:,:), ~] = T2map_Nima(temp, 'T1', T1);%, 'SetFlipAngle', obj.MyInfo.FlipAngle);
			end
			Maps.MWF = SimClass.Find_MWF(Dist, 40, 'NNLS');
		end

		function Maps = NLLS_Fitting(obj, SNR)
			if nargin < 2
				temp(:,:) = (obj.SimulatedData(:,:));
			else
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
			MyInfo.FirstTE = obj.MyInfo.Times(1);
			MyInfo.EchoSpacing = obj.MyInfo.Times(2) - obj.MyInfo.Times(1);
			MyInfo.MultiStart = true;
			for k = 1:obj.MyInfo.NumData
				[Maps.Params3PM(k,:), Maps.Res3PM(k)] = ThreePoolM_NLLS(abs(temp(k,:)),MyInfo);
				[Maps.ParamsMag3PM(k,:), Maps.ResMag3PM(k)] = Mag3PM(abs(temp(k,:)),MyInfo);
				[Maps.ParamsC3PM(k,:), Maps.ResC3PM(k)] = Complex3PM(temp(k,:),MyInfo);
				%[Maps.Params2PM(k,:), Maps.Res2PM(k)] = TwoPoolModel_NLLS(temp(k,:),MyInfo);
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

		function output = GenerateT2DecayCurves(Times,T2Val,T1Val,FA, SNR)
			if nargin < 5
				SNR = 0;
			end
			T2 = T2Val*(0.95:0.01:1.05);
			weight = zeros(1,length(T2));
			weight(T2 == T2Val) = 1;
			weight = conv(weight, gausswin(length(T2),0.3*length(T2)), 'same');
			ETL = length(Times);
			TE = Times(2) - Times(1);
			output = 0;
			for k = 1:length(T2)
				output = output + weight(k) * EPGdecaycurve(ETL,FA,TE,T2(k),T1Val,180);
			end
			% Normalize output
			output = output / sum(weight(:));
			if SNR > 0
				output = awgn(output,SNR);
			end
		end

		%function output = GenerateT2StarDecayCurves(T2SVal,Freq)

		%end
	end
end
