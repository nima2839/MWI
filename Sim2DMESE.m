classdef Sim2DMESE
% This class is dependent on the modified version of Kelly's code (must add the code to MATLAB path)
	properties
		MyInfo %
		% SeqParams: Must be provided by user
		% LookUpTable: Can be generated by this object
		% T2_Range -> to generate LookUpTable
		% B1_Range -> to generate LookUpTable
		% T2Dist: T2 distribution used to generate data
		%	% -> T2Values
		%	% -> Weights
		%
		Data	% First dim is the "NumData" variable in "MyInfo" structure, second dim is length of B1_Range, third is the number of T2Dist cells, and 4th is echo train length  
		B1Map	% This is to be used in analysis (handled by this object)
		
	end
	
	methods
		function obj = Sim2DMESE(MyInfo)
			disp('Initializing Sim2DMESE object...');
			if ~isfield(MyInfo, "SeqParams")
				error("SeqParams was not found in the input structure!");
			end
			
			if ~isfield(MyInfo, "LookUpTable")
				disp('Genrating LookUpTable using sequence paramters..')
				MyInfo.T2_Range = MyInfo.T2Dist(1).T2Values;
				obj.MyInfo.LookUpTable = Generate_LookUpTable(MyInfo);
			end
			
			obj.MyInfo = MyInfo;
			
			obj.Data = zeros(MyInfo.NumData, length(MyInfo.B1_Range), length(T2Dist), MyInfo.SeqParams.etl);
			
			disp('Simulating multi-component decay signal...')
			
			for j = 1:length(MyInfo.B1_Range)
				for k = 1:length(T2Dist)
					basis_decay = MC_MESE_SLR.Calc_basis_decay(obj.MyInfo, obj.MyInfo.B1_Range(j));
					obj.Data(1,j,k,:) = T2Dist(k).Weights * basis_decay;
				end
			end
			
			for i = 2:MyInfo.NumData
				obj.Data(i,:,:,:) = obj.Data(1,:,:,:);
			end
			
			obj.B1Map = zeros(MyInfo.NumData, length(MyInfo.B1_Range), length(T2Dist));
			
			for j = 1:length(MyInfo.B1_Range)
				obj.B1Map(:,j,:) = MyInfo.B1_Range(j);
			end
			disp('Sim2DMESE object created!')
		end
		
		function [Maps, Maps_B1] = MC_Analyzer(obj, SNR)
			dummy = MC_MESE_SLR('Dummy');
			dummy.B1Map = obj.B1Map;
			dummy.MyInfo = obj.MyInfo;
			% Add Noise
			dummy.Data = SimClass.ADD_Noise(obj.Data, SNR, obj.Data(1));
			%
			% Analyze Data
			
			[~,Maps_B1] = Calculate_Distribution(dummy);
			
			% Estimate B1 using MC analysis
			dummy.B1Map = MC_Estimate_B1Map(dummy);
			
			[~,Maps] = Calculate_Distribution(dummy);
		end
	end
	
	methods(Static = true)
		function LookUpTable = Generate_LookUpTable(MyInfo)
			if ~isfield(MyInfo, "B1_Range")
				MyInfo.B1_Range = 0.5:1e-2:1.2;
			end
			
			if ~isfield(MyInfo, "T2_Range")
				MyInfo.T2_Range = logspace(log10(8e-3), log10(2), 60);
			end
			LookUpTable = MC_MESE_SLR.Create_LookUp_Table(MyInfo);
		end
	end
end