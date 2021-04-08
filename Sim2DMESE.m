classdef Sim2DMESE
% This class is dependent on the Hargreeve's Bloch code (must add the code to MATLAB path)
	properties
		MyInfo %
		% SeqParams: Must be provided by user
		% LookUpTable: Can be generated by this object
		% T2_Range -> to generate LookUpTable
		% B1Range -> to generate LookUpTable
		% T2Dist: T2 distribution used to generate data
		%	% -> T2Values
		%	% -> Weights
		%
		Data	% First dim is the "NumData" variable in "MyInfo" structure, second dim is length of B1Range, third is the number of T2Dist cells, and 4th is echo train length  
		B1Map	% This is to be used in analysis (handled by this object)
		
	end
	
	methods
		function obj = Sim2DMESE(MyInfo)
			disp('Initializing Sim2DMESE object...');
			if ~isfield(MyInfo, "SeqParams")
				error("SeqParams was not found in the input structure!");
			end
			
			if ~isfield(MyInfo, "LookUpTable")
				disp('Generating LookUpTable using sequence parameters...')
				MyInfo.T2Range = MyInfo.T2Dist(1).T2Values;
				MyInfo.LookUpTable = Sim2DMESE.Generate_LookUpTable(MyInfo);
			end
			
			obj.MyInfo = MyInfo;
			
			obj.Data = zeros(MyInfo.NumData, length(MyInfo.B1Range), length(MyInfo.T2Dist), MyInfo.SeqParams.etl);
			
			disp('Simulating multi-component decay signal...')
			
			for j = 1:length(MyInfo.B1Range)
				for k = 1:length(MyInfo.T2Dist)
					basis_decay = MC_MESE_Nima.Calc_basis_decay(obj.MyInfo, obj.MyInfo.B1Range(j));
					obj.Data(1,j,k,:) = MyInfo.T2Dist(k).Weights * basis_decay;
				end
			end
			
			for i = 1:MyInfo.NumData
				obj.Data(i,:,:,:) = obj.Data(1,:,:,:);
			end
			
			obj.B1Map = zeros(MyInfo.NumData, numel(MyInfo.B1Range), numel(MyInfo.T2Dist));
			
			for j = 1:length(MyInfo.B1Range)
				obj.B1Map(:,j,:) = MyInfo.B1Range(j);
			end
			
			% Creating lookup table for MC_Analyzer
			dummy = MC_MESE_Nima('Dummy');
			dummy.MyInfo = obj.MyInfo;
			dummy.MyInfo.T1 = 1;
			dummy = MC_LookUpTable_Preparation(dummy, 0.5:1e-2:1.5, logspace(log10(8e-3),log10(2),60));
			obj.MyInfo.MC_Analyzer.LookUpTable = dummy.MyInfo.LookUpTable;
			disp('Sim2DMESE object created!')
		end
		
		function [Maps, Maps_B1] = MC_Analyzer(obj, SNR, B1_diff)
			if nargin > 3
				flag_b1diff =  true;
			end
			dummy = MC_MESE_Nima('Dummy');
			dummy.B1Map = obj.B1Map;
			dummy.MyInfo = obj.MyInfo;
			dummy.MyInfo.LookUpTable = obj.MyInfo.MC_Analyzer.LookUpTable;
			dummy.MyInfo.MC_Analyzer = [];
			% Normalize data
			temp = obj.Data(:,:,:,1);
			for i = 1:size(obj.Data,4)
				obj.Data(:,:,:,i) = obj.Data(:,:,:,i)./temp;
			end
			% Add Noise
			dummy.Data = SimClass.ADD_Noise(obj.Data, SNR, 1);
			%
			% Analyze Data
			
			if flag_b1diff
				for i = 1:numel(B1_diff)
					dummy.B1Map = obj.B1Map + B1_diff(i);
					[~,Maps_B1(i)] =  MESE_Calculate_Distribution(dummy, false, 0);
					Maps_B1(i).Residuals = [];
					Maps_B1(i).Distribution = [];
				end
				Maps = [];
			else
			
				[~,Maps_B1] =  MESE_Calculate_Distribution(dummy, false, 0);
			
				% Estimate B1 using MC analysis
				dummy.B1Map = MC_Estimate_B1Map(dummy, false, 0);
			
				[~,Maps] = MESE_Calculate_Distribution(dummy, false, 0);
			end
		end
	end
	
	methods(Static = true)
		function LookUpTable = Generate_LookUpTable(MyInfo)
			if ~isfield(MyInfo, "B1Range")
				MyInfo.B1Range = 0.5:1e-2:1.2;
			end
			
			if ~isfield(MyInfo, "T2Range")
				MyInfo.T2Range = logspace(log10(8e-3), log10(2), 60);
			end
			MyInfo.T1 = 1;
			LookUpTable = MC_MESE_Nima.Create_LookUp_Table(MyInfo);
		end
	end
end