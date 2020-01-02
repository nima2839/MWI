function Normalized = Normalize_FA_Map(FA_Map)
% Normalize input matrix to values under 180 degree, as EPG algorithm does
	Normalized = 180 - FA_Map;
	Normalized = 180 - abs(Normalized);
end