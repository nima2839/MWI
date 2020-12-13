function FA_Map = Create_FA_Map(Options)
% Options:
%		CenterFA: Flip angle at the center of the map;
%		Size : Matrix size;
%		Vox: Voxel size;
%		deltaFA: determins the rate that FA changes moving away from the center;
%

	if ~isstruct(Options)
		Size = Options;
		clear Options;
		Options.Size = Size;
	end
	
	if ~isfield(Options, 'CenterFA')
		Options.CenterFA = 200;
	end
	
	if ~isfield(Options, 'Vox')
		Options.Vox = floor(Options.Size / min(Options.Size)).^-1;
	end
	
	if ~isfield(Options, 'deltaFA')
		temp = (Options.Size ./ Options.Vox) .^ 2;
		Options.deltaFA = 5 * Options.CenterFA / sum(temp);
	end
	
	%disp(Options)
	FA_Map = zeros(Options.Size);
	Center = floor(Options.Size / 2);
	for i = 1:Options.Size(1)
		for j = 1:Options.Size(2)
			for k = 1:Options.Size(3)
				FA_Map(i,j,k) = Options.CenterFA - Calc_Distance_to_Center([i j k]) * Options.deltaFA;
			end
		end
	end
	function Distance = Calc_Distance_to_Center(Point)
		temp = Center - Point;
		temp = temp .* Options.Vox;
		temp = temp .^ 2;
		Distance = sqrt(sum(temp(:)));
	end
end