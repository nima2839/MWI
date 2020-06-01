function Summerize_Sim_Results(filename)
	% A costumized summary of SimClass results

	Models = load(filename, 'Models');
	Sim.Models = cell(size(Models));
	Sim.FA = Models{1}.FA;
	Sim.SNR = Models{1}.SNR;
	nSNR = length(Models{1}.SNR);
	nFA = length(Models{1}.FA);
	for M = 1:length(Models)
		
		Sim.Models{M}.mMWF = zeros(nSNR , nFA); 
		Sim.Models{M}.sdMWF = mMWF;
		Sim.Models{M}.T_mMWF = mMWF;
		Sim.Models{M}.T_sdMWF = mMWF;
		
		Sim.Models{M}.mAlpha = zeros(nSNR , nFA); 
		Sim.Models{M}.sdAlpha = mAlpha;
		Sim.Models{M}.T_mAlpha = mAlpha;
		Sim.Models{M}.T_sdAlpha = mAlpha;

		for i = 1:nSNR		
			for j = 1:nFA
				Sim.Models{M}.mMWF(i,j) = mean(Models{M}.Maps{i,j}.MWF(:));
				Sim.Models{M}.sdMWF(i,j) = std((Models{M}.Maps{i,j}.MWF(:) ));
				Sim.Models{M}.T_mMWF(i,j) = mean(Models{M}.TrueFA_Maps{i,j}.MWF(:));
				Sim.Models{M}.T_sdMWF(i,j) = std((Models{M}.TrueFA_Maps{i,j}.MWF(:)));
				
				Sim.Models{M}.mAlpha (i,j) = mean(Models{M}.Maps{i,j}.alpha(:));
				Sim.Models{M}.sdAlpha (i,j) = std((Models{M}.Maps{i,j}.alpha(:) ));
				Sim.Models{M}.T_mAlpha (i,j) = mean(Models{M}.TrueFA_Maps{i,j}.alpha(:));
				Sim.Models{M}.T_sdAlpha (i,j) = std((Models{M}.TrueFA_Maps{i,j}.alpha(:)));
			end
		end
	end
	save(strcat('Summerized_',filename),'Sim','-v7.3');
	disp('Summerized!!!');
end