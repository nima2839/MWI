function Nominal_Angle = Estimate_NominalAngle(B1, Maps)
    % Estimated nominal angle using the estimated flip angle map from the optimization step (Prasloski's method)
    % It is assumed that the optimization step is an unbiased estimator and the the nominal angle is less than 180 degrees
    %
    % 'B1' is the coregistered B1 map (Normalized)
    % 'Maps' output of the estimated method; 
    %
    % Removing voxels containing 'ggm' values over 100ms to avoid vessels and CSF affected areas
    % also areas that decay relatively fast (fat)
    idx = find (B1 > 0.95  & B1 < 1.0  & Maps.ggm < 0.1 & Maps.ggm > 0.05);
    if isempty(idx)
        error('Could not find enough voxels for the process!');
    end
    Selected_Voxels = Maps.alpha(idx)./ B1(idx);
    Nominal_Angle = mean(Selected_Voxels(:));
end
