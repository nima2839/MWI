function [maps,distributions, T1] = T2map_Nima(image,varargin)
    %
    % Nima Changes:
    %   - Some code encapsulation in the form of function definition was set for better readablity;
    %   - Changed the shape of the for loop for better parallel processing optimization and debugging;
    %	- Regularization is always enabled in this version!
    %	- A flip angle map can be given as input; 'FlipAngleMap' must be the registered B1 map in the correct scale and dimension;
    %	- Default values are changed to match the parameters used in my thesis;
    %	- Observation Weights has been added;
    %	- Added T1 map vector as an input for their corresponding T2 in T2 distribution;
    %	- T1 value is returned in case default values in case alternate values were to be used;
    %	- Residuals of fit are now saved in the maps structure; To save space this could be set to empty array on the user's end;
    %	- nCores option has been removed; MATLAB now takes care of the parallel processing optimizations;
    %	- Alternate fitting function is called to work with observation weights ("Nima_UBC_NNLS");
    %	- Threshold's default value set back to 200; Must be set to zero for simulated data;
    %	- Ratio of second echo to the first echo of the fitted curve is now recorded in the maps output;
    %   - SNR estimation was added as the ratio of signal power to residual (estimated noise) power; 
    %
    
    %
    % [maps,distributions] = T2map_SEcorr(image,...) % Nima: this was the original function
    %
    % Description:
    %   Uses NNLS to compute T2 distributions in the presence of stimulated
    %   echos by optimizing the refocusing pulse flip angle.  Records parameter
    %   maps and T2 distributions for further partitioning.
    %
    % Inputs:
    %   image: 4-D array with intensity data as (row,column,slice,echo)
    %   ...: A series of optional Property/Value pairs to modify settings.
    %     Defaults are given in brackets:
    %       'TE': Interecho spacing (0.01)
    %       'nT2': Number of T2 times to use (40)
    %       'T2Range': Min and Max T2 values ([0.015,2.000])
    %       'T1': Assumed value of T1 (1)
    %       'Threshold': First echo intensity cutoff for empty voxels (200)
    %       'Reg': Regularization routine to use, options are:
    %              'no': do not regularize the solution
    %              'chi2': use Chi2Factor based regularization (default)
    %              'lcurve': use L-Curve based regularization
    %       'Chi2Factor': Constraint on chi^2 used for regularization (Reg must
    %                     be set to 'chi2'!) (1.02)
    %       'RefCon': Refocusing Pulse Control Angle (180)
    %       'MinRefAngle': Minimum refocusing angle for EPG optimization (50)
    %       'nAngles': Number of angles used in EPG optimization (8)
    %       'SetFlipAngle': Instead of optimizing flip angle, uses this flip
    %                       angle for all voxels (not set)
    %       'nCores': Number of processor cores to use (6)
    %       'Save_regparam': yes/no option to include the regularization
    %                        paramter mu and the resulting chi^2 factor as
    %                        two outputs within the maps structure (mu=NaN and
    %                        chi2factor=1 if Reg=no) ('no')
    %       'Save_NNLS_basis': yes/no option to include a 5-D matrix of NNLS
    %                          basis matrices as another output within the maps
    %                          structure ('no')
    %
    % Ouputs:
    %   maps: Structure containing 3D maps of the following parameters
    %       -gdn, general density
    %       -ggm, general geometric mean
    %       -gva, general variance
    %       -FNR, fit to noise ratio (gdn/stdev(residuals))
    %       -alpha, refocusing pulse flip angle
    %   distributions: 4-D matrix containing T2 distributions.
    %
    % External Calls:
    %   EPGdecaycurve.m
    %   lsqnonneg_reg.m
    %   lsqnonneg_lcurve.m
    %
    
    %==========================================================================
    % Parse inputs and apply default values when necessary
    %==========================================================================
    % Make image double
    image=double(image);
    %
    % Create input parser object
    p=inputParser;
    % Define all input values
    p.addRequired('image',@(x)isa(x,'double') && ndims(x)==4);
    p.addParamValue('TE',0.01,@(x)isnumeric(x) && isscalar(x) && x>=0.001 && x<=1); % TE modification _J
    %p.addParamValue('T1',1,@(x)isnumeric(x) && isscalar(x) && x>=10 && x<=0.001); % Nima : a vector is expected now
    p.addParamValue('RefCon',       180,@(x)isnumeric(x) && isscalar(x) && x<=180 && x>=1);
    p.addParamValue('Threshold',    200,@(x)isnumeric(x) && isscalar(x)); 
    p.addParamValue('Chi2Factor',   1.02,@(x)isnumeric(x)  && isscalar(x) && x>1);
    p.addParamValue('nT2',          60,@(x)isnumeric(x) && isscalar(x) && x>=10 && x<=300);
    p.addParamValue('T2Range',      [0.008,2],@(x)isnumeric(x) && length(x)==2 && x(2)>x(1) && x(1)>=0.001 && x(2)<=10);
    p.addParamValue('MinRefAngle',  100,@(x)isnumeric(x) && isscalar(x) && x>1 && x<180);
    p.addParamValue('nAngles',      8,@(x)isnumeric(x) && isscalar(x) && x>1);
    p.addParamValue('Reg',          'lcurve',@(x)any(strcmp(x,{'no','chi2','lcurve'})));
    %p.addParamValue('SetFlipAngle',0,@(x)(isnumeric(x) && isscalar(x)));
    %p.addParamValue('nCores',4,@(x)isnumeric(x) && isscalar(x) && x>=1 && x<=8);
    p.addParamValue('Save_regparam','no',@(x)any(strcmp(x,{'yes','no'})));
    
    % Nima: Setting Observation Weights
    p.addParamValue('Observation_Weights', [], @(x)isnumeric(x) && ndims(x) == 2)
    
    % Nima: This variable determines where in the distribution we would find the myelin water T2 cut-off
    p.addParamValue('MWF_CutOff',   40e-3,@(x)isnumeric(x) && isscalar(x)); % milliseconds 
    % Nima:Set FlipAngleMap
    p.addParamValue('FlipAngleMap', [],@(x)isa(x,'double') && ndims(x)==3);
    p.addParamValue('T1',           [], @(x)isa(x,'double'));
    %
    p.addParamValue('Save_NNLS_basis','no',@(x)any(strcmp(x,{'yes','no'})));
    % Parse inputs (MATLAB will throw an error here if any variables fail validation)
    p.parse(image,varargin{:});
    % Define all variables from the inputParser Results
    TE = p.Results.TE;
    T1 = p.Results.T1; % Nima : This is a vector now
    
    RefCon = p.Results.RefCon;
    Threshold = p.Results.Threshold;
    MWF_CutOff = p.Results.MWF_CutOff;
    Chi2Factor = p.Results.Chi2Factor;
    nT2 = p.Results.nT2;
    T2Range = p.Results.T2Range;
    minangle = p.Results.MinRefAngle;
    nangles = p.Results.nAngles;
    % reg = p.Results.Reg;
    % Nima: Change this value according to the new changes
    alphamap = p.Results.FlipAngleMap;
    faset=	~isempty(alphamap);	%p.Results.SetFlipAngle;
    
    if isempty(T1)
        T1 = ones(1, nT2);
        %disp('T1 is set to 1 seconds for all!');
    end
    
    
    % nCores=p.Results.nCores;
    % nCores= 4; % in case of error message about Cores;
    
    savereg = strcmp(p.Results.Save_regparam,'yes');
    saveNNLS = strcmp(p.Results.Save_NNLS_basis,'yes');
    
    %==========================================================================
    % Initialize all the data
    %==========================================================================
    % Start the clock
    % tstart=tic;
    % Find size of the data
    [nrows,ncols,nslices,nechs] = size(image);
    
    % Nima: Setting Observation Weights
    obs_weights = p.Results.Observation_Weights;
    if isempty(obs_weights)
        obs_weights = ones(1,nechs);
    end
    
    
    % Reshaping input for parallel processing
    Map_Numel = nrows*ncols*nslices;
    image = reshape(image, [Map_Numel, nechs]);
    
    % Initialize map matrices
    gdnmap = zeros(Map_Numel, 1);
    ggmmap = zeros(Map_Numel, 1);
    gvamap = zeros(Map_Numel, 1);
    SNRmap = zeros(Map_Numel, 1); 
    FNRmap = zeros(Map_Numel, 1); 
    E2_E1map = zeros(Map_Numel, 1); % Nima: Ratio of second echo to first in the fitter curve
    
    if faset == 0
        alphamap = zeros(Map_Numel, 1);
    else
        alphamap = reshape(alphamap, [Map_Numel, 1]);
    end
    
    distributions = nan*ones(Map_Numel, nT2);
    ResMap = nan*ones(Map_Numel, nechs); % Nima
    mumap = nan*ones(Map_Numel, 1);
    chi2map = nan*ones(Map_Numel, 1);
    
    %
    % decay_basis=nan*ones(nrows,ncols,nslices,nechs,nT2);
    %
    if saveNNLS
        decay_basis = nan*ones(Map_Numel,nT2);
    end
    
    %==========================================================================
    % Find the basis matrices for each flip angle
    %==========================================================================
    % Initialize parameters and variable for angle optimization
    T2_times = logspace(log10(T2Range(1)), log10(T2Range(2)), nT2);
    if faset==0
        flip_angles=linspace(minangle, 180, nangles);
        % basis_angles is a 1xnangles cell array that will contain the decay bases of each angle
        nangles = length(flip_angles);
        basis_angles=cell(nangles);
        % Loop to compute each basis and assign them to a cell in the array
        basis_decay=zeros(nechs,nT2);
        for a=1:nangles
            for x=1:nT2
                echo_amp = EPGdecaycurve(nechs, flip_angles(a), TE, T2_times(x), T1(x), RefCon); % Nima : different T1 value could be used for each T2
                basis_decay(:,x) = echo_amp';
            end
            basis_angles{a}=basis_decay;
        end
        basis_decay_faset=[];  %ignore
    else
       basis_angles=[];  %ignore
       flip_angles=[];  %ignore
       basis_decay_faset=zeros(nechs,nT2);
       for x=1:nT2
           echo_amp = EPGdecaycurve(nechs, faset, TE, T2_times(x), T1(x), RefCon); % Nima : different T1 value could be used for each T2
           basis_decay_faset(:,x) = echo_amp';
       end
    end
    %==========================================================================
    % Process all pixels
    %==========================================================================
    % Main triple for-loop to run through each pixel in the image
    
    %if matlabpool('size')==0
    %    matlabpool('open',num2str(nCores))
    %    poolopenflag=1;
    %else
    %    poolopenflag=0;
    %end
    
    
    parfor i = 1:Map_Numel
        % Extract decay curve from the voxel
        decay_data = reshape(squeeze(image(i,:)), [nechs,1]);
        % Conditional loop to reject low signal pixels
        try
            if decay_data(1) > Threshold
                
                %obs_weights = ones(size(decay_data)); % Nima: set observation weights here
                %obs_weights(1:10) = 1.5;
                if faset == 0
                    %======================================================
                    % Find optimum flip angle
                    %======================================================
                    alphamap(i) = Estimate_Alpha(basis_angles, nangles, decay_data, reshape(obs_weights, size(decay_data)), flip_angles);
                %else
                %    alpha(col,slice) = alphamap(row,col,slice);
                end
                %======================================================
                % Fit basis matrix using alpha
                %======================================================
                basis_decay = Calc_basis_decay(nechs, nT2, alphamap(i), TE, T2_times, T1, RefCon);
    
                %==========================================================
                % Calculate T2 distribution and global parameters
                %==========================================================
                % Find distribution depending on regularization routine
                [T2_dis,mu,chi2] = Nima_UBC_NNLS(basis_decay, decay_data, reshape(obs_weights, size(decay_data)), Chi2Factor);
    
                % Compute parameters of distribution
                decay_calc = basis_decay*T2_dis;
                residuals = decay_calc-decay_data;
    
                gdnmap(i) = sum(T2_dis);
                ggmmap(i) = exp(dot(T2_dis,log(T2_times))/sum(T2_dis));
                gvamap(i) = exp(sum((log(T2_times)-log(ggmmap(i))).^2.*T2_dis')./sum(T2_dis)) - 1;
                SNRmap(i) = sum(decay_calc.^2)/sum(residuals.^2); % Nima
                FNRmap(i) = sum(T2_dis)/std(residuals);
                E2_E1map(i) = decay_calc(2)/decay_calc(1); % Nima
                distributions(i,:) = T2_dis;
                ResMap(i,:) = residuals; % Nima
                mumap(i) = mu;
                chi2map(i) = chi2;
                if saveNNLS
                    decay_basis(i,:) = basis_decay;
                end
            end
        catch ME
            disp(strcat('An error occured at index = ',string(i),', more details below:'));
            disp(ME);
        end
    end
    
    
    %if poolopenflag==1
    %    matlabpool close
    %end
    
    % Assign outputs in the correct shape
    maps.gdn = reshape(gdnmap, [nrows,ncols,nslices]);
    maps.ggm = reshape(ggmmap, [nrows,ncols,nslices]);
    maps.gva = reshape(gvamap, [nrows,ncols,nslices]);
    maps.alpha = reshape(alphamap, [nrows,ncols,nslices]);
    maps.FNR = reshape(FNRmap, [nrows,ncols,nslices]); 
    maps.SNR = reshape(SNRmap, [nrows,ncols,nslices]); 
    maps.Residuals = reshape(ResMap, [nrows,ncols,nslices,nechs]); % Nima
    maps.E2_E1 = reshape(E2_E1map, [nrows,ncols,nslices]); % Nima
    [~, idx] = min(abs(T2_times - MWF_CutOff));
    maps.MWF = reshape(squeeze(sum(distributions(:,1:idx),2)./sum(distributions,2)), [nrows,ncols,nslices]); % Nima
    distributions = reshape(distributions, [nrows,ncols,nslices,nT2]);
    if savereg
        maps.mu = reshape(mumap, [nrows,ncols,nslices]);
        maps.chi2factor = reshape(chi2map, [nrows,ncols,nslices]);
    end
    %
    % maps.NNLS_basis=decay_basis;
    %
    if saveNNLS
        maps.NNLS_basis=decay_basis;
    end
    
    
    end
    
    function alpha = Estimate_Alpha(basis_angles, nangles, decay_data, obs_weights, flip_angles)
      % Fit each basis and  find chi-squared
      
      
      chi2_alpha = zeros(1,nangles);
      
      for a=1:nangles
          [T2_dis_ls, ~, ~] = Nima_UBC_NNLS(basis_angles{a},decay_data, obs_weights); % Nima: This is to be tested!
          decay_pred=basis_angles{a}*T2_dis_ls;
          chi2_alpha(a)=sum((decay_data - decay_pred).^2);
      end
      % Find the minimum chi-squared and the corresponding angle
      alpha_spline = flip_angles(1):0.001:flip_angles(end);
      chi2_spline=interp1(flip_angles,chi2_alpha,alpha_spline,'spline');
      [~,index] = min(chi2_spline);
      alpha = alpha_spline(index);
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
    