function [maps,distributions, T1] = T2map_Nima(image,varargin)
%
% Nima Changes:
%	- Regularization is disabled in this version!
%	- A flip angle map can be given as input;
%	- Flip angles for spline function now go over 180, the results are mirrored to the other side of 180;
%	- Number of flip angles used for spline func. is kept unchanged, but it can be increased fo better accuracy!
%	- Min flip angle is increased from 50 to 60;
%	- Observation Weights has been added
% - Add T1 map vector as an input for their corresponding T2 in T2 distribution

%
% [maps,distributions] = T2map_SEcorr(image,...)
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
% Created by Thomas Prasloski
% email: tprasloski@gmail.com
% Ver. 3.2, August 2012

%==========================================================================
% Parse inputs and apply default values when necessary
%==========================================================================
% Make image double
image=double(image);
%
%image_corr = zeros(size(image));
% image_corr2 = zeros(size(image));
%
% Create input parser object
p=inputParser;
% Define all input values
p.addRequired('image',@(x)isa(x,'double') && ndims(x)==4);
p.addParamValue('TE',0.01,@(x)isnumeric(x) && isscalar(x) && x>=0.001 && x<=1); % TE modification _J
%p.addParamValue('T1',1,@(x)isnumeric(x) && isscalar(x) && x>=10 && x<=0.001); % Nima : a vector is expected now
p.addParamValue('RefCon',180,@(x)isnumeric(x) && isscalar(x) && x<=180 && x>=1);
p.addParamValue('Threshold',0,@(x)isnumeric(x) && isscalar(x)); % 200 to 0; Since we are using synthetic data
p.addParamValue('Chi2Factor',1.02,@(x)isnumeric(x)  && isscalar(x) && x>1);
p.addParamValue('nT2',200,@(x)isnumeric(x) && isscalar(x) && x>=10 && x<=300);
p.addParamValue('T2Range',[0.015,2],@(x)isnumeric(x) && length(x)==2 && x(2)>x(1) && x(1)>=0.001 && x(2)<=10);
p.addParamValue('MinRefAngle',60,@(x)isnumeric(x) && isscalar(x) && x>1 && x<180);
p.addParamValue('nAngles',8,@(x)isnumeric(x) && isscalar(x) && x>1);
p.addParamValue('Reg','lcurve',@(x)any(strcmp(x,{'no','chi2','lcurve'})));
p.addParamValue('SetFlipAngle',0,@(x)(isnumeric(x) && isscalar(x)));
p.addParamValue('nCores',4,@(x)isnumeric(x) && isscalar(x) && x>=1 && x<=8);
p.addParamValue('Save_regparam','no',@(x)any(strcmp(x,{'yes','no'})));

% Nima:Set FlipAngleMap
p.addParamValue('FlipAngleMap',[],@(x)isa(x,'double') && ndims(x)==4);
p.addParamValue('T1',ones(1,200),@(x)isa(x,'double'));
%
p.addParamValue('Save_NNLS_basis','no',@(x)any(strcmp(x,{'yes','no'})));
% Parse inputs (MATLAB will throw an error here if any variables fail validation)
p.parse(image,varargin{:});
% Define all variables from the inputParser Results
TE = p.Results.TE;
T1 = p.Results.T1; % Nima : This is a vector now
RefCon = p.Results.RefCon;
Threshold = p.Results.Threshold;
Chi2Factor = p.Results.Chi2Factor;
nT2 = p.Results.nT2;
T2Range = p.Results.T2Range;
minangle = p.Results.MinRefAngle;
nangles = p.Results.nAngles;
% reg = p.Results.Reg;
% Nima: Change this value according to the new changes
FlipAngleMap = p.Results.FlipAngleMap;
faset=	~isempty(FlipAngleMap);	%p.Results.SetFlipAngle;

% nCores=p.Results.nCores;
% nCores= 4; % in case of error message about Cores;

savereg=strcmp(p.Results.Save_regparam,'yes');
saveNNLS=strcmp(p.Results.Save_NNLS_basis,'yes');

%==========================================================================
% Initialize all the data
%==========================================================================
% Start the clock
% tstart=tic;
% Find size of the data
[nrows,ncols,nslices,nechs] = size(image);
FlipAngleMap = zeros(nrows,ncols,nslices); % This fixes the parfor problem!
% Initialize map matrices
gdnmap=nan*ones(nrows,ncols,nslices);
ggmmap=nan*ones(nrows,ncols,nslices);
gvamap=nan*ones(nrows,ncols,nslices);
SNRmap=nan*ones(nrows,ncols,nslices);
FNRmap=nan*ones(nrows,ncols,nslices);
alphamap=nan*ones(nrows,ncols,nslices);
distributions=nan*ones(nrows,ncols,nslices,nT2);
mumap=nan*ones(nrows,ncols,nslices);
chi2map=nan*ones(nrows,ncols,nslices);

%
% decay_basis=nan*ones(nrows,ncols,nslices,nechs,nT2);
%
if saveNNLS
    decay_basis=nan*ones(nrows,ncols,nslices,nechs,nT2);
end

%==========================================================================
% Find the basis matrices for each flip angle
%==========================================================================
% Initialize parameters and variable for angle optimization
T2_times=logspace(log10(T2Range(1)),log10(T2Range(2)),nT2);
if faset==0
    flip_angles=linspace(minangle, 180, nangles);
    % basis_angles is a 1xnangles cell array that will contain the decay bases of each angle
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

try

parfor row = 1:nrows
    %row
    gdn=nan*ones(ncols,nslices);
    ggm=nan*ones(ncols,nslices);
    gva=nan*ones(ncols,nslices);
    SNR=nan*ones(ncols,nslices);
    FNR=nan*ones(ncols,nslices);
    alpha=nan*ones(ncols,nslices);
    dists=nan*ones(ncols,nslices,nT2);
    mus=nan*ones(ncols,nslices);
    chi2s=nan*ones(ncols,nslices);

    if saveNNLS
        basis_matrices=nan*ones(ncols,nslices,nechs,nT2);
    else
        basis_matrices=[];
    end
    for col = 1:ncols
        for slice = 1:nslices
            % Conditional loop to reject low signal pixels
            if image(row,col,slice,1)>=Threshold
                % Extract decay curve from the pixel
                decay_data = squeeze(image(row,col,slice,:));
				        obs_weigts = exp(decay_data/max(decay_data)); % Nima: set observation weights here
                if faset == 0
                    %======================================================
                    % Find optimum flip angle
                    %======================================================
                    alpha(col,slice) = Estimate_Alpha(basis_angles, nangles, decay_data, obs_weigts, flip_angles);
                else
                    alpha(col,slice) = FlipAngleMap(row,col,slice);
                end
                %======================================================
                % Fit basis matrix using alpha
                %======================================================
                basis_decay = Calc_basis_decay(nechs, nT2, alpha(col,slice), TE, T2_times, T1, RefCon);

                if saveNNLS
                    basis_matrices(col,slice,:,:) = basis_decay;
                end
                %==========================================================
                % Calculate T2 distribution and global parameters
                %==========================================================
                % Find distribution depending on regularization routine
                [T2_dis,mu,chi2] = Nima_UBC_NNLS(basis_decay, decay_data, obs_weigts, Chi2Factor);

                dists(col,slice,:) = T2_dis;
                mus(col,slice) = mu;
                chi2s(col,slice) = chi2;
                % Compute parameters of distribution
                gdn(col,slice) = sum(T2_dis);
                ggm(col,slice) = exp(dot(T2_dis,log(T2_times))/sum(T2_dis));
                gva(col,slice) = exp(sum((log(T2_times)-log(ggm(col,slice))).^2.*T2_dis')./sum(T2_dis)) - 1;
                decay_calc = basis_decay*T2_dis;
                residuals = decay_calc-decay_data;
                FNR(col,slice) = sum(T2_dis)/sqrt(var(residuals));
                SNR(col,slice) = max(decay_data)/sqrt(var(residuals));
            end
        end
    end
    % Record temporary maps into 3D outputs
    gdnmap(row,:,:) = gdn;
    ggmmap(row,:,:) = ggm;
    gvamap(row,:,:) = gva;
    SNRmap(row,:,:) = SNR;
    FNRmap(row,:,:) = FNR;
    alphamap(row,:,:) = alpha;
    distributions(row,:,:,:) = dists;
    mumap(row,:,:) = mus;
    chi2map(row,:,:)=chi2s;
    if saveNNLS
        decay_basis(row,:,:,:,:) = basis_matrices;
    end
end

catch err
    %if poolopenflag==1
    %    matlabpool close
    %end
    rethrow(err)
end

%if poolopenflag==1
%    matlabpool close
%end

% Assign outputs
maps.gdn = gdnmap;
maps.ggm = ggmmap;
maps.gva = gvamap;
maps.alpha = alphamap;
maps.FNR = FNRmap;
maps.SNR = SNRmap;
if savereg
    maps.mu=mumap;
    maps.chi2factor=chi2map;
end
%
% maps.NNLS_basis=decay_basis;
%
if saveNNLS
    maps.NNLS_basis=decay_basis;
end


end

function alpha = Estimate_Alpha(basis_angles, nangles, decay_data, obs_weigts, flip_angles)
  % Fit each basis and  find chi-squared
  chi2_alpha = zeros(1,nangles);
  for a=1:nangles
      [T2_dis_ls, ~, ~] = Nima_UBC_NNLS(basis_angles{a},decay_data, obs_weigts); % Nima: This is to be tested!
      decay_pred=basis_angles{a}*T2_dis_ls;
      chi2_alpha(a)=sum((decay_data-decay_pred).^2);
  end
  % Find the minimum chi-squared and the corresponding angle
  alpha_spline = flip_angles(1):0.001:flip_angles(end);
  temp_Chi2 = [chi2_alpha, flip(chi2_alpha(1:end-1))];
  chi2_spline = interp1([flip_angles, (360 - flip(flip_angles(1:end-1)))], temp_Chi2, alpha_spline,'spline'); % Nima : Here is a trick for ya!
  [~,index] = min(chi2_spline);
  alpha = alpha_spline(index);
  if alpha_spline(index) > 180.5 % Nima: maps everything to under 180; Threshold value is a bit arbitrary
    alpha = 360 - alpha_spline(index);
  end
end

function basis_decay = Calc_basis_decay(nechs, nT2, alpha, TE, T2_times, T1, RefCon)
  basis_decay=zeros(nechs,nT2);
  % Compute the NNLS basis over T2 space

  for x=1:nT2
      echo_amp = EPGdecaycurve(nechs, alpha, TE, T2_times(x), T1(x), RefCon); % Nima : T1 vector is used
      basis_decay(:,x) = echo_amp';
  end
end
