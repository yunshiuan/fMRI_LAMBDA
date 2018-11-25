function [DCM] = DCM_specify_model(GLM_mat,interested_inputs,VOIs,TE,A_matrix,B_matrix,C_matrix)
% Specify a DCM.mat to be estimated which SPM recognizes.
% function [DCM] = DCM_specify_model(GLM_mat,interested_inputs, VOIs,TE,
%                            A_matrix,B_matrix,C_matrix,options,D_matrix)
%
% Inputs:
% GLM_mat   the file name of SPM.mat which stores the 1st level GLM results. This
%           detrmine the design matrix (onset/duration of each
%           stimuli). Note that covariates of no interests could or could
%           not be included, since their inclusion/exclusion depends solely
%           on the parameter - interested_inputs.
% interested_inputs
%           a vector with numeric indices which determines which inputs in the
%           GLM_mat(indeices correspond to the order in the contrast manager)
%           should be included in the DCM analysis.
% VOIs      The absolute paths of VOIs, which are the time series for DCM
%           to fit. VOI could be creatd via SPM built-in module.
% TE        The time echo. This could be found via spm_vol(‘img.nii’).
%           The TE value is stored in the V.desp.
% A_matrix  The intrinsic connectivity binary matrix. The size should be
%           n(VOIs) x n(VOIs).
% B_matrix The modulating input binary matrix. The size should be n(VOIs)
%           x n(VOIs) x n(interested_inputs).
% C_matrix The driving input binary matrix. The size should be
%          n(VOIs)x n(interested_inputs).
%
% Optional inputs:
% options a binary vector with length of 5. The order is as followed:
%         nonlinearity (0: bilinear; 1:non-linear),states per region (0: one state;
%         1: two states), stochastic effect (0: no, 1: yes), center
%         input (0: no, 1: yes), usging CSD data features (0:time series, 1:yes).
%         This vecter is set to zeros(1,5) by default.
% D_matrix non-linear input if thr nonlinearity option is switched on.
%
%
% Example:
% For a fully-specified models with 3 VOIs and 2 interested inputs:
% [DCM] = DCM_specify_model('my_SPM.mat',[1,2],{'VOI_1.mat','VOI_2.mat','VOI_3.mat'},
%                   0.04,ones(3,3),ones(3,3,2),ones(3,2))
% Yun-Shiuan Chuang - Feb 9, 2018
% Contact: yunshiuan.chuang@gmail.com

%% Default Argument check
%----------------------------------------------------------
if nargin < 7
    help DCM_specify_model; error('Necessary parameters are missing');
end
if  isempty(GLM_mat)|| isempty(interested_inputs)|| isempty(VOIs)||isempty(TE)||...
        isempty(A_matrix)||isempty(B_matrix)||isempty(C_matrix)
    help DCM_specify_model; error('Necessary parameters are missing');
end
if ~exist('options','var'); options = zeros(1,5); end
%% Load regions of interest
%--------------------------------------------------------------------------
n_VOI=size(VOIs,1);
if ~exist('D_matrix','var'); D_matrix=zeros(n_VOI,n_VOI,0);end

for v=1:n_VOI
    load(VOIs{v},'xY');
    DCM.xY(v) = xY;
end

DCM.n = length(DCM.xY);      % number of regions
DCM.v = length(DCM.xY(1).u); % number of time points

%% Time series
%--------------------------------------------------------------------------
load(GLM_mat); % Load in GLM result matrix.
DCM.Y.dt  = SPM.xY.RT;% TR of the EPI run
DCM.Y.X0  = DCM.xY(1).X0; % time-series of the covariates of the GLM (nVolumes x nConfounds))
for i = 1:DCM.n % Loop through all the VOIs
    DCM.Y.y(:,i)  = DCM.xY(i).u; %time course of ROI i (105 x 1)
    DCM.Y.name{i} = DCM.xY(i).name; %name of ROI i
end
DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v); %Error covariance constraints (cell 1 x nROIs)

%% Experimental inputs (per run)
%--------------------------------------------------------------------------
n_interested_inputs = size(interested_inputs,2); % Amount of interested inputs.
DCM.U.dt   =  SPM.Sess.U(1).dt; % TR(s)/nslices
DCM.U.name = [SPM.Sess.U(interested_inputs).name]; % Name of each interested inputs.
% Onset time
U.u=[];
for i=1:n_interested_inputs
    U.u=[U.u SPM.Sess.U(interested_inputs(i)).u(33:end,1)];  % Time course of a stimulus being on and off
end
DCM.U.u    = U.u;
%% DCM parameters and options
%--------------------------------------------------------------------------
DCM.delays =repmat(SPM.xY.RT/2,DCM.n,1);
DCM.TE=TE;
DCM.options.nonlinear  = options(1);
DCM.options.two_state  = options(2);
DCM.options.stochastic = options(3);
DCM.options.centre  = options(4);
DCM.options.induced  = options(5);

%% Connectivity matrices for model
%--------------------------------------------------------------------------
DCM.a = A_matrix; %Fully specified intrinsic connecitons (between 8 VOIs)
DCM.b = B_matrix; %Fully specified modulating inputs
DCM.c = C_matrix;
DCM.d = D_matrix;

end