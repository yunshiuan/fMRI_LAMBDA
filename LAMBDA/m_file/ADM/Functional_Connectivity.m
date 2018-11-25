%% Functional Connectivity

%% Step1: Preprocessing
%Motion correction : reslice needed instead of only estimate (regardless of which approach,specific to functoinal connectivity)

%Approach one(Ziling's approach) (Fox, 2006)
% 1 Traditional Preprocessing
%- Native space (specific to each person)
%- So no nomarlization needed (GLM: will do nomarlizaiton in the first level)
%- Neither do smooth or coregistration

% 2 Remove Noises: Detrend and Bandpass (Use "Rest" software made by China) 
%- Detrend: Remove the increasing fMRI signal during the task
%- Remove Noise by Bandpass Filter: Only attain High(0.08)/ Low(0.009)
%  Freq signals(will also remove the effects of motion)
%   the two functons: rest_detrend; rest_bandpass

% 3 Partial out to remain only the Residuals (Use SPM "First Level" function)
% If Task : Should remove the HRF shape--> only remain the signals altered by connectivity among
%  regions
% If Resting: Don't need to remove

% Only remain the grey matters(Regardless of Task or Resting)
% white matters and CSF serve as the covariates which measure the whole
% brain noises(e.g., breathing)--> Use "First level" to partial out the
% white matters and CSF signals, only left the "residuals" to output; 
% If it is task, then should also partial out the 15 conditions, outcomes, and motions' HRF 
%addpath('/bml/Data/Bank2/Scripts/Matlab/imaging/')
%help extract_voxel_values_ADMEmod (Use to extract CSF and white matter voxel-wise averaged 218*5 images)
%SPM model estimate: Write Residual=T (output: Res_nii: 218*5 slides)

%(Optional: Deal with motions)
% 4-a Scrubbing
% Remove the slide which motions are severe
% 4-b ICA
% Remove motions' conponent

% First Level---------------------------
% Still use native space
% ROI: Backwards deformation (Convert MNI ROI back to Native Space)
% (Each Voxel contain a residual value)
% Nomarlize before second level
% Correlation Computation: ROIs(seeds, e.g., striatum) to other  voxels (e.g., cortex)
% and then fisher Z-transformation (each voxel has a z value linked to the seed)
% Each seed has one z map which could be further transformed to nii.
% Second Level--------------------------

%Approah two========================================
% 1 Generally same as the Approach one
%   Only add back normalize and smooth during preprocessing (Convert to MNI space)
% 2 Remove Noises: Detrend and Bandpass (Use "Rest" software made by China)
% 3 Partial out to remain only the Residuals (Use SPM "First Level" function)
% 4 (Optional: Deal with motions)

% First Level---------------------------
% Use MNI space
% (Each Voxel contain a residual value)

%Approach three=====================================

%% Degree centrality(Buckner,2009)==================================
%-Graph Thoery: Each node (voxel) link to other nodes by edges 
%-How many voxles link to the voxel? -- Criteria needed

%Preprocessing:
% The same (Approach 2)
%Correlation Matrix (nxn matrix, each person has a unique nxn matrix): 
% (Only specify those voxles in the grey matter map)

% Residual of the voxel correlates with all the other
% voxles' correlations (total : n voxels, each contains 218*5 images)
% Convert Correlation Matrix(nxn) to 1(pass the criteria)/ 0 matrix:
%Creteria types:
% > r = 0.25
% R > PR90(across the person's voxels)
% Dealing with Negative Correlations: --> Ignore or not ?

%Degree Centrality Computation: The sum of the passed-criteria connectivities 
% (each person has an image with n centralities)
% Z-transformation: scaling the n centralities across the person's voxels
% Output as nii.(SPM_write_vol) file, each voxel represents as a centrality

%Second Level:
% Group sample t-test: overall group hub
% Regression: which voxels correlate with SVS
% Or predefined ROI: Correlate centralityies in the ROI with SVS

%Interpretation:

%% Step2: First Level
%% Step3: Second Level