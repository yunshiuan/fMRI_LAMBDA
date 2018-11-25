%% Part 7-C(With SPM): Degree Centrality: Regress Degree Centrality by SVS(D~SVS+sex) (whole cerebral grey matter)
%% Note: Since this is a whole-brian voxels appraoch, one could actually conduct second-level by SPM (See Part7-C).
%% By SPM, tidy result table could be output (peak location, k size, p_FDR etc.).
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd(path);
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "cor2mni"

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044

%% Utalise SPM second-level analysis(Untidify-directly from batch)
matlabbatch{1}.spm.stats.factorial_design.dir = {'/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized_with_SPM'};
%%
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS023.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS024.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS025.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS026.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS027.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS028.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS029.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS030.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS032.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS033.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS034.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS035.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS036.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS037.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS039.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS041.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS042.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS043.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS045.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS052.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS054.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS055.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS056.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS057.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS058.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS059.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS060.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS061.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS062.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS064.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS067.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS072.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS078.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS081.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS082.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS083.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS084.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS087.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS092.nii,1'
                                                          '/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii/smoothed_Z_matrixS095.nii,1'
                                                          };
%%
%%
matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [0
                                                      0
                                                      0
                                                      0
                                                      0
                                                      1
                                                      1
                                                      0
                                                      0
                                                      1
                                                      0
                                                      0
                                                      0
                                                      1
                                                      0
                                                      1
                                                      0
                                                      0
                                                      0
                                                      1
                                                      0
                                                      0
                                                      1
                                                      1
                                                      1
                                                      0
                                                      0
                                                      0
                                                      0
                                                      0
                                                      1
                                                      0
                                                      1
                                                      1
                                                      1
                                                      0
                                                      0
                                                      0
                                                      1
                                                      0];
%%
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Gender';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
%%
matlabbatch{1}.spm.stats.factorial_design.cov(2).c = [-0.92673662
                                                      0.26053698
                                                      1.83935826
                                                      -0.10574955
                                                      0.32368983
                                                      -0.18153298
                                                      1.42254945
                                                      0.85417378
                                                      0.18475356
                                                      0.85417378
                                                      0.58893181
                                                      0.02055615
                                                      -1.62141798
                                                      0.50051782
                                                      -1.0022945
                                                      -0.44677495
                                                      -1.29302316
                                                      0.5510401
                                                      0.82891264
                                                      -0.56045008
                                                      -0.32046925
                                                      -0.19416355
                                                      -1.1035646
                                                      -0.96462833
                                                      2.02881682
                                                      -0.93936719
                                                      0.85417378
                                                      -0.21942469
                                                      -2.89710558
                                                      -1.24250088
                                                      -0.74990864
                                                      0.32368983
                                                      0.57630124
                                                      -0.37099153
                                                      1.06889348
                                                      0.98047949
                                                      -0.3836221
                                                      0.90469607
                                                      -0.74990864
                                                      -1.1035646];
%%
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Hedonism';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
%%
matlabbatch{1}.spm.stats.factorial_design.cov(3).c = [0.65266002
                                                      0.52945604
                                                      -0.74878521
                                                      -0.40997428
                                                      0.95040296
                                                      -0.17383332
                                                      1.91550078
                                                      0.37031757
                                                      -0.01469485
                                                      -0.89765668
                                                      -0.16869982
                                                      -0.54344525
                                                      1.67935982
                                                      -0.05576285
                                                      -0.67132438
                                                      0.7501965
                                                      2.05410525
                                                      0.63212602
                                                      1.09927443
                                                      -1.52907706
                                                      1.00687145
                                                      0.87340048
                                                      -0.58451325
                                                      -1.27753561
                                                      -1.63174704
                                                      0.04177364
                                                      1.2481459
                                                      1.21221141
                                                      -0.13276533
                                                      -1.3545381
                                                      1.49968736
                                                      -0.51264426
                                                      -0.2919038
                                                      -0.26623631
                                                      1.00173795
                                                      -1.3237371
                                                      -0.19436732
                                                      -1.08759615
                                                      -0.25596931
                                                      -0.77958621];
%%
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'Security';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/1_raw_residual_whole_brain/cerebral_grey_mask.nii,1'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'overall_positive';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'overall_negative';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'hedonism_positive';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'hedonism_negative';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 -1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'security_positive';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'security_negative';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 -1];
matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.005;
matlabbatch{4}.spm.stats.results.conspec.extent = 15;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.csv = true;

%% Change Result table name
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized_with_SPM')
path_old=pwd;
path_new=strcat(path_old,'/result_table');
filename={'Overall_positive.csv','Overall_negative.csv',...
    'Hedonism_positive.csv','Hedonism_negative.csv',...
    'Security_positive.csv','Security_negative.csv'};
for i=1:6
movefile(strcat('spm_2017Aug15_00',num2str(i),'.csv'), strcat('./result_table/',filename{i}));
end

cluster_threshold_beta(64,64,38,3.4375,4,8,3,'none',0,0,0.05,0.005,10000,'test');

