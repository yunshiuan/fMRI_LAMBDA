%% Part 2-a - Second Level (Parametric): Collapse all age and set age as a covariate in the second levle GLM
spm_jobman('initcfg');

%% Constants
%Paths
PATH_ROOT='D:\Yun-Shiuan_LAMBDA';
PATH_TOOL_CODE='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code';
addpath(PATH_TOOL_CODE);% enable read_mixed_csv()
PATH_FIRST_GLM=fullfile(PATH_ROOT,'First_level_parametric_design_matrix_and_result','Trial_2','with_reslice');
PATH_SECOND_GLM=fullfile(PATH_ROOT,'Second_level_parametric_new','Trial_4','with_reslice');
%Files
FILE_VALID_RUN=fullfile(PATH_ROOT,'Run_inclusion_info','inclusive_runs_indexes_new_April_5.csv');
FILE_DEMOGRAPHIC=fullfile(PATH_ROOT,'demographic','tidy_demographic.csv');
FILE_COMMONE_MASK=fullfile(PATH_ROOT,'template','for_this_study','binary_overlapped_WM_GM_common_mask_n40_with_signal_with_reslice.nii,1');

%Parameters
%Note that T contrast are named from 002~014 (001 is an F contrast)
T_Parameters_name={'B0_LL' 'B0_LF' 'B0_FF'... 
                   'Dist_All'...
                   'Dist_LL' 'Dist_LF' 'Dist_FF'... 
                   't_LF-LL' 't_LL-LF' ...
                   't_FF-LL' 't_LL-FF' ...
                   't_FF-LF' 't_LF-FF'};
NUMERIC_CONTRAST_INDEX=pad(strtrim(cellstr((num2str([2:14]')))),4,'left','0');       
%% Derive subjects with valid runs
run_inclusion_index=read_mixed_csv_to_table(FILE_VALID_RUN);
subject_list=unique(run_inclusion_index.sub_id);
%% Create a common mask based on first level results fisrt, to ensure there's signal within the mask for each particiapants

%% Read in Gender and Age as covariates
table_demopraphic=read_mixed_csv_to_table(FILE_DEMOGRAPHIC);

table_demopraphic=sortrows(table_demopraphic,'sub_id','ascend'); % Need to sort by subID to make sure the order is ascending
valid_Gender_Dummy=str2double(table2array(table_demopraphic(ismember(lower(table_demopraphic.sub_id),subject_list),'gender_dummy')));
valid_Age_Precise=str2double(table2array(table_demopraphic(ismember(lower(table_demopraphic.sub_id),subject_list),'scan_age_precise')));

%% Second Level Matrix: T_Parameter(from first level)~ B0 + Age + Gender
       
for con=1:numel(T_Parameters_name) %NOTE: con=1 refers to spmT_0002
    clear matlabbatch
    path_second_result = fullfile(PATH_SECOND_GLM,T_Parameters_name{con});
    mkdir(path_second_result);
    
    matlabbatch{1}.spm.stats.factorial_design.dir = {path_second_result};    
    first_T_Contrast=fullfile(PATH_FIRST_GLM,subject_list,['spmT_',NUMERIC_CONTRAST_INDEX{con},'.nii,1']);
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = first_T_Contrast;
    % Covariate
    % Age
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = valid_Age_Precise;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Age';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    % Gender
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = valid_Gender_Dummy;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Gender';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
    % Other parameters
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {FILE_COMMONE_MASK};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % Estimation------------------------------
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    % Contrast Manager------------------------
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = ['pos_' T_Parameters_name{con}];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = ['neg_' T_Parameters_name{con}];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'pos_Age';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'neg_Age';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'pos_Gender';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 1];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'neg_Gender';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 -1];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1;
    save(fullfile(path_second_result,'second_glm_batch.mat'),'matlabbatch')
end

%% Run Second Level
for con=1:numel(T_Parameters_name)
    path_second_result = fullfile(PATH_SECOND_GLM,T_Parameters_name{con});
    load(fullfile(path_second_result,'second_glm_batch.mat'));
    spm_jobman('run',matlabbatch);
end

%% Output Result Table
for first_con=1:numel(T_Parameters_name) 
    clear matlabbatch
    spm_File=fullfile(PATH_SECOND_GLM,T_Parameters_name{first_con},'SPM.mat');
    matlabbatch{1}.spm.stats.results.spmmat = {spm_File};

    for second_con=1:6
        threshold=(second_con<=2)*(0.005)+(second_con>=3)*(0.05); 
        % Set a higher threshold for age and gender, 
        % as they're exploratoty.
        matlabbatch{1}.spm.stats.results.conspec(second_con).titlestr = '';
        matlabbatch{1}.spm.stats.results.conspec(second_con).contrasts = second_con;
        matlabbatch{1}.spm.stats.results.conspec(second_con).threshdesc = 'none';
        matlabbatch{1}.spm.stats.results.conspec(second_con).thresh = threshold;
        matlabbatch{1}.spm.stats.results.conspec(second_con).extent = 15;
        matlabbatch{1}.spm.stats.results.conspec(second_con).conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec(second_con).mask.none = 1;
    end
    
    matlabbatch{1}.spm.stats.results.units = 1;
    matlabbatch{1}.spm.stats.results.export{1}.csv = true;

    spm_jobman('run',matlabbatch);
    
    
    % Rename and move the csv to result_table folder
    path_result_table=fullfile(PATH_SECOND_GLM,T_Parameters_name{first_con},'Result_table');
    mkdir(path_result_table);
    csv_old_name=dir(fullfile(PATH_SECOND_GLM,T_Parameters_name{first_con},'spm*.csv'));
    csv_old_name={csv_old_name.name}';
    csv_new_name={['pos_' T_Parameters_name{first_con},'_p005_k15'] ['neg_' T_Parameters_name{first_con},'_p005_k15'],...
                  'pos_Age_p05_k15' 'neg_Age_p05_k15',...
                  'pos_Gender_p05_k15' 'neg_Gender_p05_k15'};    
    csv_old_full_name=fullfile(PATH_SECOND_GLM,T_Parameters_name{first_con},csv_old_name);
    csv_new_full_name=fullfile(path_result_table,strcat(csv_new_name','.csv'));

    for csv=1:6
        movefile(csv_old_full_name{csv,1},csv_new_full_name{csv,1});
    end
end

