Note:
    (1)For scripts that are not in the 'For_adult_analysis' folder,
       they are compatible for both the adults and the children with adjustion of constants,
       or, they are analysing by the combination of both the adult and child data(e.g., RSA analysis).
    (2)For scripts in the 'For_adult_analysis' folder,
       they are designed exclusively for the adult (usually due to distinct raw data type).
    (3)helper_function:
        Helper functions designed for higher-level analysis (usually used along with other package).
    (4)tool_code
        Helper functions designed for lower-level analysis (usually used for basic I/O control).
For each project
    (1)RSA: The RRSA project.
    Pipeline order:
    As in trial 14 of the RSA. 2018/4/18:
         (1)RSA_GLM_estimate_t_value_RSA_ID_discrete_18_ID
         (2)RSA_VOI_Extraction_t_value_RSA_ID_discrete_18_ID
         (3)RSA_fMRIDataPreparation_t_value_no_run_info_RSA_ID_discrete_18_id_include_adult
         (4)RSA_1st_order_analysis_t_value_discrete_18_id_include_adult
         (5)RSA_2nd_order_analysis_t_value_discrete_18_include_adult_mds
         (6)RSA_2nd_order_analysis_t_value_discrete_18_include_adult_relatedness_test
         (7)RSA_2nd_order_analysis_t_value_discrete_18_include_adult_hybrid_base_models
                
