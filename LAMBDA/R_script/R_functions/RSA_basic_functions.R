#RSA basic finctions===============================================================
#The script contains general basic functions that is frequently used in multiple specific steps in RSA.
#These functions are used in other helper function script (e.g,RSA_conceptual_models_functions.R)
#Function to vectorize a RDM--------------------------------------------------
#Only retain the upper part of a symmetric matrix, and output as a vector
#INPUE:
#-matrix: the input matrix
#-return_indexes: If true, then return a matrix containing the cell indexes along with the vectorized matrix.
#VALUE:
#-vector.
#-$vector: the vectorized matrix.
#-$indexes: the cell index of each element in the output vector.
vectorize_symmetric_matrix=function(matrix,return_indexes=F){
  n_row=dim(matrix)[1]
  n_col=dim(matrix)[2]
  vector=c()#The vector to store the output
  if(!return_indexes){
    for(row_index in 1:(n_row-1)){
      for(col_index in (row_index+1):(n_col)){
        vector=c(vector,matrix[row_index,col_index])
      }
    }
    return(vector)
  }else{
    cell_indexes=matrix(ncol = 2,nrow = ((n_row)*(n_col-1))/2)
    r_index_for_cell_indexes=1
    for(row_index in 1:(n_row-1)){
      for(col_index in (row_index+1):(n_col)){
        vector=c(vector,matrix[row_index,col_index])
        cell_indexes[r_index_for_cell_indexes,]=c(row_index,col_index)
        r_index_for_cell_indexes=r_index_for_cell_indexes+1;
      }
    }
    return(list(vector=vector,indexes=cell_indexes))
  }
}
#Function to scale RDM with any unit to the scale of a given range---------------------------
#Scale and shift (linear transformation) the input variable to the desired range
#Parameter
#variable: could be a vector, a matrix , or an dataframe
#Value
#scaled_variable: the rescaled variable
scale_RSA=function(variable,new_min,new_max){
  
    old_min=min(variable)
    old_max=max(variable)
    old_range=old_max-old_min
    new_range=new_max-new_min
    
    scaling=new_range/old_range
    
    if(old_range!=0){
      scaled_variable=(variable*scaling)
      shift=new_min-min(scaled_variable)
      scaled_variable=scaled_variable+shift
      #Special case if all the inputs share the same value (i.e., range==0)
    }else{
      shift=new_min-old_min
      scaled_variable=(variable+shift)
    }
    return(scaled_variable)
  }
#Function to read in average neural RDMs----------------------------------------------------------------------
#Extract the average neual RDMs values from a list of average neural RDMs read in from matlab file
#Parameter:
#-list_collect_RDMs_neural: the list of average neural RDMs read in from matlab file
#-list_interested_neural_RDMs: the list of interested brain region to read in

RSA_read_avg_nerual_RDMS=function(list_collect_RDMs_neural,list_interested_neural_RDMs){
  list_collect_neural_RDMs=list()
  list_collect_neural_names=c()
  RDM_index=1
  for (age_group_index in 1:length(list_collect_RDMs_neural)){
    #Name the dim names of the array by corresponding VOI names to facilitate indexing in the inner loop
    VOI_names_in_mat=unlist(list_collect_RDMs_neural[[age_group_index ]]["name",,])
    dimnames(list_collect_RDMs_neural[[age_group_index]])[[2]]=VOI_names_in_mat
    age_group_name=names(list_collect_RDMs_neural)[age_group_index]
    
    for(VOI_index in 1:length(list_interested_neural_RDMs)){
      VOI_name_interested=list_interested_neural_RDMs[VOI_index]
      RDM=list_collect_RDMs_neural[[age_group_index]]["RDM",VOI_name_interested,1]
      list_collect_neural_RDMs[[RDM_index]]=RDM
      list_collect_neural_names[RDM_index]=paste0(age_group_name,"_",VOI_name_interested)
      RDM_index=RDM_index+1
    }
  }
  names(list_collect_neural_RDMs)=list_collect_neural_names
return(list_collect_neural_RDMs)
}

#Function to read in individual neural RDMs----------------------------------------------------------------------
#Extract the individual neual RDMs values from a list of average neural RDMs read in from matlab file
#Parameter:
#-list_collect_RDMs_neural: the list of average neural RDMs read in from matlab file
#-list_interested_neural_RDMs: the list of interested brain region to read in

RSA_read_ind_nerual_RDMS=function(list_collect_RDMs_neural,list_interested_neural_RDMs){
  list_collect_neural_RDMs=list()
  list_collect_neural_names=c()
  RDM_index=1
  num_subject=length(list_collect_RDMs_neural[[1]][1,1,])
  for (sub_id_index in 1:num_subject){
    #Name the dim names of the array by corresponding VOI names to facilitate indexing in the inner loop
    VOI_names_in_mat=unlist(list_collect_RDMs_neural[[1]]["name",,sub_id_index])
    VOI_names=str_extract(string = VOI_names_in_mat,pattern = "^(R|L)-(V1|IPS|DLPFC)")
    dimnames(list_collect_RDMs_neural[[1]])[[2]]=VOI_names
    sub_name=unique(str_extract(string = VOI_names_in_mat,pattern = "(?<=(R|L)-(V1|IPS|DLPFC)).*"))
    
    for(VOI_index in 1:length(list_interested_neural_RDMs)){
        VOI_name_interested=list_interested_neural_RDMs[VOI_index]
        RDM=list_collect_RDMs_neural[[1]]["RDM",VOI_name_interested,sub_id_index]
        list_collect_neural_RDMs[[RDM_index]]=RDM
        list_collect_neural_names[RDM_index]=paste0(VOI_name_interested,"_",sub_name)
        RDM_index=RDM_index+1
    }
  }
  names(list_collect_neural_RDMs)=list_collect_neural_names
  return(list_collect_neural_RDMs)
}