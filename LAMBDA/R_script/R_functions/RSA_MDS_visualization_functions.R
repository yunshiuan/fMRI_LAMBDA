#Helper Functions for Representational Similarity Analysis(RSA)[MDS visualization]----------------------------------
#Note:
#Helper functions for:
#1st and 2nd order MDS visualization
#1st and 2nd order Shepard plots
library("grid")
#Constants=====================================)========
RATIO_GEOM_TEXT_FONT_SIZE=(0.352777778/1)
FILE_RSA_BASIC_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_basic_functions.R"
#Helper helper function========================)===================
source(FILE_RSA_BASIC_FUNCTIONS)
#Shepard plot===================================)------------------------------------------------------------
#Shepard plot geom function==========================================================================
geom_shepard=function(plot_title,text_font,df_text_corr,label){
  list(
    geom_point(size=2,aes(color=dis_measure_type,shape=dis_measure_type)),
    geom_line(aes(linetype=dis_measure_type,color=dis_measure_type)),
    
    scale_linetype_manual(name="Measure type",
                          values = PAR_TYPE_LINETYPE),
    scale_color_manual(name="Measure type",
                       values = PAR_MEASURE_COLOR),
    scale_shape_manual(name="Measure type",
                       values=PAR_TYPE_SHAPE),
    labs(title=plot_title,x="Dissimilarity",y="Distance \\\ Disparity"),
    theme_bw(),
    geom_text(data=df_text_corr,
              mapping = aes(x = -Inf,y = Inf,label = label),
              hjust=0,vjust=1,size=text_font)
  )
}
#Shepard for 1st order MDS (mat to df)-----------------------------------------------------------
#Read in MDS result from .mat file, and convert into a tidy data frame for further visualization.
#This handle .mat with 3 layers: (1)stress type, (2)age group, (3) VOI
shepard_mat_to_df_1_order=function(MDS_result,list_stress,list_age_group,list_VOI){
  df_shepard_collect=data.frame()
  #Save the coordinates (df_mds_shepard_coordinates)
  for(stress_type_index in 1:length(list_stress)){
    for(age_group_index in 1:length(list_age_group)){
      for(VOI_index in 1:length(list_VOI)){
        
        dissimilarities=vectorize_symmetric_matrix(MDS_result[[stress_type_index]][[age_group_index]][[VOI_index]][[INDEX_DISSIMILARITIES]],
                                                   return_indexes = T)
        #Retrieve the original cell index (in the matrix) of each vetorized element
        cell_indexes=dissimilarities$indexes
        point1_name=cell_indexes[,1]
        point2_name=cell_indexes[,2]
        
        dissimilarities=dissimilarities$vector
        distances=vectorize_symmetric_matrix(MDS_result[[stress_type_index]][[age_group_index]][[VOI_index]][[INDEX_DISTANCES]])
        disparities=vectorize_symmetric_matrix(MDS_result[[stress_type_index]][[age_group_index]][[VOI_index]][[INDEX_DISPARITIES]])
        stress=MDS_result[[stress_type_index]][[age_group_index]][[VOI_index]][[INDEX_STRESS]]
        
        df_shepard=data.frame(dissimilarities=dissimilarities,
                              distances=distances,
                              disparities=disparities,
                              stress=stress,
                              point1_name=point1_name,
                              point2_name=point2_name)
        
        df_shepard=
          df_shepard%>%
          mutate(MDS_stress_type=list_stress[[stress_type_index]],
                 age_group=list_age_group[[age_group_index]],
                 VOI_name=list_VOI[[VOI_index]])
        
        df_shepard_collect=
          df_shepard_collect%>%
          bind_rows(df_shepard)
      }
    }
  }
  return(df_shepard_collect)
}
#Shepard for 1st order MDS (visualization)--------------------------------------------------------
#Visualization
#This plot out Shepard plot with 2 options: (1)stress type, (2)VOI list types
#Facet: VOI~age
shepard_plot_1_order=function(df_shepard_collect,
                              list_interested_age_groups,
                              list_interested_stress_type,
                              list_interested_VOI_scope,
                              list_pdf_size_for_VOI_scope,
                              text_font_size_for_VOI_scope,
                              file_type_suffix=".pdf",
                              path_figures_output,
                              name_plot_title_base="Shepard plot \nCriterion = ",
                              name_file_name_base="Shepard_plot_agexVOI_",
                              save_output=FALSE,
                              peek_first_plot=TRUE,
                              return_plots=FALSE){
  #Collect the ggplots
  ggplot_collect=list()
  plot_index=1
  for(list_VOI_index in 1:length(list_interested_VOI_scope)){
    VOI_list_type=list_interested_VOI_scope[[list_VOI_index]]
    VOI_list_name=names(list_interested_VOI_scope)[list_VOI_index]
    
    pdf_width=list_pdf_size_for_VOI_scope[[list_VOI_index]][1]
    pdf_height=list_pdf_size_for_VOI_scope[[list_VOI_index]][2]
    #text font
    text_font=text_font_size_for_VOI_scope[list_VOI_index]
    
    for(stress_type_index in 1:length(list_interested_stress_type)){
      stress_type=list_interested_stress_type[[stress_type_index]]
      path_output=file.path(path_figures_output,stress_type,"Shepard_plot")
      dir.create(path_output,recursive = T)
      file_output=file.path(path_output,paste0(name_file_name_base,VOI_list_name,file_type_suffix))
      
      #Filter in  the interested scope
      df_shepard_visualization=
        df_shepard_collect%>%
        filter(MDS_stress_type %in% stress_type,
               age_group%in%list_interested_age_groups,
               VOI_name%in%VOI_list_type)
      #To make sure both axes share the same limits in the plot
      # limit_x_range=range(df_shepard_visualization$x)
      # limit_y_range=range(df_shepard_visualization$y)
      # limit_common_range=c(min(limit_x_range[1],limit_y_range[1]),
      #                      max(limit_x_range[2],limit_y_range[2]))
      
      #plot title
      plot_title=paste0(name_plot_title_base,stress_type)
      
      #Generate Correlation measures between dissimilarity and distance
      df_shepard_corr=
        df_shepard_visualization%>%
        select(dissimilarities,distances,VOI_name,age_group)%>%
        group_by(VOI_name,age_group)%>%
        summarise(Pearson_distance=cor(x=dissimilarities,y=distances,method = "pearson"),
                  Spearman_distance=cor(x=dissimilarities,y=distances,method = "spearman"))
      df_text_corr= #To be add onto each panel
        data.frame(label=paste0("Corr(Dissimilarity, Distance) \n",
                                round(df_shepard_corr$Pearson_distance,2),"(Pearson)\n",
                                round(df_shepard_corr$Spearman_distance,2),"(Spearman)"),
                   VOI_name=df_shepard_corr$VOI_name,
                   age_group=df_shepard_corr$age_group)
      
      #plot
      ggplot_collect[[plot_index]]=
        df_shepard_visualization%>%
          gather(dis_measure_type,dis_measure,matches("distances|disparities"))%>%
          ggplot(aes(x=dissimilarities,y=dis_measure))+
          # geom_point(size=2,aes(color=dis_measure_type,shape=dis_measure_type))+
          # geom_line(aes(linetype=dis_measure_type,color=dis_measure_type))+
          facet_wrap(VOI_name~age_group,ncol = length(list_interested_age_groups),
                     scales = "free")+
          geom_shepard(plot_title=plot_title,text_font=text_font,
                     df_text_corr=df_text_corr,label=label)+
          # scale_linetype_manual(name="Measure type",
          #                       values = PAR_TYPE_LINETYPE)+
          # scale_color_manual(name="Measure type",
          #                    values = PAR_MEASURE_COLOR)+
          # scale_shape_manual(name="Measure type",
          #                    values=PAR_TYPE_SHAPE)+
          # labs(title=plot_title,x="Dissimilarity",y="Distance \\\ Disparity")+
          # theme_bw()+
          # geom_text(data=df_text_corr,
          #           mapping = aes(x = -Inf,y = Inf,label = label),
          #           hjust=0,vjust=1,size=text_font)+
          theme(plot.title = element_text(hjust=0.5,size = 18,face = "bold"),
                axis.title = element_text(size=16),
                legend.title  = element_text(size=16),
                legend.text = element_text(size=14),
                # strip.background=element_rect(fill=NA, color=NA),
                strip.text = element_text(size=14),
                #Remove grid
                panel.grid.major =  element_blank(),
                panel.grid.minor =  element_blank())
          # ggrepel::geom_text_repel(aes(label=RSA_id)) #Mark the RSA ID
      
      #Save the plots to a file (false by default)
      if(save_output){
          ggplot_collect[[plot_index]]+
            ggsave(filename = file_output,
                   width = pdf_width,height = pdf_height)
        }
      plot_index=plot_index+1
      
    }
  }
  #Peek the first plot (true by default)
  if(peek_first_plot){
    print(ggplot_collect[[1]])
  }
  #Return the ggplot list (false by default)
  if(return_plots){
    plot_names=do.call(paste,
                  c(expand.grid(list_interested_stress_type,names(list_interested_VOI_scope)),
                    sep="_"))
    names(ggplot_collect)=plot_names
    return(ggplot_collect)
  }
}
#Shepard for 2nd order MDS (mat to df)---------------------------------------------------------
#Read in MDS result from .mat file, and convert into a tidy data frame for further visualization.
#This handle .mat with 2 layers: (1)stress type, (2)age group
shepard_mat_to_df_2_order=function(list_MDS_result,
                                   list_stress,
                                   list_age_group,
                                   list_ACC_version){
  df_shepard_collect=data.frame()
  for(acc_version_index in 1:length(list_ACC_version)){
    MDS_result=list_mds_result[[acc_version_index]]
    mds_version=list_ACC_version[acc_version_index]
    
    for(stress_type_index in 1:length(list_stress)){
      
      for(age_group_index in 1:length(list_age_group)){
  
          dissimilarities=vectorize_symmetric_matrix(MDS_result[[stress_type_index]][[age_group_index]][[INDEX_DISSIMILARITIES]],
                                                     return_indexes = T)
          #Retrieve the original cell index (in the matrix) of each vetorized element
          cell_indexes=dissimilarities$indexes
          point1_name=cell_indexes[,1]
          point2_name=cell_indexes[,2]
          
          dissimilarities=dissimilarities$vector
          distances=vectorize_symmetric_matrix(MDS_result[[stress_type_index]][[age_group_index]][[INDEX_DISTANCES]])
          disparities=vectorize_symmetric_matrix(MDS_result[[stress_type_index]][[age_group_index]][[INDEX_DISPARITIES]])
          stress=MDS_result[[stress_type_index]][[age_group_index]][[INDEX_STRESS]]
          
          df_shepard=data.frame(dissimilarities=dissimilarities,
                                distances=distances,
                                disparities=disparities,
                                stress=stress,
                                point1_name=point1_name,
                                point2_name=point2_name,
                                acc_version=mds_version,stringsAsFactors = F)
          
          df_shepard=
            df_shepard%>%
            mutate(MDS_stress_type=list_stress[[stress_type_index]],
                   age_group=list_age_group[[age_group_index]])
          
          df_shepard_collect=
            df_shepard_collect%>%
            bind_rows(df_shepard)
        }
      }
  }
  return(df_shepard_collect)
}

#Shepard for 2nd order MDS (visualization)--------------------------------------------------------
#Visualization
#This plot out Shepard plot with 2 options: (1)stress type, (2)VOI list types(currently length=1)
#Facet: VOI~age
shepard_plot_2_order=function(df_shepard_collect,
                              list_interested_age_groups,
                              list_interested_stress_type,
                              list_interested_VOI_scope,
                              list_ACC_version,
                              list_pdf_size_for_VOI_scope,
                              text_font_size_for_VOI_scope,
                              file_type_suffix=".pdf",
                              path_figures_output,
                              name_plot_title_base="Shepard plot \nCriterion = ",
                              name_file_name_base="Shepard_plot_agexVOI_",
                              save_output=FALSE,
                              peek_first_plot=TRUE,
                              return_plots=FALSE){
  #Collect the ggplots
  ggplot_collect=list()
  plot_index=1
  for(acc_version_index in 1:length(list_ACC_version)){
      name_acc_version=list_ACC_version[acc_version_index]
  
      df_shepard_visualization_ACC=
        df_shepard_collect%>%
        filter(acc_version%in%name_acc_version)
      
      for(list_VOI_index in 1:length(list_interested_VOI_scope)){
        VOI_list_type=list_interested_VOI_names[[list_VOI_index]]
        VOI_list_name=names(list_interested_VOI_names)[list_VOI_index]
        
        pdf_width=list_pdf_size_for_VOI_scope[[list_VOI_index]][1]
        pdf_height=list_pdf_size_for_VOI_scope[[list_VOI_index]][2]
        #text font
        text_font=text_font_size_for_VOI_scope[list_VOI_index]
        
        for(stress_type_index in 1:length(list_interested_stress_type)){
          stress_type=list_interested_stress_type[[stress_type_index]]
          path_output=file.path(path_figures_output,name_acc_version,stress_type,"Shepard_plot")
          dir.create(path_output,recursive = T)
          file_output=file.path(path_output,paste0(name_file_name_base,VOI_list_name,file_type_suffix))
          
          #Filter in  the interested scope
          df_shepard_visualization=
            df_shepard_visualization_ACC%>%
            filter(MDS_stress_type %in% stress_type,
                   age_group%in%list_interested_age_groups)
          #To make sure both axes share the same limits in the plot
          # limit_x_range=range(df_shepard_visualization$x)
          # limit_y_range=range(df_shepard_visualization$y)
          # limit_common_range=c(min(limit_x_range[1],limit_y_range[1]),
          #                      max(limit_x_range[2],limit_y_range[2]))
          
          #plot title
          plot_title=paste0(name_plot_title_base,stress_type)
          
          #Generate Correlation measures between dissimilarity and distance
          df_shepard_corr=
            df_shepard_visualization%>%
            select(dissimilarities,distances,age_group)%>%
            group_by(age_group)%>%
            summarise(Pearson_distance=cor(x=dissimilarities,y=distances,method = "pearson"),
                      Spearman_distance=cor(x=dissimilarities,y=distances,method = "spearman"))
          df_text_corr= #To be add onto each panel
            data.frame(label=paste0("Corr(Dissimilarity, Distance) \n",
                                    round(df_shepard_corr$Pearson_distance,2),"(Pearson)\n",
                                    round(df_shepard_corr$Spearman_distance,2),"(Spearman)"),
                       age_group=df_shepard_corr$age_group)
          
          #plot
          ggplot_collect[[plot_index]]=
            df_shepard_visualization%>%
            gather(dis_measure_type,dis_measure,matches("distances|disparities"))%>%
            ggplot(aes(x=dissimilarities,y=dis_measure))+
            geom_shepard(plot_title=plot_title,text_font=text_font,
                         df_text_corr=df_text_corr,label=label)+
            # geom_point(size=2,aes(color=dis_measure_type,shape=dis_measure_type))+
            # geom_line(aes(linetype=dis_measure_type,color=dis_measure_type))+
            facet_wrap(~age_group,ncol = length(list_interested_age_groups),
                       scales = "free")+
            # scale_linetype_manual(name="Measure type",
            #                       values = PAR_TYPE_LINETYPE)+
            # scale_color_manual(name="Measure type",
            #                    values = PAR_MEASURE_COLOR)+
            # scale_shape_manual(name="Measure type",
            #                    values=PAR_TYPE_SHAPE)+
            # labs(title=plot_title,x="Dissimilarity",y="Distance \\\ Disparity")+
            # theme_bw()+
            # geom_text(data=df_text_corr,
            #           mapping = aes(x = -Inf,y = Inf,label = label),
            #           hjust=0,vjust=1,size=text_font)+
            theme(plot.title = element_text(hjust=0.5,size = 18,face = "bold"),
                  axis.title = element_text(size=16),
                  legend.title  = element_text(size=16),
                  legend.text = element_text(size=14),
                  # strip.background=element_rect(fill=NA, color=NA),
                  strip.text = element_text(size=14),
                  #Remove grid
                  panel.grid.major =  element_blank(),
                  panel.grid.minor =  element_blank())
          # ggrepel::geom_text_repel(aes(label=RSA_id)) #Mark the RSA ID
          
          #Save the plots to a file (false by default)
          if(save_output){
            ggplot_collect[[plot_index]]+
              ggsave(filename = file_output,
                     width = pdf_width,height = pdf_height)
          }
          plot_index=plot_index+1
          
        }
      }
    }
  #Peek the first plot (true by default)
  if(peek_first_plot){
    print(ggplot_collect[[1]])
  }
  
  #Return the ggplot list (false by default)
  if(return_plots){
    plot_names=do.call(paste,
                       c(expand.grid(list_interested_stress_type,names(list_interested_VOI_scope)),
                         sep="_"))
    names(ggplot_collect)=plot_names
    return(ggplot_collect)
  }
}

#MDS plot ==============================)============================================
#MDS plot for 1st order MDS (mat to df)-----------------------------------------------------------
mds_mat_to_df_1_order=function(MDS_result,
                               list_stress,
                               list_age_group,
                               list_VOI){
  #And convert to a data frame (df_mds_coordinates)
  df_mds_coordinates=data.frame()
  #Save the coordinates (df_mds_coordinates)
  for(stress_type_index in 1:length(list_stress)){
    MDS_stress_type=list_stress[[stress_type_index]]
    for(age_group_index in 1:length(list_age_group)){
      age_group=list_age_group[[age_group_index]]
      for(VOI_index in 1:length(list_VOI)){
        VOI_name=list_VOI[[VOI_index]]
        df_coordinate=data.frame(MDS_result[[stress_type_index]][[age_group_index]]
                                 [[VOI_index]][[INDEX_COORDINATE]])
        names(df_coordinate)=c("x","y")
        
        df_coordinate=
          df_coordinate%>%
          mutate(MDS_stress_type=MDS_stress_type,
                 age_group=age_group,
                 VOI_name=VOI_name,
                 RSA_id=1:NUM_CONDITION,
                 RSA_format=rep(NAMES_RSA_FORMAT,each=6),
                 RSA_magnitude=rep(rep(NAMES_RSA_MAGNITUDE,each=2),times=3),
                 RSA_sign=rep(c(NAMES_RDS_SIGNS),times= 9))
        df_mds_coordinates=
          df_mds_coordinates%>%
          bind_rows(df_coordinate)
        print(paste0("Done: Stress type:",MDS_stress_type,
                     "; Age group: ",age_group,
                     "; VOI:",VOI_name))
      }
    }
  }
  return(df_mds_coordinates)
}
#MDS plot for 1st order MDS (visualization)-----------------------------------------------------------
mds_plot_1_order=function(df_mds_coordinates,
                          list_interested_age_groups,
                          list_interested_stress_type,
                          list_interested_VOI_names,
                          list_pdf_size,
                          file_type_suffix=".pdf",
                          par_mag_size,
                          path_figures_output,
                          name_plot_title_base="Conditions MDS for each neural RDM (18conditions) \nCriterion = ",
                          name_file_name_base="First_order_MDS_agexVOI_",
                          remove_plot_title=F,
                          name_plot_title_full=NULL,
                          with_kendall_tau=FALSE,
                          save_output=FALSE,
                          peek_first_plot=TRUE,
                          return_plots=FALSE){
  #Collect the ggplots
  ggplot_collect=list()
  plot_index=1
  
  for(list_VOI_index in 1:length(list_interested_VOI_names)){
    VOI_list_type=list_interested_VOI_names[[list_VOI_index]]
    VOI_list_name=names(list_interested_VOI_names)[list_VOI_index]
    pdf_width=list_pdf_size[[list_VOI_index]][1]
    pdf_height=list_pdf_size[[list_VOI_index]][2]
    
    for(stress_type_index in 1:length(list_interested_stress_type)){
      stress_type=list_interested_stress_type[[stress_type_index]]
      path_output=file.path(PATH_FIGURES_OUTPUT,stress_type,"MDS_plot")
      file_output=file.path(path_output,paste0(name_file_name_base,VOI_list_name,file_type_suffix))
      #Filter in  the interested scope
      df_mds_visualization=
        df_mds_coordinates%>%
        filter(MDS_stress_type %in% stress_type,
               age_group%in%list_interested_age_groups,
               VOI_name%in%VOI_list_type)
      #To make sure both axes share the same limits in the plot
      limit_x_range=range(df_mds_visualization$x)
      limit_y_range=range(df_mds_visualization$y)
      limit_common_range=c(min(limit_x_range[1],limit_y_range[1]),
                           max(limit_x_range[2],limit_y_range[2]))
      
      #Plot title
      if(remove_plot_title){
        plot_title=NULL
        #If plot title specified
      }else if(!is.null(name_plot_title_full)){
        plot_title=name_plot_title_full
      #If not specified, create the plot title by combinging the base with stress type 
      }else if(name_plot_title_base=="Conditions MDS for each neural RDM (18 conditions) \nCriterion = "){
        plot_title=paste0(name_plot_title_base,stress_type)
      # If stress type is not important, then ignore it.
      }else{
        plot_title=name_plot_title_base
      }
      
      # # Add Kendall tau to the plots
      # grob <- grobTree(textGrob(c("a","b"), x = 0.1,y = 0.9, hjust=0,
      #                           gp=gpar(col="red", fontsize=13, fontface="italic")))
      #plot
      ggplot_collect[[plot_index]]=
        df_mds_visualization%>%
          ggplot(aes(x=x,y=y))+
          geom_point(aes(color=RSA_format,size=RSA_magnitude,shape=RSA_sign),stroke=1.5)+ #Plot the format and mag.
          # facet_wrap(VOI_name~age_group,ncol = length(list_interested_age_groups))+
          facet_grid(VOI_name~age_group,switch = "y")+
          # geom_point(aes(shape=RSA_sign),color="black")+ #Overlay the sign
          scale_x_continuous(limits = limit_common_range)+
          scale_y_continuous(limits = limit_common_range)+
          scale_color_manual(name="Format",
                             values = PAR_FORMAT_COLOR)+
          scale_size_manual(values = par_mag_size,
                            name="Magnitude",labels = c("Near","Medium","Far"))+
          scale_shape_manual(name="Sign",values=PAR_SIGN_SHAPE,
                             labels=c("Right > Left","Right < Left"))+
          labs(title=plot_title)+
          theme_bw()+
          coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)+ #Fix the x/y axis as having same units
          theme(axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                plot.title = element_text(hjust=0.5,size = 14,face = "bold"),
                strip.background=element_rect(fill=NA, color=NA),
                strip.text = element_text(size=12),
                #Remove grid
                panel.grid.major =  element_blank(),
                panel.grid.minor =  element_blank(),
                #legend key space
                legend.key.size = unit(2, 'lines'),
                #text size
                legend.title = element_text(size=14),
                legend.text  = element_text(size=12)
                )
      # ggrepel::geom_text_repel(aes(label=c(45,43)),size=4)+ #Mark the sign
      # ggrepel::geom_text_repel(aes(label=RSA_id)) #Mark the RSA ID
      
      # Add the Kendall's tau if desired
      if(with_kendall_tau){
          ggplot_collect[[plot_index]]=
            ggplot_collect[[plot_index]]+
            # annotation_custom(grob)
            geom_text(mapping = aes(x = Inf,y = Inf,
                                    label = paste0("tau == ",
                                                   gsub(x = round(kendall_tau,2),
                                                        pattern = "^0",replacement = ""),
                                                   "~' '")),
                      hjust="inward",vjust=1.2,size=12*RATIO_GEOM_TEXT_FONT_SIZE,parse = T)
      }
      #Save the plots to a file (false by default)
      if(save_output){
        dir.create(path_output,recursive = T)
        ggplot_collect[[plot_index]]+
          ggsave(filename = file_output,
                 width = pdf_width,height = pdf_height)
      }
      plot_index=plot_index+1
    }
  }
  
  #Peek the first plot (true by default)
  if(peek_first_plot){
    print(ggplot_collect[[1]])
  }
  
  #Return the ggplot list (false by default)
  if(return_plots){
    plot_names=do.call(paste,
                       c(expand.grid(list_interested_stress_type,
                                     names(list_interested_VOI_names)),
                         sep="_"))
    names(ggplot_collect)=plot_names
    return(ggplot_collect)
  }
}
#MDS plot for 2nd order MDS (mat to df)-----------------------------------------------------------
mds_mat_to_df_2_order=function(list_MDS_result,
                               list_stress,
                               list_age_group,
                               list_ACC_version){
  
  #And convert to a data frame (df_mds_coordinates)
  df_mds_coordinates=data.frame() #Save the coordinates
  for(acc_version_index in 1:length(list_ACC_version)){
    mds_result=list_mds_result[[acc_version_index]]
    mds_version=NAMES_ACC_VERSION[acc_version_index]
    
    for(stress_type_index in 1:length(list_stress)){
      for(age_group_index in 1:length(list_age_group)){
        df_coordinate=data.frame(mds_result[[stress_type_index]][[age_group_index]][[INDEX_COORDINATE]])
        dot_names=gsub(mds_result[[stress_type_index]][[age_group_index]][[INDEX_DOT_NAME]],pattern = " ",replacement = "")
        names(df_coordinate)=c("x","y")
        
        df_coordinate=
          df_coordinate%>%
          mutate(MDS_stress_type=list_stress[[stress_type_index]],
                 age_group=list_age_group[[age_group_index]],
                 acc_version=mds_version,
                 RDM_name=dot_names)
        
        df_mds_coordinates=
          df_mds_coordinates%>%
          bind_rows(df_coordinate)
        print(paste0("Done: ",
                     "ACC version",acc_version_index,
                     "; Stress type:",stress_type_index,
                     "; Age group: ",age_group_index))
      }
    }
  } 
  return(df_mds_coordinates)
}
#MDS plot for 2nd order MDS (visualization)-----------------------------------------------------------
mds_plot_2_order=function(df_mds_coordinates,
                          list_interested_age_groups,
                          list_interested_stress_type,
                          list_interested_VOI_names,
                          list_interested_conceptual_models,
                          list_ACC_version,
                          list_pdf_size,
                          # list_pdf_size_for_VOI_scope,
                          # text_font_size_for_VOI_scope,
                          file_type_suffix=".pdf",
                          path_figures_output,
                          name_plot_title_base="Second order MDS for RDMs \nCriterion = ",
                          name_file_name_base="Second_order_MDS_",
                          remove_plot_title=F,
                          name_plot_title_full=NULL,
                          df_rdm_color,
                          save_output=FALSE,
                          peek_first_plot=TRUE,
                          return_plots=FALSE){
  #Collect the ggplots
  ggplot_collect=list()
  plot_index=1
  for(acc_version_index in 1:length(list_ACC_version)){
    
    name_acc_version=list_ACC_version[acc_version_index]
    
    df_mds_visualization_ACC=
      df_mds_coordinates%>%
      filter(acc_version%in%name_acc_version)
    
    for(list_VOI_index in 1:length(list_interested_VOI_names)){
      
      VOI_list_type=list_interested_VOI_names[[list_VOI_index]]
      VOI_list_name=names(list_interested_VOI_names)[list_VOI_index]
      
      #The interested RDMs for this iteration BASE MODELS + (with/without ACC x with/without DLPFC)
      list_RDM_interested=c(list_interested_conceptual_models,VOI_list_type) 
      
      #Include ACC if intrested
      if (name_acc_version ==  "with_ACC"){
        list_RDM_interested=c(list_RDM_interested,"ACC")
        #No need to include ACC if not intrested
      }
      
      pdf_width=list_pdf_size[[list_VOI_index]][1]
      pdf_height=list_pdf_size[[list_VOI_index]][2]
      
      for(stress_type_index in 1:length(list_interested_stress_type)){
        
        stress_type=list_interested_stress_type[[stress_type_index]]
        path_output=file.path(path_figures_output,name_acc_version,stress_type,"MDS_plot")
        dir.create(path_output,recursive = T)
        file_output=file.path(path_output,paste0(name_file_name_base,VOI_list_name,file_type_suffix))
        #Filter in  the interested scope of RDMs
        df_mds_visualization=
          df_mds_visualization_ACC%>%
          filter(MDS_stress_type %in% stress_type,
                 age_group %in% list_interested_age_groups,
                 RDM_name %in% list_RDM_interested)
        #To make sure both axes share the same limits in the plot
        limit_x_range=range(df_mds_visualization$x)
        limit_y_range=range(df_mds_visualization$y)
        limit_common_range=c(min(limit_x_range[1],limit_y_range[1]),
                             max(limit_x_range[2],limit_y_range[2]))
        #Plot title
        if(remove_plot_title){
          plot_title=NULL
          #If plot title specified
        }else if(!is.null(name_plot_title_full)){
          plot_title=name_plot_title_full
          #If not specified, create the plot title by combinging the base with stress type 
        }else if(name_plot_title_base=="Second order MDS for RDMs \nCriterion = "){
          plot_title=paste0(name_plot_title_base,stress_type)
          # If stress type is not important, then ignore it.
        }else{
          plot_title=name_plot_title_base
        }
        #RDM colors
        plot_RDM_color=
          df_rdm_color%>%
          filter(RDM_name%in%list_RDM_interested)%>%
          pull(RDM_color)
        
        #plot
        ggplot_collect[[plot_index]]=
          df_mds_visualization%>%
            ggplot(aes(x=x,y=y))+
            geom_point(aes(color=RDM_name),size=4.5)+ #Plot the format and mag.
            facet_wrap(~age_group,ncol = length(list_interested_age_groups))+
            coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)+ #Fix the x/y axis as having same units
            # geom_point(aes(shape=RSA_sign),color="black")+ #Overlay the sign
            scale_x_continuous(limits = limit_common_range)+
            scale_y_continuous(limits = limit_common_range)+
            scale_color_manual(name="RDM Name",
                               values =  plot_RDM_color)+
            labs(title=plot_title)+
            theme_bw()+
            theme(axis.title = element_blank(),
                  axis.ticks = element_blank(),
                  axis.text = element_blank(),
                  plot.title = element_text(hjust=0.5,size = 14,face = "bold"),
                  strip.background=element_rect(fill=NA, color=NA),
                  strip.text = element_text(size=14),
                  #Remove legend
                  legend.position = "none",
                  #Remove grid
                  panel.grid.major =  element_blank(),
                  panel.grid.minor =  element_blank())+
            # ggrepel::geom_text_repel(aes(label=c(45,43)),size=4)+ #Mark the sign
            ggrepel::geom_text_repel(aes(label=RDM_name_as_label),size=4,
                                     force = 10,
                                     box.padding = unit(0.35, "lines"), 
                                     point.padding = unit(0.35, "lines")#,direction="x"
                                     ) #Mark the RDM names
          
            #Save the plots to a file (false by default)
            if(save_output){
              ggplot_collect[[plot_index]]+
                ggsave(filename = file_output,
                       width = pdf_width,height = pdf_height)
            }
        
            # cat(plot_index,"\n acc: ",name_acc_version,
            #        "\n VOI_list_name: ",VOI_list_name,
            #        "\n stress_type: ",stress_type,
            #        "\n list_RDM_interested:" ,list_RDM_interested,
            #        "\n plot_RDM_color: ",plot_RDM_color,"\n")
            print(ggplot_collect[[plot_index]])
            plot_index=plot_index+1
         }
      }
  }
  #Peek the first plot (true by default)
  if(peek_first_plot){
    print(ggplot_collect[[1]])
  }
  
  #Return the ggplot list (false by default)
  if(return_plots){
    plot_names=do.call(paste,
                       c(expand.grid(list_interested_stress_type,names(list_interested_VOI_names),NAMES_ACC_VERSION),
                         sep="_"))
    names(ggplot_collect)=plot_names
    
    return(ggplot_collect)
  }
}
#LOO:MDS plot for 1st order MDS (mat to df)---------------------------------
mds_mat_to_df_1_order_LOO=function(MDS_result,
                                   list_stress,
                                   list_age_group,
                                   list_VOI){
  #And convert to a data frame (df_mds_coordinates)
  df_mds_coordinates=data.frame()
  #Save the coordinates (df_mds_coordinates)
  for(stress_type_index in 1:length(list_stress)){
    MDS_stress_type=list_stress[[stress_type_index]]
    for(age_group_index in 1:length(list_age_group)){
        list_sub_name_LOO=dimnames(MDS_result[[stress_type_index]][[age_group_index]])[[1]]
        age_group=list_age_group[[age_group_index]]
        for(sub_LOO_index in 1:length(list_sub_name_LOO)){
            name_sub_LOO=list_sub_name_LOO[sub_LOO_index]
            for(VOI_index in 1:length(list_VOI)){
              VOI_name=list_VOI[[VOI_index]]
              df_coordinate=data.frame(MDS_result[[stress_type_index]][[age_group_index]]
                                       [[sub_LOO_index]][[VOI_index]][[INDEX_COORDINATE]])
              names(df_coordinate)=c("x","y")
              
              df_coordinate=
                df_coordinate%>%
                mutate(MDS_stress_type=MDS_stress_type,
                       age_group=age_group,
                       VOI_name=list_VOI[[VOI_index]],
                       sub_name_LOO=name_sub_LOO,
                       RSA_id=1:NUM_CONDITION,
                       RSA_format=rep(NAMES_RSA_FORMAT,each=6),
                       RSA_magnitude=rep(rep(NAMES_RSA_MAGNITUDE,each=2),times=3),
                       RSA_sign=rep(c(NAMES_RDS_SIGNS),times= 9))
              df_mds_coordinates=
                df_mds_coordinates%>%
                bind_rows(df_coordinate)
              print(paste0("Done: Stress type:",MDS_stress_type,
                           "; Age group: ",age_group,
                           "; Sub name (LOO):",name_sub_LOO,
                           "; VOI:",VOI_name))
            }
        }
    }
  }
  return(df_mds_coordinates)
}
#Rotation and relfection the 1st order MDS plot====================)========================
#Rotation: find the rotation parameter - the axis radian--------------------------------------------------
#Find the axis in the 1st order MDS plot which indicates the dimension of magnitude.
#The axis is identified by maximizing the Kendall's tau-b (order projected on the axis, the real order definded by magnitude).
#Kendall's tau is chosen because magnitude is an ordinary variable (N,M,F) and it deals with ties better than Spearman's rho.
#Parameter:
#-cord_x: A numeric vector of the original x coordinate values.
#-cord_y: A numeric vector of the original y coordinate values.
#-magnitude: A numeric vector of the magnitude value of the points.
#Return:
#-axis_radian: The radian of the magnitude-axis identified. The vector points to the increasing magnitude direction.
#-kendall_tau: Spearman's R that indicates how obvious/dominant the magnitude-axis is.
mds_find_mag_axis=function(cord_x,cord_y,magnitude){
  
  #Shift the coordinates to the 1st quadrant to facilitate projeciton
  min_x=min(cord_x)
  min_y=min(cord_y)
  if(min_x<0){
    cord_x=cord_x-min(cord_x)
  }
  if(min_y<0){
    cord_y=cord_y-min(cord_y)
  }
  matrix_coordinates=rbind(cord_x,cord_y) # To be applied by the rotation matrix
  
  #Rotate the points for 180 (radian: 0~pi) degree and compute the according Kendall's tau
  #Rotating the points is equivalent to rotating the axis, yet easier (can use the rotation matrix)
  radian_search_space=seq(from=0,to = pi,length.out = 180)
  collect_kendall_tau=c()
  for (radian_index in 1:length(radian_search_space)){
    radian=radian_search_space[radian_index]#The radian at this iteration
    real_order=magnitude #The real order definded by magnitude
    #Rotation
    rotation_matrix=matrix(data = c(cos(radian),-sin(radian),
                                    sin(radian),cos(radian)),
                           ncol = 2,
                           byrow = T)
    cord_x_rotated=c(rotation_matrix[1,]%*% matrix_coordinates) #Only care about the x-axis value
    collect_kendall_tau[radian_index]=
      cor(x = cord_x_rotated,y = real_order,
          method = "kendal")
  }
  
  #Return the radian of the axis which maximize the Kendall's tau
  max_radian_index=which.max(abs(collect_kendall_tau))  
  max_radian=radian_search_space[max_radian_index]
  max_kendall_tau=collect_kendall_tau[max_radian_index]
  
  #Add pi if the axis is going to the decreasing magnitude direction
  if(max_kendall_tau<0){
    max_radian=max_radian+pi
    max_kendall_tau=-max_kendall_tau
  }
  return(list(axis_radian=max_radian,
              kendall_tau=max_kendall_tau))
}

#Rotation: rotate the MDS plot according to the axis radian-----------------------------------------------
#Rotate the mds coordinates by a given radian
#Parameter:
#-cord_x: A numeric vector of the original x coordinate values.
#-cord_y: A numeric vector of the original y coordinate values.
#-radian: The radian to rotate.
#Return:
#-cord_x_rotated: The rotated x coordinate values.
#-cord_y_rotated: The rotated y coordinate values.
mds_rotate=function(cord_x,cord_y,radian){
  coordinates=rbind(cord_x,cord_y) # To be applied by the rotation matrix
  rotation_matrix=matrix(data = c(cos(radian),-sin(radian),
                                  sin(radian),cos(radian)),
                         ncol = 2,
                         byrow = T)
  coordinates_rotated=rotation_matrix%*% coordinates
  
  return(list(cord_x_rotated=coordinates_rotated[1,],
              cord_y_rotated=coordinates_rotated[2,]))
}
#Reflection: identify if the MDS needs to be reflected or not------------------
#Decide whether the format axis of the MDS plot needs to be reflected,
#so that the format is aligned in the desired order.
#Note that reflection could only be done vertically in this specific case since
#the x-axis has already beed used as the magnitude axis 
#(with direction, pointing at the increasing magnitude)
#Parameter:
#-cord_x: A numeric vector of the original x coordinate values.
#-cord_y: A numeric vector of the original y coordinate values.
#-RSA_format: A character vector of the RSA format of each point.
#-RSA_order_desired: The desired top and the botton of the format vertical axis.
#Note that the position of CN is less important since it should be in the middle layer.
#Return:
#-need_reflection: Whether the plot needs to be reflected.
mds_format_need_reflection_or_not=function(cord_x,cord_y,RSA_format,
                                           RSA_order_desired=c("FF","LL")){
  df_mds_format=data.frame(cord_x,cord_y,RSA_format,stringsAsFactors = F)
  need_reflection=
    df_mds_format%>%
      filter(RSA_format!="CN")%>%
      group_by(RSA_format)%>%
      summarise(mean_cord_y=mean(cord_y))%>%
      arrange(desc(mean_cord_y))%>%
      pull(RSA_format)%>%
      as.character()%>%
      identical(RSA_order_desired)
  #No need if already in the desired order
  need_reflection=!need_reflection
  return(need_reflection)
}
#Reflection: reflect (mirror image) the MDS plot----------------------------
#Create the mirror image by reflecting the mds plot.
#Parameter:
#-cord_x: A numeric vector of the original x coordinate values.
#-cord_y: A numeric vector of the original y coordinate values.
#-mirror_vertical: Boolean value indicating whether to reflect vertically or not.
#-mirror_horizontal: Boolean value indicating whether to reflect horizontally or not.
#Return:
#-cord_x_reflected: The reflected x coordinate values.
#-cord_y_reflected: The reflected y coordinate values.
mds_reflect=function(cord_x,cord_y,
                     mirror_vertical,
                     mirror_horizontal){
  #Relfection
  if(mirror_horizontal){
    min_x=min(cord_x)
    max_x=max(cord_x)
    cord_x_reflected=scale_RSA(-cord_x,new_min =min_x,new_max = max_x)
  }else{
    cord_x_reflected=cord_x
  }
  
  if(mirror_vertical){
    min_y=min(cord_y)
    max_y=max(cord_y)
    cord_y_reflected=scale_RSA(-cord_y,new_min =min_y,new_max = max_y)
  }else{
    cord_y_reflected=cord_y
  }
  return(list(cord_x_reflected=cord_x_reflected,
              cord_y_reflected=cord_y_reflected))
}

#Rotation and relfection the 2nd order MDS plot====================)========================
#Rotation: find the rotation parameter - the axis radian--------------------------------------------------
#Find the axis in the 2nd order MDS plot which set the G-format on the left and the middle point of log-A and log-S on the right.
#Parameter:
#-cord_x: A numeric vector of the original x coordinate values.
#-cord_y: A numeric vector of the original y coordinate values.
#-RDM_name_as_label: A character vector of the names of the RDM points.
#Return:
#-axis_radian: The radian of the axis identified. The vector points at the format-to-magnitude direction.
mds_find_2nd_order_axis=function(cord_x,cord_y,RDM_name_as_label){
  
  #Shift the coordinates to the 1st quadrant to facilitate projection
  min_x=min(cord_x)
  min_y=min(cord_y)
  if(min_x<0){
    cord_x=cord_x-min(cord_x)
  }
  if(min_y<0){
    cord_y=cord_y-min(cord_y)
  }
  # matrix_coordinates=rbind(cord_x,cord_y) #To be applied by the rotation matrix
  df_coordinates=data.frame(cord_x=cord_x,
                            cord_y=cord_y,
                            RDM_name_as_label,
                            stringsAsFactors = F)
    
  #Identify the axis by finding the two points
  #Starting point
  starting_point=
    df_coordinates%>%
      filter(RDM_name_as_label=="G-format")%>%
      select(cord_x,cord_y)%>%
      c()%>%
      unlist()
  
  #Ending point
  #log-A
  log_A_point=
    df_coordinates%>%
    filter(RDM_name_as_label=="log-A")%>%
    select(cord_x,cord_y)%>%
    c()%>%
    unlist()
  #log-S
  log_S_point=
    df_coordinates%>%
    filter(RDM_name_as_label=="log-S")%>%
    select(cord_x,cord_y)%>%
    c()%>%
    unlist()
  
  #Ending point (in the middle)
  ending_point=colMeans(matrix(c(log_A_point,log_S_point),byrow = T,nrow = 2))
    
  #Compute the radian of the slope of the axis
  delta_x=ending_point[1]-starting_point[1]
  delta_y=ending_point[2]-starting_point[2]
  slope=unname(delta_y/delta_x)
  radian=atan(slope)
  
  #If the direction is opposite (correct direction: left: G-format, right: log-magnitude)
  if(delta_x<0){
    radian=radian+pi
  }
  #Note the negative sign, since the output is the radian to be rotated (clock-wise).
  return(list(axis_radian=-radian))
}

#Reflection: identify if the MDS needs to be reflected or not------------------
#Decide whether the format axis of the MDS plot needs to be reflected,
#so that the format is aligned in the desired order.
#Note that reflection could only be done vertically in this specific case since
#the x-axis has already beed used as the magnitude axis 
#(with direction, pointing at the increasing magnitude)
#Parameter:
#-cord_x: A numeric vector of the original x coordinate values.
#-cord_y: A numeric vector of the original y coordinate values.
#-RSA_format: A character vector of the RSA format of each point.
#-RSA_order_desired: The desired order of the log-A and log-S (from top to the button).
#Note that the position of CN is less important since it should be in the middle layer.
#Return:
#-need_reflection: Whether the plot needs to be reflected.
mds_2nd_need_reflection_or_not=function(cord_x,cord_y,RDM_names,
                                           RSA_order_desired=c("log-A","log-S")){
  df_mds_RDM=data.frame(cord_x,cord_y,RDM_names,stringsAsFactors = F)
  need_reflection=
    df_mds_RDM%>%
    filter(RDM_names=="log-A"|RDM_names=="log-S")%>%
    distinct()%>%
    arrange(desc(cord_y))%>%
    pull(RDM_names)%>%
    as.character()%>%
    identical(RSA_order_desired)
  #No need if already in the desired order
  need_reflection=!need_reflection
  return(need_reflection)
}
