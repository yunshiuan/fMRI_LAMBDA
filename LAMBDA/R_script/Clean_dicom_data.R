library(pbapply)

temp=list.files(path = c("\\\\wcs-cifs\\wc\\edneuro\\LAMBDA"),
                pattern="*\\.dcm",
                recursive=T,
                full.names=T)
temp=grep(temp,pattern = "\\d+\\.dcm$",value = T)

pblapply(temp,
         function(x){
           file.remove(x)
         }
       )

# For testing if this works
# test=list.files(path = c("\\\\wcs-cifs\\wc\\edneuro\\LAMBDA\\DF1004\\dicoms\\00200.SRC__MPnRAGE_MOTION_CORRECTED_0"),
#                 pattern="*\\.dcm",
#                 recursive=T,
#                 full.names=T)
# test=grep(test=,pattern = "\\d+\\.dcm$",value = T)
# pblapply(test,
#          function(x){
#            file.remove(x)
#          }
# )

