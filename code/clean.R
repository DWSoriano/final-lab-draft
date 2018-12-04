# clean and combine data:

combineRespFeat = function(resp_dat, fit_feat) {
  # Arguments:
  #   resp_dat: input matrix
  #   fit_feat: response variable
  #
  # Returns:
  #   Combined and cleaned response and features dataset
  
# clean response:
resp_dat <- data.frame(resp_dat)
names(resp_dat) <- paste("voxel_resp", 1:ncol(resp_dat), sep = "")

# clean features:
fit_feat <- data.frame(fit_feat)
names(fit_feat) <- paste("feat", 1:ncol(fit_feat), sep = "")

# Combine responses and features
resp_feat_df <- cbind(resp_dat, fit_feat)

return(resp_feat_df)
}
