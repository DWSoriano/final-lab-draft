# Stability and boot

bootVoxel = function(resp_feat_df) {
  # Arguments:
  #   resp_feat_df: combined response and feature data frame.
  #                 output of combineRespFeat function from clean.R
  #
  # Returns:
  #   [1]: boot_df: data frame with lammbdas and correlations for all 100 
  #        bootstrap samples for the voxel 7 and voxel 20 models
  #   [2]: voxel7_feat_df: data frame with summary of features included across
  #        the 100 bootstrap samples for the voxel 7 model
  #   [3]: voxel20_feat100_df: data frame with summary of top 100 features
  #        based on coefficient values across the 100 bootstrap samples for 
  #        the voxel 20 model
  #   [4]: voxel7_feat_boot: data frame with summary of number of features 
  #        across bootstrap samples and number of orig features
 
  
# get voxel 7 feature matrix and response vector:
voxel_df_7 <-
  select(resp_feat_df, 
         names(resp_feat_df)[grepl(paste0("(voxel_resp", 7, "$|feat)"),
                                   names(resp_feat_df))])
names(voxel_df_7)[1] = "voxel_resp"
x_7 <- model.matrix(voxel_resp~.,voxel_df_7)[,-1]
y_7 <- voxel_df_7$voxel_resp

# get voxel 20 feature matrix and response vector:
voxel_df_20 <-
  select(resp_feat_df, 
         names(resp_feat_df)[grepl(paste0("(voxel_resp", 20, "$|feat)"),
                                   names(resp_feat_df))])
names(voxel_df_20)[1] = "voxel_resp"
x_20 <- model.matrix(voxel_resp~.,voxel_df_20)[,-1]
y_20 <- voxel_df_20$voxel_resp

#############################################################
# get voxel 7 features from training on original training set
#############################################################
# Fit model for voxel
model_7 <- glmnet(x_7[train,], y_7[train], alpha = 1)

# ESCV
set.seed(1)
escv_glm_7 <- escv.glmnet(x_7[train,], y_7[train], alpha = 1, 
                          nfolds = 10)
lambda_opt_7 <- escv_glm_7$lambda.escv

# save features with non-zero coefficients
coef_7 <- data.frame(
  matrix(predict(model_7, s = lambda_opt_7, type="coefficients")))
colnames(coef_7) <- "coef"
coef_7$feat_num <- as.numeric(row.names(coef_7)) - 1
coef_7$feat <- paste0("feat", coef_7$feat_num)
coef_7 <- select(coef_7, feat, coef)
coef_7 <- filter(coef_7, !feat == "feat0")
voxel7_orig_feats <- coef_7 %>%
  filter(coef != 0) %>%
  select(feat)

rm(model_7, escv_glm_7, lambda_opt_7, coef_7)

################################
# run on 100 bootstrap samples:
################################
for (i in 1:100) {
  
  # bootstrap sample of images in training set
  set.seed(i)
  train_boot <- sample(train, length(train), replace = T)
  
  # voxel 7
  #############
  # Lasso
  #############
  
  # Fit model for voxel
  model_7 <- glmnet(x_7[train_boot,], y_7[train_boot], alpha = 1)
  
  # ESCV
  set.seed(1)
  escv_glm_7 <- escv.glmnet(x_7[train_boot,], y_7[train_boot], alpha = 1, 
                            nfolds = 10)
  lambda_opt_7 <- escv_glm_7$lambda.escv
  # predict on validation set using optimal lambda
  pred_7 <- predict(model_7, s = lambda_opt_7, newx = x_7[validation,])
  # compute validation set cor
  cor_7 <- cor(pred_7, y_7[validation])
  mse_7 <- mean((pred_7 - y_7[validation])^2)
  
  # save features with non-zero coefficients
  coef_7 <- data.frame(
    matrix(predict(model_7, s = lambda_opt_7, type="coefficients")))
  colnames(coef_7) <- "coef"
  coef_7$feat_num <- as.numeric(row.names(coef_7)) - 1
  coef_7$feat <- paste0("feat", coef_7$feat_num)
  coef_7 <- select(coef_7, feat, coef)
  coef_7 <- filter(coef_7, !feat == "feat0")
  coef_7_non_zero <- coef_7 %>%
    filter(coef != 0) %>%
    select(feat)
  # record boot num
  if(nrow(coef_7_non_zero) != 0) {
  coef_7_non_zero$boot = i }
  
  # combine across bootstrap samples
  if(i == 1) {
    coef_7_non_zero_feat <- coef_7_non_zero
  }
  if (i != 1) {
    coef_7_non_zero_feat <- rbind(
      coef_7_non_zero_feat,
      coef_7_non_zero)
  }
  
  # voxel 20
  ########
  # Ridge
  ########
  
  # Fit model for voxel
  model_20 = glmnet(x_20[train_boot,], y_20[train_boot], alpha = 0)
  
  # aicc
  set.seed(1)
  aicc_glm_20 <- ic.glmnet(x_20[train_boot,], y_20[train_boot], alpha = 0, 
                           crit = "aicc")
  lambda_opt_20 <- aicc_glm_20$lambda
  # predict on validation set using optimal lambda
  pred_20 <- predict(model_20, s = lambda_opt_20, newx = x_20[validation,])
  # compute validation set cor
  cor_20 <- cor(pred_20, y_20[validation])
  mse_20 <- mean((pred_20 - y_20[validation])^2)
  
  # save features
  coef_20 <- data.frame(
    matrix(predict(model_20, s = lambda_opt_20, type="coefficients")))
  colnames(coef_20) <- "coef"
  coef_20$feat_num <- as.numeric(row.names(coef_20)) - 1
  coef_20$feat <- paste0("feat", coef_20$feat_num)
  coef_20 <- select(coef_20, feat, coef)
  coef_20 <- filter(coef_20, !feat == "feat0")
  coef_20_top_100 <- coef_20 %>%
    arrange(desc(abs(coef))) %>%
    select(feat) %>%
    slice(1:100)
  
  # combine across bootstrap samples
  if(i == 1) {
    coef_20_top_100_feat <- coef_20_top_100
  }
  if (i != 1) {
    coef_20_top_100_feat <- rbind(
      coef_20_top_100_feat,
      coef_20_top_100)
  }
  
  if(i == 1) {
    boot_df <- data.frame(boot = i,
                          lambda7 = lambda_opt_7,
                          cor7 = cor_7,
                          mse7 = mse_7, 
                          lambda20 = lambda_opt_20,
                          cor20 = cor_20,
                          mse20 = mse_20)
  }
  if (i != 1) {
    boot_df <- rbind(
      boot_df,
      data.frame(boot = i,
                 lambda7 = lambda_opt_7,
                 cor7 = cor_7,
                 mse7 = mse_7, 
                 lambda20 = lambda_opt_20,
                 cor20 = cor_20,
                 mse20 = mse_20)
    )
  }
}

# summarize features that appear across bootstrap samples for voxel 7
voxel7_feat_df <-
  coef_7_non_zero_feat %>%
  select(feat) %>%
  group_by(feat) %>%
  summarise(n_boot = n(),
            boot_perc = n_boot/100) %>%
  arrange(desc(n_boot))

# summarize number of features across boot samples and number of orig features
voxel7_feat_boot <-
  coef_7_non_zero_feat %>%
  mutate(orig = 
           ifelse(feat %in% voxel7_orig_feats$feat, 1, 0)) %>%
  select(boot, orig) %>%
  group_by(boot) %>%
  summarise(n_feat = n(),
            n_feat_orig = sum(orig))

for (i in 1:100) {
  if(!(i %in% voxel7_feat_boot$boot)) {
    voxel7_feat_boot <-
      rbind(voxel7_feat_boot,
            data.frame(boot = i, n_feat = 0, n_feat_orig = 0))
  }
}
voxel7_feat_boot <- arrange(voxel7_feat_boot, boot)

# summarize features that appear across bootstrap samples for voxel 20
voxel20_feat100_df <-
  coef_20_top_100_feat %>%
  select(feat) %>%
  group_by(feat) %>%
  summarise(n_boot = n()) %>%
  arrange(desc(n_boot))

return(list(boot_df, voxel7_feat_df, voxel20_feat100_df, voxel7_feat_boot))
}
