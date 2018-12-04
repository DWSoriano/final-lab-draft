# Find optimal lambda

lambdaOpt = function(x, y, criteria, alph) {
  # Arguments:
  #   x: input matrix
  #   y: response variable
  #   criteria: model selection technique. 
  #             Options: "cv", "escv", "aic", "aicc", "bic"
  #   alph: The elasticnet mixing parameter, with 0 <= alph <= 1.
  #         alph = 1: lasso penalty, alph = 0: ridge penalty
  #
  # Returns:
  #   [1] optimal lambda value
  #   [2] validation set correlation using optimal lambda value
  
  # Fit model for voxel
  model <- glmnet(x[train,], y[train], alpha = alph)
  
  if(criteria == "cv") {
    # CV
    set.seed(1)
    cv_glm <- cv.glmnet(x[train,], y[train], alpha = alph, nfolds = 10)
    lambda_opt <- cv_glm$lambda.min
  } else if(criteria == "escv"){
    # ESCV
    set.seed(1)
    escv_glm <- escv.glmnet(x[train,], y[train], alpha = alph, nfolds = 10)
    lambda_opt <- escv_glm$lambda.escv
  } else if(criteria == "aic"){
    # AIC
    set.seed(1)
    aic_glm <- ic.glmnet(x[train,], y[train], alpha = alph, crit = "aic")
    lambda_opt <- aic_glm$lambda
  } else if(criteria == "aicc"){
    # AICc
    set.seed(1)
    aicc_glm <- ic.glmnet(x[train,], y[train], alpha = alph, crit = "aicc")
    lambda_opt <- aicc_glm$lambda
  } else if(criteria == "bic"){
    # BIC
    set.seed(1)
    bic_glm <- ic.glmnet(x[train,], y[train], alpha = alph, crit = "bic")
    lambda_opt <- bic_glm$lambda
  }
  
  # predict on validation set using optimal lambda
  pred <- predict(model, s = lambda_opt, newx = x[validation,])
  # compute validation set cor
  cor <- cor(pred, y[validation])
  return(c(lambda_opt, cor))
}

bestVoxel = function(resp_feat_df) {
  # Arguments:
  #   resp_feat_df: combined response and feature data.
  #                 Output of combineRespFeat function from clean.R
  #
  # Returns:
  #   [1]: optimal model for each voxel for each model selection technique and
  #        model type
  #   [2]: optimal model for each voxel based on correlation of observed and
  #        predicted responses.
  
  # set number of cores
  ncores <- 4
  registerDoParallel(ncores)
  
  optimal_lambda <- 
    foreach(i = 1:ncol(resp_dat), .combine = rbind) %dopar% {
      
      rm(list=ls(pattern="lasso|ridge|enet"))
      # Get df with response for voxel i:
      voxel_df <-
        select(resp_feat_df, 
               names(resp_feat_df)[grepl(paste0("(voxel_resp", i, "$|feat)"),
                                         names(resp_feat_df))])
      names(voxel_df)[1] = "voxel_resp"
      x <- model.matrix(voxel_resp~.,voxel_df)[,-1]
      y <- voxel_df$voxel_resp
      
      ########
      # Lasso
      ########
      # CV
      cv_lasso <- lambdaOpt(x = x, y = y, criteria = "cv", alph = 1)
      # ESCV
      escv_lasso <- lambdaOpt(x = x, y = y, criteria = "escv", alph = 1)
      # AIC
      aic_lasso <- lambdaOpt(x = x, y = y, criteria = "aic", alph = 1)
      # AICc
      aicc_lasso <- lambdaOpt(x = x, y = y, criteria = "aicc", alph = 1)
      # BIC
      bic_lasso <- lambdaOpt(x = x, y = y, criteria = "bic", alph = 1)
      
      ########
      # Ridge
      ########
      # CV
      cv_ridge <- lambdaOpt(x = x, y = y, criteria = "cv", alph = 0)
      # ESCV
      escv_ridge <- lambdaOpt(x = x, y = y, criteria = "escv", alph = 0)
      # AIC
      aic_ridge <- lambdaOpt(x = x, y = y, criteria = "aic", alph = 0)
      # AICc
      aicc_ridge <- lambdaOpt(x = x, y = y, criteria = "aicc", alph = 0)
      # BIC
      bic_ridge <- lambdaOpt(x = x, y = y, criteria = "bic", alph = 0)
      
      #############
      # Elastic Net
      #############
      # CV
      cv_enet <- lambdaOpt(x = x, y = y, criteria = "cv", alph = 0.5)
      # ESCV
      escv_enet <- lambdaOpt(x = x, y = y, criteria = "escv", alph = 0.5)
      # AIC
      aic_enet <- lambdaOpt(x = x, y = y, criteria = "aic", alph = 0.5)
      # AICc
      aicc_enet <- lambdaOpt(x = x, y = y, criteria = "aicc", alph = 0.5)
      # BIC
      bic_enet <- lambdaOpt(x = x, y = y, criteria = "bic", alph = 0.5)
      
      rm(voxel_df)
      
      data.frame(voxel = i,
                 # Lasso
                 cv_lasso_lam = cv_lasso[1],
                 escv_lasso_lam = escv_lasso[1],
                 aic_lasso_lam = aic_lasso[1],
                 aicc_lasso_lam = aicc_lasso[1],
                 bic_lasso_lam = bic_lasso[1],
                 cv_lasso_cor = cv_lasso[2],
                 escv_lasso_cor = escv_lasso[2],
                 aic_lasso_cor = aic_lasso[2],
                 aicc_lasso_cor = aicc_lasso[2],
                 bic_lasso_cor = bic_lasso[2],
                 # Ridge
                 cv_ridge_lam = cv_ridge[1],
                 escv_ridge_lam = escv_ridge[1],
                 aic_ridge_lam = aic_ridge[1],
                 aicc_ridge_lam = aicc_ridge[1],
                 bic_ridge_lam = bic_ridge[1],
                 cv_ridge_cor = cv_ridge[2],
                 escv_ridge_cor = escv_ridge[2],
                 aic_ridge_cor = aic_ridge[2],
                 aicc_ridge_cor = aicc_ridge[2],
                 bic_ridge_cor = bic_ridge[2],
                 # Elastic Net
                 cv_enet_lam = cv_enet[1],
                 escv_enet_lam = escv_enet[1],
                 aic_enet_lam = aic_enet[1],
                 aicc_enet_lam = aicc_enet[1],
                 bic_enet_lam = bic_enet[1],
                 cv_enet_cor = cv_enet[2],
                 escv_enet_cor = escv_enet[2],
                 aic_enet_cor = aic_enet[2],
                 aicc_enet_cor = aicc_enet[2],
                 bic_enet_cor = bic_enet[2])
    }
  
  # gather data
  model_df <- optimal_lambda %>%
    gather(cor_type, cor, contains("cor")) %>%
    filter(!is.na(cor))
  
  # get lambda corresponding to row
  model_df$lam <- 
    apply(model_df, 1, function(x) x[sub("cor", "lam", x["cor_type"])])
  
  # get max cor by voxel
  model_df_agg <- aggregate(cor ~ voxel, model_df, max)
  # merge back
  model_df_opt <- merge(model_df_agg, 
                        select(model_df, voxel, cor_type, cor, lam))
  # If same cor, choose model with largest lambda (greatest penalty)
  model_df_opt_agg <- aggregate(lam ~ voxel, model_df_opt, max)
  # merge back
  model_df_opt_max <- merge(model_df_opt_agg, model_df_opt)
  
  return(list(optimal_lambda, model_df_opt_max))
}

# Make predictions on test set
predictTest = function(model_df_opt_max) {
  # Arguments:
  #   model_df_opt_max: 2nd output of bestVoxel function from lambda_opt_corr.R
  #
  # Returns:
  #   voxel_corr_test: prediction correlations for each voxel's model on 
  #                    validation and test sets.
  
  # get alphas
  model_df_opt_max$alph <-
    with(model_df_opt_max,
         ifelse(grepl("ridge", cor_type),
                0,
                ifelse(grepl("enet", cor_type),
                       0.5,
                       1)))
  
  model_df_opt_max_test <- distinct(model_df_opt_max, voxel, lam, alph)
  
  test_performance <-
    foreach(i = 1:nrow(model_df_opt_max_test), .combine = rbind) %dopar% {
      
      alph <- filter(model_df_opt_max_test, voxel == i)$alph
      lam <- as.numeric(filter(model_df_opt_max_test, voxel == i)$lam)
      
      voxel_df <-
        select(resp_feat_df,
               names(resp_feat_df)[grepl(paste0("(voxel_resp", i, "$|feat)"),
                                         names(resp_feat_df))])
      
      names(voxel_df)[1] = "voxel_resp"
      x <- model.matrix(voxel_resp~.,voxel_df)[,-1]
      y <- voxel_df$voxel_resp
      
      # Fit model for voxel
      model <- glmnet(x[-test,], y[-test], alpha = alph)
      
      # predict on validation set using optimal lambda
      pred <- predict(model, s = lam, newx = x[test,])
      # compute validation set cor
      cor <- cor(pred, y[test])
      
      data.frame(voxel = i,
                 test_cor = cor)
    }
  
  voxel_corr_test <- distinct(model_df_opt_max, voxel, cor) %>%
    inner_join(test_performance, by = "voxel") %>%
    arrange(voxel) %>%
    mutate(diff = cor - test_cor)
  
  return(voxel_corr_test)
}

# Make predictions on val_feat
predictValFeat = function(model_df_opt_max) {
  # Arguments:
  #   model_df_opt_max: 2nd output of bestVoxel function from lambda_opt_corr.R
  #
  # Returns:
  #   pred_val_feat: response predictions 
  
  # clean val_feat:
  val_feat <- data.frame(val_feat)
  names(val_feat) <- paste("feat", 1:ncol(val_feat), sep = "")
  x_val <- model.matrix(~.,val_feat)[,-1]
  
  # get alphas
  model_df_opt_max$alph <-
    with(model_df_opt_max,
         ifelse(grepl("ridge", cor_type),
                0,
                ifelse(grepl("enet", cor_type),
                       0.5,
                       1)))
  
  # get lambda and alpha for voxel 1 model
  model_df_opt_max_val <- distinct(model_df_opt_max, voxel, lam, alph) %>%
    filter(voxel == 1)

      alph <- model_df_opt_max_val$alph
      lam <- as.numeric(model_df_opt_max_val$lam)
      
      voxel_df <-
        select(resp_feat_df,
               names(resp_feat_df)[grepl(paste0("(voxel_resp", 1, "$|feat)"),
                                         names(resp_feat_df))])
      
      names(voxel_df)[1] = "voxel_resp"
      x <- model.matrix(voxel_resp~.,voxel_df)[,-1]
      y <- voxel_df$voxel_resp
      
      # Fit model for voxel
      model <- glmnet(x, y, alpha = alph)
      
      # predict on validation set using optimal lambda
      pred <- predict(model, s = lam, newx = x_val)

      return(pred)
}
