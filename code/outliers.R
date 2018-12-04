# Find outliers for voxel 7 and voxel 20 models:

voxel7Accuracy = function(accuracy) {
  # Arguments:
  #   accuracy: specify most or least accurate predictions
  #
  # Returns:
  #   images with the most or least accurate predictions for voxel 7

# voxel 7:
i = 7
voxel_df <-
  select(resp_feat_df, 
         names(resp_feat_df)[grepl(paste0("(voxel_resp", i, "$|feat)"),
                                   names(resp_feat_df))])
names(voxel_df)[1] = "voxel_resp"
x <- model.matrix(voxel_resp~.,voxel_df)[,-1]
y <- voxel_df$voxel_resp

#############
# Elastic Net
#############

# Fit model for voxel
model <- glmnet(x[train,], y[train], alpha = 0.5)

# ESCV
set.seed(1)
escv_glm <- escv.glmnet(x[train,], y[train], alpha = 0.5, nfolds = 10)
lambda_opt <- escv_glm$lambda.escv
# predict on validation set using optimal lambda
pred <- predict(model, s = lambda_opt, newx = x[validation,])
# compute res
res <- data.frame(y[validation] - pred)
colnames(res) <- "res"
res$image <- row.names(res)

if (accuracy == "least") {
  # Least accurate
  return(select(slice(arrange(res,desc(abs(res))), 1:6), image))
} else if (accuracy == "most") {
  # Most accurate
  return(select(slice(arrange(res,abs(res)), 1:6), image))
}
}

# # Least accurate
# voxel7Accuracy("least")
# 
# # Read in a raw image.
# img1 <- ReadImage(1313)
# image(img1, col = gray((1:500) / 501))
# # no people
# img1 <- ReadImage(680)
# image(img1, col = gray((1:500) / 501))
# # no people
# img1 <- ReadImage(1309)
# image(img1, col = gray((1:500) / 501))
# # people
# img1 <- ReadImage(1294)
# image(img1, col = gray((1:500) / 501))
# # people
# img1 <- ReadImage(622)
# image(img1, col = gray((1:500) / 501))
# # people
# img1 <- ReadImage(547)
# image(img1, col = gray((1:500) / 501))
# # people
# 
# # 4 of top 6 least accurate are images with people
# 
# # Most accurate
# voxel7Accuracy("most")
# 
# # Read in a raw image
# img1 <- ReadImage(1314)
# image(img1, col = gray((1:500) / 501))
# # no people
# img1 <- ReadImage(1085)
# image(img1, col = gray((1:500) / 501))
# # no people
# img1 <- ReadImage(1637)
# image(img1, col = gray((1:500) / 501))
# # no people
# img1 <- ReadImage(941)
# image(img1, col = gray((1:500) / 501))
# # no people
# img1 <- ReadImage(1138)
# image(img1, col = gray((1:500) / 501))
# # no people
# img1 <- ReadImage(101)
# image(img1, col = gray((1:500) / 501))
# # no people
# 
# # 0 of top 6 most accurate are images with people


voxel20Accuracy = function(accuracy) {
  # Arguments:
  #   accuracy: specify most or least accurate predictions
  #
  # Returns:
  #   images with the most or least accurate predictions for voxel 20
  
# voxel 20:

i = 20
voxel_df <-
  select(resp_feat_df, 
         names(resp_feat_df)[grepl(paste0("(voxel_resp", i, "$|feat)"),
                                   names(resp_feat_df))])
names(voxel_df)[1] = "voxel_resp"
x <- model.matrix(voxel_resp~.,voxel_df)[,-1]
y <- voxel_df$voxel_resp

########
# Ridge
########

# Fit model for voxel
model = glmnet(x[train,], y[train], alpha = 0)

# aicc
set.seed(1)
aicc_glm <- ic.glmnet(x[train,], y[train], alpha = 0, crit = "aicc")
lambda_opt <- aicc_glm$lambda
# predict on validation set using optimal lambda
pred <- predict(model, s = lambda_opt, newx = x[validation,])
# compute res
res <- data.frame(y[validation] - pred)
colnames(res) <- "res"
res$image <- row.names(res)

if (accuracy == "least") {
  # Least accurate
  return(select(slice(arrange(res,desc(abs(res))), 1:10), image))
} else if (accuracy == "most") {
  # Most accurate
  return(select(slice(arrange(res,abs(res)), 1:10), image))
}
}

# # Least accurate
# voxel20Accuracy("least")
# 
# # Read in a raw image.
# img1 <- ReadImage(430)
# image(img1, col = gray((1:500) / 501))
# # monkey
# img1 <- ReadImage(1409)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(1046)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(1590)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(133)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(512)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(282)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(124)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(574)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(106)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# 
# # 1 of 10 least accurate are images with monkeys
# 
# # Most accurate
# voxel20Accuracy("most")
# 
# # Read in a raw image
# img1 <- ReadImage(748)
# image(img1, col = gray((1:500) / 501))
# # monkey
# img1 <- ReadImage(1451)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(1013)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(5)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(948)
# image(img1, col = gray((1:500) / 501))
# # monkey
# img1 <- ReadImage(862)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(1082)
# image(img1, col = gray((1:500) / 501))
# # monkey
# img1 <- ReadImage(593)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(584)
# image(img1, col = gray((1:500) / 501))
# # no monkey
# img1 <- ReadImage(234)
# image(img1, col = gray((1:500) / 501))
# # monkey

# 4 of 10 most accurate are images with monkeys