# Plot functions:

###########################
# Response Correlation Plot
###########################
corrPlot <- function(resp_dat) {
  # Arguments:
  #   resp_dat: responses from the 20 voxels to 1750 images
  #
  # Returns:
  #   correlation plot of voxel responses
  corrplot(cor(resp_dat), tl.col="black",
           title = "Figure 1: correlation of responses for the 20 voxels",
           mar=c(0,0,1.5,0),
           diag = FALSE)
  
}

###########################
# Model Correlation Plot
###########################
modelCorrPlot <- function(optimal_lambda_df) {
  # Arguments:
  #   optimal_lambda_df: 1st output of bestVoxel function from 
  #                      lambda_opt_corr.R
  #
  # Returns:
  #   correlation plot of model performances
  
# get only correlations
model_corr <- optimal_lambda_df %>%
  select(names(optimal_lambda_df)[grepl("cor|voxel", names(optimal_lambda_df))])

# clean names
names(model_corr) <- gsub("_", " ", gsub("_cor", "", names(model_corr)))
# make matrix
model_corr_matrix <- as.matrix(select(model_corr, -voxel))

corrplot(model_corr_matrix, 
         method = "color", is.corr = FALSE,
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         tl.cex = 0.8, na.label = "o", cl.ratio = 0.1, cl.align = "l",
         mar=c(0,0,1.5,0),
         title = "Figure 2: correlation by model, criteria, and voxel")
}

corBar <- function(voxel_corr_test, vox) {
  # Arguments:
  #   voxel: voxel corresponding to model to plot
  #   voxel_corr_test: prediction correlations for each voxel's model on 
  #                    validation and test sets.
  #                    Output of predictTest function from lambda_opt_corr.R
  #
  # Returns:
  #   a bar chart comparing validation vs test correlation metric
  
  voxel_plot <- filter(voxel_corr_test, voxel == vox) %>%
    rename(validation = cor, test = test_cor)
  voxel_plot_gather <- gather(voxel_plot, cor_type, cor, validation, test)
  
  if(vox == 7) {
    return(
  ggplot(voxel_plot_gather, aes(x=cor_type, y = cor)) + 
    geom_bar(stat = "identity", fill = "darkblue") +
    geom_text(aes(label = round(cor, 2)), vjust=-0.5, color="black", 
              size= 3) +
    ggtitle("voxel 7 elastic net") +
    ylab("correlation") +
      ylim(0,0.6) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10))
    )
  }
  if(vox == 20) {
    return(
    ggplot(voxel_plot_gather, aes(x=cor_type, y = cor)) + 
      geom_bar(stat = "identity", fill = "darkblue") +
      geom_text(aes(label = round(cor, 2)), vjust=-0.7, color="black", 
                size= 3) +
      ggtitle("voxel 20 ridge regression") +
      ylab("correlation") +
      ylim(-0.03,0.6) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 10))
    )
  }
    
}

# Box plots by voxel:
voxelBootBox <- function(vox, value, title, ymin, ymax) {
  # Arguments:
  #   vox: voxel number for model
  #   value: lambda or cor
  #   title: include title? yes or no
  #   ymin: min y-axis value
  #   ymax: max y-axis value
  # Returns:
  #   Box plot of value by voxel across bootstrap samples
  
  # select voxel and value
  boot_plot <- boot_df %>%
    select(paste0(value, vox)) %>%
    na.omit
  
  if(value == "cor") {value = "correlation"}
  
  # gather
  boot_plot_gather <- gather(boot_plot)
  boot_plot_gather$key = value
  names(boot_plot_gather) <- c(value, "value")
  
  # for aes string
  y_val <- "value"
  
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  if(title == "yes") {
  return(
  # make box plot
  ggplot(boot_plot_gather, aes_string(x = value, y = y_val)) + 
    geom_boxplot() +
    ggtitle(paste0("voxel ", vox)) + 
    theme_bw() +
    theme(legend.position="none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10)) +
    scale_y_continuous(labels=scaleFUN, limits = c(ymin, ymax)))
    } else if (title == "no") {
      return(
        # make box plot
        ggplot(boot_plot_gather, aes_string(x = value, y = y_val)) + 
          geom_boxplot() +
          theme_bw() +
          theme(legend.position="none",
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                plot.title = element_text(hjust = 0.5)) +
         scale_y_continuous(labels=scaleFUN, limits = c(ymin, ymax)))
    }
}

# plot correlation for each voxel:
voxelCorr <- function(model_df_opt_max) {
  # Arguments:
  #   model_df_opt_max: 2nd output of bestVoxel function from lambda_opt_corr.R
  #
  # Returns:
  #   Bar chart of correlations for each voxel
  
  # alpha
  model_df_opt_max$alph <-
    with(model_df_opt_max,
         ifelse(grepl("ridge", cor_type),
                0,
                ifelse(grepl("enet", cor_type),
                       0.5,
                       1)))
  # label by model
  model_df_opt_max$model <-
    with(model_df_opt_max,
         ifelse(alph == 0,
                "ridge",
                ifelse(alph == 0.5,
                       "elastic net",
                       "lasso")))
  
  # get plotting df
  model_corr_plot <- model_df_opt_max %>%
    distinct(voxel, lam, cor, model) %>%
    rename(correlation = cor) %>%
    arrange(voxel)
  
  # make voxel factor
  model_corr_plot$voxel <- as.factor(model_corr_plot$voxel)
  
  # plot
  ggplot(model_corr_plot, aes(x = voxel, y= correlation, fill = model)) +
   geom_bar(stat="identity") +
   ggtitle("Figure 3: correlation between fitted and observed values by voxel")+
   geom_text(aes(label= round(correlation,2)), vjust=-0.3, size = 2) +
   theme_minimal() +
   theme(plot.title = element_text(size = 11))
}

# plot features across boot samples:
bootFeatureLine <- function(voxel7_feat_boot) {
  # Arguments:
  #   voxel7_feat_boot: 4th output of bestVoxel function from lambda_opt_corr.R
  #
  # Returns:
  #   Line chart of number of features across bootstrap sample and orig
  
  # rename
  voxel7_feat_boot <- voxel7_feat_boot %>%
    rename("total features" = n_feat, original = n_feat_orig)
  
  # gather
  boot_gather_feat <- gather(voxel7_feat_boot, feat, `number of features`,
                             `total features`, original)
  boot_gather_feat$`bootstrap sample` <- as.factor(boot_gather_feat$boot)
  
  # plot
  ggplot(boot_gather_feat, aes(x = `bootstrap sample`, y = `number of features`, 
                               group = feat)) +
    ggtitle(paste0("Figure 8: total features and original features by ",
                   "bootstrap\nsample for voxel 7 model")) +
    geom_line(aes(color=feat)) +
    geom_point(aes(color=feat)) +
    ylim(0, 80) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.title = element_blank()) +
    scale_x_discrete(breaks=seq(0, 100, 10))
}