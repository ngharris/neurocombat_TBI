
# libraries for analysis/harmonization
library(dplyr)
library(stringr)
library(neuroCombat)
#library(rstatix)
#library(effectsize)
#library(raster)
#library(caret)
#library(compute.es)

# libraries for visualization
library(tidyr)
#library(tidyverse)
library(ggpubr)
#library(ggsci)
#library(gridExtra)
#library(cluster)
library(readr)
#library(Rtsne)
library(ggplot2)
#library(ggfortify)
#library(patchwork)
#library(dabestr)

#MMRM
library(lmerTest)
library(afex)
library(lme4)
library(multcomp)
#library(EnvStats)


volumes_global = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/SUMMARY_corrected_data_v2/PREharm/common_space_common_shams/all3site_DATA_VOLUME_zmap_thresh_z3_1_lower_upper.txt")
volumes_global_4site = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/Analyze07/GlobalSham_Univariate/DATA_VOLUME_zmap_thresh_z3_1_lower_upper.txt")
AD_global = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/rethresh_original_univ_data/UCLA/DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/rethresh_original_univ_data/UJHU/DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/rethresh_original_univ_data/USUH/DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt")
MD_global = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/rethresh_original_univ_data/UCLA/DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/rethresh_original_univ_data/UJHU/DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/rethresh_original_univ_data/USUH/DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt")
RD_global = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/rethresh_original_univ_data/UCLA/DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/rethresh_original_univ_data/UJHU/DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/rethresh_original_univ_data/USUH/DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt")

volumes_postvoxelharmonization_global = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/SUMMARY_corrected_data_v2/PSTharm/common_space_common_shams/all3site_DATA_VOLUME_zmap_thresh_z3_1_lower_upper.txt")
AD_post_global = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/univariate_harmonized_data_other/UCLA/DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/univariate_harmonized_data_other/UJHU/DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/univariate_harmonized_data_other/USUH/DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt")
MD_post_global = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/univariate_harmonized_data_other/UCLA/DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/univariate_harmonized_data_other/UJHU/DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/univariate_harmonized_data_other/USUH/DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt")
RD_post_global = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/univariate_harmonized_data_other/UCLA/DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/univariate_harmonized_data_other/UJHU/DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt", "/Users/rachelfox/R_Projects/DataHarmonizationFiles/univariate_harmonized_data_other/USUH/DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt")

locations <- function(volumes, sitenum){
  batch = c()
  i <- 1
  if (sitenum == "Threesite"){
    batch = rep(c(1, 2, 3), times=c(99, 89, 96)) # 3 sites
  }
  else {
    batch = rep(c(1, 2, 3, 4), times=c(99, 89, 96, 59)) # 4 sites
  }
  return(batch)
}

# filter out only high/low FA volume and covariates, combine into one dataframe
allvolumes <- function(volumes, commonspace){
  allv <- data.frame()
  for (file in volumes){
    a <- read.table(file)
    b <- dplyr::select(a, 1, 3, 5)
    allv <- rbind(allv, b)
  }
  
  #adds columns with identifiers (treatment, number, day, sex) and removes original ID column
  if (commonspace=="Commonspace" | commonspace=="Threesite") {
    #if (commonspace=="Threesite") {
    allv$V1 <- substr(allv$V1, 6, nchar(allv$V1))
    #allv[c('treatment','number','day','sex')] <- str_split_fixed(allv$V1, '_', 5)
    
  }
  #if (commonspace=="Commonspace") {
  #allv$V1 <- substr(allv$V1, 6, nchar(allv$V1))
  #allv[c('treatment','number','day','sex')] <- str_split_fixed(allv$V1, '_', 5)
  
  # }
  
  
  allv[c('treatment','number','day','sex')] <- str_split_fixed(allv$V1, '_', 4)
  allv <- subset(allv, select = -c(V1))
  
  print(allv)
  colnames(allv)[1:2] <- c("low_FA_volume", "high_FA_volume")
  
  # converts to cubic mm
  allv$low_FA_volume <- (as.numeric(allv$low_FA_volume) / 1000)
  allv$high_FA_volume <- (as.numeric(allv$high_FA_volume) / 1000)
  
  print(commonspace)
  
  if (commonspace =="Commonspace") {
    allv$location <- rep(c(1, 2, 3, 4), times=c(99, 89, 96, 59)) # four site
  }
  else if (commonspace=="Threesite") {
    allv <- allv[1:284, ]
    allv$location <- rep(c(1, 2, 3), times=c(99, 89, 96))
  }
  else if (commonspace=="Newmeas_voxel") {
    allv <- allv[1:281, ]
    allv$location <- rep(c(1, 2, 3), times=c(99, 86, 96))
  }
  else if (commonspace=="Threesite_ss") {
    allv$location <- rep(c(1, 2, 3), times=c(99, 89, 96))
  }
  else {
    allv$location <- rep(c(1, 2, 3, 4), times=c(99, 89, 96, 59))
  }
  print(allv)
  return(allv)
}

harmonize_split_together <- function(allvs, sitenum){
  big <- list()
  if (sitenum=="Threesite") {
    #pooled data harmonization - first entry in list of data
    mod <- model.matrix(~allvs$treatment+allvs$day+allvs$sex+allvs$atrophy)
    
    reduced_allvs <- t(dplyr::select(allvs, 1, 2))
    print(reduced_allvs)
    harm_v <- neuroCombat(dat = reduced_allvs, batch = locations(allvs, "Threesite"), mod = mod, mean.only = TRUE)
    big_v <- cbind(t(harm_v$dat.combat), dplyr::select(allvs, 3, 4, 5, 6, 7,8))
    
    #pooled original data - second entry in list of data 
    orig_big_v <- cbind(t(harm_v$dat.original), dplyr::select(allvs, 3, 4, 5, 6, 7,8))
    
    big <- append(big, list(big_v))
    big <- append(big, list(orig_big_v))
    
    return(big)
  }
  
  else { # four site
    #pooled data harmonization - first entry in list of data
    mod <- model.matrix(~allvs$treatment+allvs$day+allvs$sex)
    
    reduced_allvs <- t(dplyr::select(allvs, 1, 2))
    print(reduced_allvs)
    harm_v <- neuroCombat(dat = reduced_allvs, batch = locations(allvs, "Foursite"), mod = mod, mean.only = TRUE)
    big_v <- cbind(t(harm_v$dat.combat), dplyr::select(allvs, 3, 4, 5, 6, 7)) # four site
    
    #pooled original data - second entry in list of data
    orig_big_v <- cbind(t(harm_v$dat.original), dplyr::select(allvs, 3, 4, 5, 6, 7)) # four site
    
    big <- append(big, list(big_v))
    big <- append(big, list(orig_big_v))
    
    return(big)
  }

}

# get effect sizes
get_e_size <- function(sepoutput){
  d_values <- matrix(ncol = 2)
  colnames(d_values) <- c("low_FA_volume", "high_FA_volume")
  
  for(df in sepoutput){
    for (site in 1:3) {
      # effect size equation
      dlow_3 <- (mean(df[,1][which(df$treatment == "CCI" & df$day == "d03" & df$location == site)]) - mean(df[,1][which(df$treatment == "SHM" & df$day == "d03" & df$location == site)]))/sd(df[,1][which(df$day == "d03" & df$location == site)])
      dlow_30 <- (mean(df[,1][which(df$treatment == "CCI" & df$day == "d30" & df$location == site)]) - mean(df[,1][which(df$treatment == "SHM" & df$day == "d30" & df$location == site)]))/sd(df[,1][which(df$day == "d30" & df$location == site)])
      dhigh_3 <- (mean(df[,2][which(df$treatment == "CCI" & df$day == "d03" & df$location == site)]) - mean(df[,2][which(df$treatment == "SHM" & df$day == "d03" & df$location == site)]))/sd(df[,2][which(df$day == "d03" & df$location == site)])
      dhigh_30 <- (mean(df[,2][which(df$treatment == "CCI" & df$day == "d30" & df$location == site)]) - mean(df[,2][which(df$treatment == "SHM" & df$day == "d30" & df$location == site)]))/sd(df[,2][which(df$day == "d30" & df$location == site)])
      d_values <- rbind(d_values, c(dlow_3, dhigh_3))
      d_values <- rbind(d_values, c(dlow_30, dhigh_30))
    }
  }
  d_values = d_values[-1,]
  rownames(d_values) <- c("orig_03_1", "orig_30_1", "harm_03_1", "harm_30_1", "orig_03_2", "orig_30_2", "harm_03_2", "harm_30_2", "orig_03_3", "orig_30_3", "harm_03_3", "harm_30_3")
  return(d_values)
}

# construct dataset and perform harmonization
allvs <- allvolumes(volumes_global, "Threesite") # for four-site, AD/MD/RD data use "Commonspace", and "Newmeas_voxel" respectively
allvs_4site <- allvolumes(volumes_global_4site, "Commonspace") # 4 site
addmeas <- allvolumes(MD_global, "Threesite_ss")

atrophy = c("/Users/rachelfox/R_Projects/DataHarmonizationFiles/Analyze07/atrophy-expansionData_uh3/ATROPHY_all.txt")
a <- read.table(atrophy)
a[c('treatment','number','day','sex', 'm')] <- str_split_fixed(a$V1, '_', 5)
a <- subset(a, select = -c(V1, V2, m))
colnames(a)[1] <- c("atrophy")
a$atrophy = a$atrophy / 1000

merged_df <- merge(allvs, a, by = c("treatment", "number", "day", "sex"))
merged_df <- merge(addmeas, a, by = c("treatment", "number", "day", "sex")) # replace to use for AD/MD/RD
atrophy <- merged_df[, "atrophy"]
merged_df <- merged_df[, -which(names(merged_df) == "atrophy")]
merged_df <- cbind(merged_df, atrophy)
merged_df <- merged_df[,c(5,6,1,2,3,4,7,8)]

togoutput <- harmonize_split_together(merged_df, "Threesite")
togoutput_4site <- harmonize_split_together(allvs_4site, "Foursite")
valscombined <- togoutput[[1]]
origvalscombined <- togoutput[[2]]
valscombined_4site <- togoutput_4site[[1]]
origvalscombined_4site <- togoutput_4site[[2]]
V_effect_size_together <- get_e_size(togoutput[1:2])

allvs_postharmonized <- allvolumes(volumes_postvoxelharmonization_global, "Threesite")
allvs_postharmonized <- allvolumes(AD_post_global, "Newmeas_voxel") # replace with MD/RD/AD
merged_post <- merge(allvs_postharmonized, a, by = c("treatment", "number", "day", "sex"))
atrophy <- merged_post[, "atrophy"]
merged_post <- merged_post[, -which(names(merged_post) == "atrophy")]
merged_post <- cbind(merged_post, atrophy)
allvs_postharmonized <- merged_post[,c(5,6,1,2,3,4,7,8)]

# plots to produce figure 2
shamsharmonized_4site <- valscombined_4site[valscombined_4site$treatment=="SHM", ]
shamsoriginal_4site <- origvalscombined_4site[origvalscombined_4site$treatment == "SHM", ]
model <- aov(high_FA_volume ~ location, data = shamsharmonized_4site)
df_plot <- data.frame(
  FittedValues = model$fitted.values,
  Residuals = model$residuals,
  Location = as.factor(shamsharmonized_4site$location)
)
ggplot(df_plot, aes(x = Location, y = Residuals, color = Location)) + geom_point(size=3, alpha=0.8) + 
  geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  labs(x = expression("Site"), y = "Residuals", title = expression("Residuals versus Fitted Values \nby Site," ~ FA[low]), color="Site") + 
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") + 
  theme_pubr(base_size = 14) + theme(legend.position = "bottom", legend.direction = "horizontal") + 
  theme(legend.position = "bottom", legend.direction = "horizontal", panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("#1F78B4", "#33A02C", "#6A3D9A", "#E31A1C"))


# linear mixed model ANOVA to get effect of covariates - used in figure 3
get_mmrm <- function(measure, df) {
  lmeModel = lm(measure ~ treatment*sex*location*atrophy, data=df)
  mmrm <- anova(lmeModel)
  return(mmrm) 
}
# original
low_volume_original_d03 <- get_mmrm(origvalscombined[origvalscombined$day == "d03", ]$low_FA_volume, origvalscombined[origvalscombined$day == "d03", ])
low_volume_original_d30 <- get_mmrm(origvalscombined[origvalscombined$day == "d30", ]$low_FA_volume, origvalscombined[origvalscombined$day == "d30", ])
high_volume_original_d03 <- get_mmrm(origvalscombined[origvalscombined$day == "d03", ]$high_FA_volume, origvalscombined[origvalscombined$day == "d03", ])
high_volume_original_d30 <- get_mmrm(origvalscombined[origvalscombined$day == "d30", ]$high_FA_volume, origvalscombined[origvalscombined$day == "d30", ])

# univariate harmonization
low_volume_harmonized_d03 <- get_mmrm(valscombined[valscombined$day == "d03", ]$low_FA_volume, valscombined[valscombined$day == "d03", ])
low_volume_harmonized_d30 <- get_mmrm(valscombined[valscombined$day == "d30", ]$low_FA_volume, valscombined[valscombined$day == "d30", ])
high_volume_harmonized_d03 <- get_mmrm(valscombined[valscombined$day == "d03", ]$high_FA_volume, valscombined[valscombined$day == "d03", ])
high_volume_harmonized_d30 <- get_mmrm(valscombined[valscombined$day == "d30", ]$high_FA_volume, valscombined[valscombined$day == "d30", ])

# voxel-level harmonization
low_volume_postharmonized_d03 <- get_mmrm(allvs_postharmonized[allvs_postharmonized$day == "d03", ]$low_FA_volume, allvs_postharmonized[allvs_postharmonized$day == "d03", ])
low_volume_postharmonized_d30 <- get_mmrm(allvs_postharmonized[allvs_postharmonized$day == "d30", ]$low_FA_volume, allvs_postharmonized[allvs_postharmonized$day == "d30", ])
high_volume_postharmonized_d03 <- get_mmrm(allvs_postharmonized[allvs_postharmonized$day == "d03", ]$high_FA_volume, allvs_postharmonized[allvs_postharmonized$day == "d03", ])
high_volume_postharmonized_d30 <- get_mmrm(allvs_postharmonized[allvs_postharmonized$day == "d30", ]$high_FA_volume, allvs_postharmonized[allvs_postharmonized$day == "d30", ])


# plots to produce figure 3
modify_df <- function(bxp_df){
  bxp_df$treatLoc[bxp_df$treatment=="SHM" & bxp_df$location=="1"] <- "1_SHM.Site1"
  bxp_df$treatLoc[bxp_df$treatment=="SHM" & bxp_df$location=="2"] <- "1_SHM.Site2"
  bxp_df$treatLoc[bxp_df$treatment=="SHM" & bxp_df$location=="3"] <- "1_SHM.Site3"
  bxp_df$treatLoc[bxp_df$treatment=="SHM" & bxp_df$location=="4"] <- "1_SHM.Site4"
  bxp_df$treatLoc[bxp_df$treatment=="CCI" & bxp_df$location=="1"] <- "2_CCI.Site1"
  bxp_df$treatLoc[bxp_df$treatment=="CCI" & bxp_df$location=="2"] <- "2_CCI.Site2"
  bxp_df$treatLoc[bxp_df$treatment=="CCI" & bxp_df$location=="3"] <- "2_CCI.Site3"
  bxp_df$treatLoc[bxp_df$treatment=="CCI" & bxp_df$location=="4"] <- "2_CCI.Site4"
  return(bxp_df)
}

boxplot_df_orig <- modify_df(origvalscombined)
boxplot_df_03_orig <- boxplot_df_orig[boxplot_df_orig$day == "d03", ]
boxplot_df_30_orig <- boxplot_df_orig[boxplot_df_orig$day == "d30", ]

boxplot_df_har <- modify_df(valscombined)
boxplot_df_03_har <- boxplot_df_har[boxplot_df_har$day == "d03", ]
boxplot_df_30_har <- boxplot_df_har[boxplot_df_har$day == "d30", ]

boxplot_df_posth <- modify_df(allvs_postharmonized)
boxplot_df_03_posth <- boxplot_df_posth[boxplot_df_posth$day == "d03", ]
boxplot_df_30_posth <- boxplot_df_posth[boxplot_df_posth$day == "d30", ]

make_violin_plots <- function(dataf, yval, ylab, title, fillcol, ylm) {
  violinp <- ggplot(dataf, aes(x=treatment, y=yval, fill=treatLoc)) + ylim(0, ylm) + 
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = fillcol, alpha = 0.1) +
    geom_violin(position=position_dodge(0.8), width=1.1, alpha=0.75) + 
    geom_boxplot(position=position_dodge(0.8), width=0.2) + 
    stat_summary(fun=median, show.legend = FALSE, geom="crossbar", position=position_dodge(0.8), width=0.5) +
    labs(x="", y=ylab, title=title) +
    theme_pubr(base_size = 11) + scale_x_discrete(limits=c("SHM", "CCI"), labels=c("Sham", "CCI")) +
    #scale_fill_manual(values=c("#A6CEE3", "#B2DF8A","#CAB2D6", "#FB9A99", "#1F78B4", "#33A02C", "#6A3D9A", "#E31A1C"), name="Injury and Location", labels=c("Sham Site 1", "Sham Site 2", "Sham Site 3", "Sham Site 4", "CCI Site 1", "CCI Site 2", "CCI Site 3", "CCI Site 4")) # uncomment for 4 sites
    scale_fill_manual(values=c("#A6CEE3", "#B2DF8A","#CAB2D6", "#1F78B4", "#33A02C", "#6A3D9A"), name="Injury and Location", labels=c("Sham Site 1", "Sham Site 2", "Sham Site 3", "CCI Site 1", "CCI Site 2", "CCI Site 3"))
  return (violinp)
}

orig_d3_low <- make_violin_plots(boxplot_df_03_orig, boxplot_df_03_orig$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 03 \nBefore Harmonization", "#FFFFFF", 10)
orig_d3_high <- make_violin_plots(boxplot_df_03_orig, boxplot_df_03_orig$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 03 \nBefore Harmonization", "#FFFFFF", 50)
orig_d30_low <- make_violin_plots(boxplot_df_30_orig, boxplot_df_30_orig$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 30 \nBefore Harmonization", "#FFFFFF", 10)
orig_d30_high <- make_violin_plots(boxplot_df_30_orig, boxplot_df_30_orig$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 30 \nBefore Harmonization", "#FFFFFF", 50)

har_d3_low <- make_violin_plots(boxplot_df_03_har, boxplot_df_03_har$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 03 \nUnivariate Harmonization", "#f2f2f2", 10)
har_d3_high <- make_violin_plots(boxplot_df_03_har, boxplot_df_03_har$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 03 \nUnivariate Harmonization", "#f2f2f2", 50)
har_d30_low <- make_violin_plots(boxplot_df_30_har, boxplot_df_30_har$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 30 \nUnivariate Harmonization", "#f2f2f2", 10)
har_d30_high <- make_violin_plots(boxplot_df_30_har, boxplot_df_30_har$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 30 \nUnivariate Harmonization", "#f2f2f2", 50)

posth_d3_low <- make_violin_plots(boxplot_df_03_posth, boxplot_df_03_posth$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 03 \nVoxel-Level Harmonization", "#d9d9d9", 10)
posth_d3_high <- make_violin_plots(boxplot_df_03_posth, boxplot_df_03_posth$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 03 \nVoxel-Level Harmonization", "#d9d9d9", 50)
posth_d30_low <- make_violin_plots(boxplot_df_30_posth, boxplot_df_30_posth$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 30 \nVoxel-Level Harmonization", "#d9d9d9", 10)
posth_d30_high <- make_violin_plots(boxplot_df_30_posth, boxplot_df_30_posth$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 30 \nVoxel-Level Harmonization", "#d9d9d9", 50)

