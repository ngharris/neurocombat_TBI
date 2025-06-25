
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


volumes_global = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\all3site_DATA_VOLUME_zmap_thresh_z3_1_lower_upper.txt")
volumes_global_4site = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\4siteDATA_VOLUME_zmap_thresh_z3_1_lower_upper.txt")
AD_global = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD\\UCLA\\DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD\\UJHU\\DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD\\USUH\\DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt")
MD_global = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD\\UCLA\\DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD\\UJHU\\DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD\\USUH\\DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt")
RD_global = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD\\UCLA\\DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD\\UJHU\\DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD\\USUH\\DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt")

volumes_postvoxelharmonization_global = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\voxelwise_harm_all3site_DATA_VOLUME_zmap_thresh_z3_1_lower_upper.txt")
AD_post_global = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD_post\\UCLA\\DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD_post\\UJHU\\DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD_post\\USUH\\DATA_VOLUME_AD_zmap_thresh_z3_1_lower_upper.txt")
MD_post_global = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD_post\\UCLA\\DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD_post\\UJHU\\DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD_post\\USUH\\DATA_VOLUME_MD_zmap_thresh_z3_1_lower_upper.txt")
RD_post_global = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD_post\\UCLA\\DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD_post\\UJHU\\DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt", "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\AD_MD_RD_post\\USUH\\DATA_VOLUME_RD_zmap_thresh_z3_1_lower_upper.txt")

locations <- function(volumes, sitenum){
  batch = c()
  i <- 1
  if (sitenum == "Threesite"){
    batch = rep(c(1, 2, 3), times=c(99, 86, 96)) # 3 sites
  }
  else {
    batch = rep(c(1, 2, 3, 4), times=c(99, 86, 96, 59)) # 4 sites
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
    allv <- allv[1:280, ]
    allv$location <- rep(c(1, 2, 3), times=c(99, 85, 96))
  }
  else if (commonspace=="Newmeas_voxel") {
    allv <- allv[1:281, ]
    allv$location <- rep(c(1, 2, 3), times=c(99, 86, 96))
  }
  else if (commonspace=="Threesite_ss") {
    allv$location <- rep(c(1, 2, 3), times=c(99, 89, 96))
  }
  else {
    allv$location <- rep(c(1, 2, 3, 4), times=c(99, 86, 96, 59))
  }
  print(allv)
  return(allv)
}

harmonize_split_together <- function(allvs1, sitenum){
  big <- list()
  if (sitenum=="Threesite") {
    #pooled data harmonization - first entry in list of data
    #mod <- model.matrix(~allvs1$treatment+allvs1$day+allvs1$sex+allvs1$atrophy)
    mod <- model.matrix(~allvs1$day+allvs1$sex+allvs1$atrophy)
    reduced_allvs <- t(dplyr::select(allvs1, 1, 2))
    print(table(allvs1$location))
    print(length(reduced_allvs))
    print(dim(mod))
    #harm_v <- neuroCombat(dat = reduced_allvs, batch = allvs1$location, mod = mod, mean.only = F)
    harm_v <- sva::ComBat_seq(counts=reduced_allvs,batch=allvs1$location,covar_mod=mod,group = allvs1$treatment)
    #big_v <- cbind(t(harm_v$dat.combat), dplyr::select(allvs1, 3, 4, 5, 6, 7))
    print(dim(harm_v))
    big_v <- cbind(t(harm_v), dplyr::select(allvs1, 3, 4, 5, 6, 7))
    print('bigv ok')
    #pooled original data - second entry in list of data 
    #orig_big_v <- cbind(t(harm_v$dat.original), dplyr::select(allvs1, 3, 4, 5, 6, 7))
    orig_big_v <- cbind(t(reduced_allvs), dplyr::select(allvs1, 3, 4, 5, 6, 7))
    big <- append(big, list(big_v))
    big <- append(big, list(orig_big_v))
    
    return(big)
    print("done")
  }
  
  else { # four site
    #pooled data harmonization - first entry in list of data
    mod <- model.matrix(~allvs1$treatment+allvs1$day+allvs1$sex)
    
    reduced_allvs <- t(dplyr::select(allvs1, 1, 2))
    print(reduced_allvs)
    harm_v <- neuroCombat(dat = reduced_allvs, batch = allvs1$location, mod = mod, mean.only = TRUE)
    big_v <- cbind(t(harm_v$dat.combat), dplyr::select(allvs1, 3, 4, 5, 6, 7)) # four site
    
    #pooled original data - second entry in list of data
    orig_big_v <- cbind(t(harm_v$dat.original), dplyr::select(allvs1, 3, 4, 5, 6, 7)) # four site
    
    big <- append(big, list(big_v))
    big <- append(big, list(orig_big_v))
    
    return(big)
    print("done")
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
addmeas <- allvolumes(AD_global, "Threesite_ss")

atrophy = c("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\ATROPHY_alltxt.txt")
a <- read.table(atrophy)
a[c('treatment','number','day','sex', 'm')] <- str_split_fixed(a$V1, '_', 5)
a <- subset(a, select = -c(V1, V2, m))
colnames(a)[1] <- c("atrophy")
a$atrophy = a$atrophy / 1000
a$number <- str_remove(a$number,"^0")
#merged_df <- merge(allvs, a, by = c("treatment", "number", "day", "sex")) # comment/uncomment this line and the next to switch between FA & AD/MD/RD
merged_df <- merge(addmeas, a, by = c("treatment", "number", "day", "sex")) # replace to use for AD/MD/RD
atrophy <- merged_df[, "atrophy"]
merged_df <- merged_df[, -which(names(merged_df) == "atrophy")]
merged_df <- cbind(merged_df, atrophy)
merged_df <- merged_df[,c(5,6,1,2,3,4,7,8)]

#merged_df$low_FA_volume[merged_df$low_FA_volume == 0] <- 1e-8
#merged_df$low_FA_volume <- log(subset_d$low_FA_volume)

#merged_df$high_FA_volume[merged_df$high_FA_volume == 0] <- 1e-8
#merged_df$high_FA_volume <- log(subset_d$high_FA_volume)

togoutput <- harmonize_split_together(merged_df, "Threesite")
togoutput_4site <- harmonize_split_together(allvs_4site, "Foursite")
valscombined <- togoutput[[1]]
origvalscombined <- togoutput[[2]]
valscombined_4site <- togoutput_4site[[1]]
origvalscombined_4site <- togoutput_4site[[2]]
V_effect_size_together <- get_e_size(togoutput[1:2])

combat_seq_harm <- harmonize_split_together(merged_df, 'Threesite')
#allvs_postharmonized <- allvolumes(volumes_postvoxelharmonization_global, "Threesite") #comment/uncomment between this line and the next to switch between FA & AD/MD/RD
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

resplot <- ggviolin(df_plot, x = "Location", y = "Residuals", fill = "Location",
         palette =c("#1F78B4", "#33A02C", "#6A3D9A", "#E31A1C"),
         xlab = "Site", ylab = "Residuals", alpha=0.85, width=1, linetype = 'blank', 
         add = c('point','jitter'), add.params = list(color="Location", jitter=0.1, alpha=0.25)) +
  #geom_point(size=3, alpha=0.8) + geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+theme(text = element_text(size = 20), aspect.ratio = 5/3) +
  scale_x_discrete(expand = c(0, 0))
ggsave("Fig2B.tiff", plot = resplot, width=8, units='in')

# linear mixed model ANOVA to get effect of covariates - used in figure 3
get_mmrm_anova <- function(measure, df) {
  lmeModel = lm(measure ~ treatment*sex*location*atrophy, data=df)
  mmrm <- anova(lmeModel)
  return(mmrm) 
}

library(glmmTMB)
library(car)
get_mmrm <- function(measure, df) {
  df <- as.data.frame(df)
  df$low_FA_volume[df$low_FA_volume <= 0] <- 1e-8
  df$low_FA_volume <- round(df$low_FA_volume) #need to do this for poisson
  df$high_FA_volume[df$high_FA_volume <= 0] <- 1e-8
  df$high_FA_volume <- round(df$high_FA_volume)
  df <- as.data.frame(df)
  print(head(df))
  formula_str <- paste0(measure, " ~ treatment * sex * location")
  #formula_str <- paste0(measure, " ~ treatment * sex * location + (1 | number)")
  glmModel <- glmmTMB(
    formula = as.formula(formula_str),
    data = df,
    family = poisson(link = "log")
  )
  mmrm <- Anova(glmModel, test = "Chisq")
  #mmrm <- summary(glmModel)
  return(mmrm)
  #return(summary(glmModel))
}

for (v in list(origvalscombined, valscombined, allvs_postharmonized)) {
  print("New dataframe")
  for (d in c("d03", "d30")) {
    subset_d <- v[v$day == d, ]
    
    print("New day")
    print(d)
    print("Low volume")
    
    subset_d$low_FA_volume[subset_d$low_FA_volume == 0] <- 1e-4
    subset_d <- subset_d[subset_d$low_FA_volume > 0, ]
    #print(get_mmrm("low_FA_volume", subset_d))
    subset_d$low_FA_volume_log <- log(subset_d$low_FA_volume)
    print(get_mmrm_anova(subset_d$low_FA_volume_log, subset_d))
    
    print("High volume")
    subset_d$high_FA_volume[subset_d$high_FA_volume == 0] <- 1e-4
    subset_d <- subset_d[subset_d$high_FA_volume > 0, ]
    #print(get_mmrm("high_FA_volume", subset_d))
    subset_d$high_FA_volume_log <- log(subset_d$high_FA_volume)
    print(get_mmrm_anova(subset_d$high_FA_volume_log, subset_d))
  }
}

origvalscombined$low_FA_volume[origvalscombined$low_FA_volume <= 0] <- 1e-8
low_volume_original_d03 <- get_mmrm(measure = "low_FA_volume", df = origvalscombined[origvalscombined$day == "d03", ])

valscombined$low_FA_volume[valscombined$low_FA_volume <= 0] <- 1e-8
low_volume_harmonized_d03 <- get_mmrm(measure = "low_FA_volume", df = valscombined[valscombined$day == "d03", ])



# original
low_volume_original_d03 <- get_mmrm("low_FA_volume", origvalscombined[origvalscombined$day == "d03", ])
low_volume_original_d30 <- get_mmrm("low_FA_volume", origvalscombined[origvalscombined$day == "d30", ])
high_volume_original_d03 <- get_mmrm("high_FA_volume", origvalscombined[origvalscombined$day == "d03", ])
high_volume_original_d30 <- get_mmrm("high_FA_volume", origvalscombined[origvalscombined$day == "d30", ])

# univariate harmonization
low_volume_harmonized_d03 <- get_mmrm("low_FA_volume", valscombined[valscombined$day == "d03", ])
low_volume_harmonized_d30 <- get_mmrm("low_FA_volume", valscombined[valscombined$day == "d30", ])
high_volume_harmonized_d03 <- get_mmrm("high_FA_volume", valscombined[valscombined$day == "d03", ])
high_volume_harmonized_d30 <- get_mmrm("high_FA_volume", valscombined[valscombined$day == "d30", ])

# voxel-level harmonization
low_volume_postharmonized_d03 <- get_mmrm("low_FA_volume", allvs_postharmonized[allvs_postharmonized$day == "d03", ])
low_volume_postharmonized_d30 <- get_mmrm("low_FA_volume", allvs_postharmonized[allvs_postharmonized$day == "d30", ])
high_volume_postharmonized_d03 <- get_mmrm("high_FA_volume", allvs_postharmonized[allvs_postharmonized$day == "d03", ])
high_volume_postharmonized_d30 <- get_mmrm("high_FA_volume", allvs_postharmonized[allvs_postharmonized$day == "d30", ])


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
  violinp <- ggplot(dataf, aes(x=treatment, y=yval, fill=treatLoc)) + 
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = fillcol, alpha = 0.1) +
    geom_violin(position=position_dodge(0.8), width=1, alpha=0.75) + 
    geom_boxplot(position=position_dodge(0.8), width=0.05) + 
    geom_point(position = position_jitterdodge(seed = 1, jitter.width = 0.01),alpha=0.5) + 
    #stat_summary(fun=median, show.legend = FALSE, geom="crossbar", position=position_dodge(0.8), width=0.05) +
    labs(x="", y=ylab, title=title) +
    scale_x_discrete(limits=c("SHM", "CCI"), labels=c("Sham", "CCI")) +
    #scale_fill_manual(values=c("#A6CEE3", "#B2DF8A","#CAB2D6", "#FB9A99", "#1F78B4", "#33A02C", "#6A3D9A", "#E31A1C"), name="Injury and Location", labels=c("Sham Site 1", "Sham Site 2", "Sham Site 3", "Sham Site 4", "CCI Site 1", "CCI Site 2", "CCI Site 3", "CCI Site 4")) # uncomment for 4 sites
    scale_fill_manual(values=c("#A6CEE3", "#B2DF8A","#CAB2D6", "#1F78B4", "#33A02C", "#6A3D9A"), name="Injury and Location", labels=c("Sham Site 1", "Sham Site 2", "Sham Site 3", "CCI Site 1", "CCI Site 2", "CCI Site 3"))+
    ggpubr::theme_pubr(base_size = 11) #+ ylim(0, ylm)
  return (violinp)
}

orig_d3_low <- make_violin_plots(boxplot_df_03_orig, boxplot_df_03_orig$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 03 \nBefore Harmonization", "#FFFFFF", 30)
orig_d3_high <- make_violin_plots(boxplot_df_03_orig, boxplot_df_03_orig$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 03 \nBefore Harmonization", "#FFFFFF", 50)
orig_d30_low <- make_violin_plots(boxplot_df_30_orig, boxplot_df_30_orig$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 30 \nBefore Harmonization", "#FFFFFF", 30)
orig_d30_high <- make_violin_plots(boxplot_df_30_orig, boxplot_df_30_orig$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 30 \nBefore Harmonization", "#FFFFFF", 50)

har_d3_low <- make_violin_plots(boxplot_df_03_har, boxplot_df_03_har$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 03 \nUnivariate Harmonization", "#f2f2f2", 30)
har_d3_high <- make_violin_plots(boxplot_df_03_har, boxplot_df_03_har$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 03 \nUnivariate Harmonization", "#f2f2f2", 50)
har_d30_low <- make_violin_plots(boxplot_df_30_har, boxplot_df_30_har$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 30 \nUnivariate Harmonization", "#f2f2f2", 30)
har_d30_high <- make_violin_plots(boxplot_df_30_har, boxplot_df_30_har$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 30 \nUnivariate Harmonization", "#f2f2f2", 50)

posth_d3_low <- make_violin_plots(boxplot_df_03_posth, boxplot_df_03_posth$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 03 \nVoxel-Level Harmonization", "#d9d9d9", 30)
posth_d3_high <- make_violin_plots(boxplot_df_03_posth, boxplot_df_03_posth$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 03 \nVoxel-Level Harmonization", "#d9d9d9", 50)
posth_d30_low <- make_violin_plots(boxplot_df_30_posth, boxplot_df_30_posth$low_FA_volume, expression("Low Volume [mm]" ^ 3), "Low Volume, Day 30 \nVoxel-Level Harmonization", "#d9d9d9", 30)
posth_d30_high <- make_violin_plots(boxplot_df_30_posth, boxplot_df_30_posth$high_FA_volume, expression("High Volume [mm]" ^ 3), "High Volume, Day 30 \nVoxel-Level Harmonization", "#d9d9d9", 50)

#kruskal wallis test - group effects
low_d3_treatment <- kruskal.test(data=origvalscombined[origvalscombined$day == "d03", ], low_FA_volume~treatment)
low_d30_treatment <- kruskal.test(data=origvalscombined[origvalscombined$day == "d30", ], low_FA_volume~treatment)
high_d3_treatment <- kruskal.test(data=origvalscombined[origvalscombined$day == "d03", ], high_FA_volume~treatment)
high_d30_treatment <- kruskal.test(data=origvalscombined[origvalscombined$day == "d30", ], high_FA_volume~treatment)

uh_low_d3_treatment <- kruskal.test(data=valscombined[valscombined$day == "d03", ], low_FA_volume~treatment)
uh_low_d30_treatment <- kruskal.test(data=valscombined[valscombined$day == "d30", ], low_FA_volume~treatment)
uh_high_d3_treatment <- kruskal.test(data=valscombined[valscombined$day == "d03", ], high_FA_volume~treatment)
uh_high_d30_treatment <- kruskal.test(data=valscombined[valscombined$day == "d30", ], high_FA_volume~treatment)

vh_low_d3_treatment <- kruskal.test(data=allvs_postharmonized[allvs_postharmonized$day == "d03", ], low_FA_volume~treatment)
vh_low_d30_treatment <- kruskal.test(data=allvs_postharmonized[allvs_postharmonized$day == "d30", ], low_FA_volume~treatment)
vh_high_d3_treatment <- kruskal.test(data=allvs_postharmonized[allvs_postharmonized$day == "d03", ], high_FA_volume~treatment)
vh_high_d30_treatment <- kruskal.test(data=allvs_postharmonized[allvs_postharmonized$day == "d30", ], high_FA_volume~treatment)

#kruskal wallist test - location effects
low_d3_location <- kruskal.test(data=origvalscombined[origvalscombined$day == "d03", ], low_FA_volume~location)
low_d30_location <- kruskal.test(data=origvalscombined[origvalscombined$day == "d30", ], low_FA_volume~location)
high_d3_location <- kruskal.test(data=origvalscombined[origvalscombined$day == "d03", ], high_FA_volume~location)
high_d30_location <- kruskal.test(data=origvalscombined[origvalscombined$day == "d30", ], high_FA_volume~location)

uh_low_d3_location <- kruskal.test(data=valscombined[valscombined$day == "d03", ], low_FA_volume~location)
uh_low_d30_location <- kruskal.test(data=valscombined[valscombined$day == "d30", ], low_FA_volume~location)
uh_high_d3_location <- kruskal.test(data=valscombined[valscombined$day == "d03", ], high_FA_volume~location)
uh_high_d30_location <- kruskal.test(data=valscombined[valscombined$day == "d30", ], high_FA_volume~location)

vh_low_d3_location <- kruskal.test(data=allvs_postharmonized[allvs_postharmonized$day == "d03", ], low_FA_volume~location)
vh_low_d30_location <- kruskal.test(data=allvs_postharmonized[allvs_postharmonized$day == "d30", ], low_FA_volume~location)
vh_high_d3_location <- kruskal.test(data=allvs_postharmonized[allvs_postharmonized$day == "d03", ], high_FA_volume~location)
vh_high_d30_location <- kruskal.test(data=allvs_postharmonized[allvs_postharmonized$day == "d30", ], high_FA_volume~location)

#4 site
low_d3_treatment_4 <- kruskal.test(data=origvalscombined_4site[origvalscombined_4site$day == "d03", ], low_FA_volume~treatment)
low_d30_treatment_4 <- kruskal.test(data=origvalscombined_4site[origvalscombined_4site$day == "d30", ], low_FA_volume~treatment)
high_d3_treatment_4 <- kruskal.test(data=origvalscombined_4site[origvalscombined_4site$day == "d03", ], high_FA_volume~treatment)
high_d30_treatment_4 <- kruskal.test(data=origvalscombined_4site[origvalscombined_4site$day == "d30", ], high_FA_volume~treatment)

uh_low_d3_treatment_4 <- kruskal.test(data=valscombined_4site[valscombined_4site$day == "d03", ], low_FA_volume~treatment)
uh_low_d30_treatment_4 <- kruskal.test(data=valscombined_4site[valscombined_4site$day == "d30", ], low_FA_volume~treatment)
uh_high_d3_treatment_4 <- kruskal.test(data=valscombined_4site[valscombined_4site$day == "d03", ], high_FA_volume~treatment)
uh_high_d30_treatment_4 <- kruskal.test(data=valscombined_4site[valscombined_4site$day == "d30", ], high_FA_volume~treatment)

low_d3_location_4 <- kruskal.test(data=origvalscombined_4site[origvalscombined_4site$day == "d03", ], low_FA_volume~location)
low_d30_location_4 <- kruskal.test(data=origvalscombined_4site[origvalscombined_4site$day == "d30", ], low_FA_volume~location)
high_d3_location_4 <- kruskal.test(data=origvalscombined_4site[origvalscombined_4site$day == "d03", ], high_FA_volume~location)
high_d30_location_4 <- kruskal.test(data=origvalscombined_4site[origvalscombined_4site$day == "d30", ], high_FA_volume~location)

uh_low_d3_location_4 <- kruskal.test(data=valscombined_4site[valscombined_4site$day == "d03", ], low_FA_volume~location)
uh_low_d30_location_4 <- kruskal.test(data=valscombined_4site[valscombined_4site$day == "d30", ], low_FA_volume~location)
uh_high_d3_location_4 <- kruskal.test(data=valscombined_4site[valscombined_4site$day == "d03", ], high_FA_volume~location)
uh_high_d30_location_4 <- kruskal.test(data=valscombined_4site[valscombined_4site$day == "d30", ], high_FA_volume~location)

#adjust p values with FDR correction
loc_sites4 <- ls(pattern='location_4')
treat_sites_4 <- ls(pattern = 'treatment_4') 
loc_sites <- ls(pattern='*_*d*_location$')
treat_sites <- ls(pattern = '*d*_treatment$')

correct <- function(vecofdf){
  corpvec <- c()
  for(i in vecofdf){
    v <- eval(parse(text=i))
    print(v)
    p <- v$p.value
    corpvec <- c(corpvec, p)
  }
  corpvec <- p.adjust(corpvec)
  corpvec <- rbind(vecofdf, corpvec)
}

corloc_sites4 <- correct(loc_sites4)
cortreat_sites4 <- correct(treat_sites_4)

corloc <- correct(loc_sites)
cortreat <- correct(treat_sites)
