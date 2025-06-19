# imports
library(oro.nifti)
library(neurobase)
library(stringr)
library(neuroCombat)
library(dataPreparation)
library(dplyr)
library(stringr)
library(parallel)
library(effectsize)
library(foreach)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggfortify)
library(lme4)
library(grid)
library(matrixStats)
library(stats)
library(tidyverse)
library(lsr)
library(pwr)

setwd("")

#   list of files * excludes SHM mean/SD files from FA folder
files <- list.files(path = "", pattern=".nii", full.names = FALSE, no.. = TRUE)
#   read NIFTI files, turn them into vectors (72 x 49 x 96 voxels = vector of length 338688)
NIFTI_df <- data.frame()
# read in NIFTI files
files <- files[!str_detect(files, "_Harmonized")] #if harmonized files exist
files <- files[!str_detect(files, "significant")]
files <- files[!str_detect(files, "outlier")]
my.cluster <- makeCluster(detectCores() - 3)
clusterEvalQ(my.cluster, library(oro.nifti))
clusterExport(cl = my.cluster, varlist = c('readnii', 'files'))
doParallel::registerDoParallel(cl = my.cluster)

st <- Sys.time()
df <- foreach(i = 1:length(files), .combine = rbind, .multicombine = TRUE) %dopar% {
  #print(i)
  readnii(files[i])
  #assign(paste0("X_",files[i]), x , .GlobalEnv)
  #NIFTI_df <- rbind(NIFTI_df, x)
}
df
df <- as.data.frame(df)
et <- Sys.time()
print(et-st)
stopCluster(my.cluster)

#add identifiers for harmonization: treatment, day, sex, ID (number is not super relevant here due to duplicates)
NIFTI_df <- data.frame(files, df)
NIFTI_df[c('location', 'ID')] <- str_split_fixed(NIFTI_df$files, '_', n=2)
NIFTI_df <- NIFTI_df %>% dplyr::select('location', 'ID', everything())
NIFTI_df$ID <- substring(NIFTI_df$ID, 18)
substring(NIFTI_df$ID[NIFTI_df$location == 'USUH'], 11, 11) = '-'
NIFTI_df[c('treatment', 'number', 'day', 'sex')] <- str_split_fixed(substring(NIFTI_df$ID, 1, nchar(NIFTI_df$ID)-7), '_', 4)
NIFTI_df <- NIFTI_df %>% dplyr::select('treatment', 'number', 'day', 'sex', everything())
numbers <- NIFTI_df[, -1:-7]

#perform kalmogorov smirnov test 
ks <- function(location_list, plugin){
  combos <- combn(location_list, 2)
  plugin$location <- rep(c("UCLA", "UJHU", "USUH", "UNFL"), times=c(99, 89, 96, 53))
  print(ncol(plugin))
  lpoint <- ncol(plugin) -1
  print(lpoint)
  plugin$avgks <- rowMeans(plugin[,1:lpoint])
  ksstats <- data.frame()
  for(i in 1:ncol(combos)){
    #print(i)
    #print(df[df$location ==combos[1,i],338690])
    ksres <- ks.test(plugin[plugin$location==combos[1,i],ncol(plugin)], plugin[plugin$location==combos[2,i],ncol(plugin)])
    ksfilt <- unlist(ksres)[c("p.value")]
    ksstats <- rbind(ksstats, ksfilt)
  }
  return(ksstats)
}
ks_test <- ks(c("UCLA", "UJHU", "USUH", "UNFL"), numbers)
#harmonizing NIFTI data, voxel-wise (each column is a voxel)
harmonize_tog <- function(numbers, NIFTI_df, log = TRUE){
  colnames(numbers) <- (1:ncol(numbers))
  # remove columns that are 0 and use for harmonization
  numbers1 <- t(numbers[,-which_are_constant(numbers, verbose = FALSE)])
  if(log == TRUE){
    numbers1 <- numbers1 + 1
    numbers1 <- log(numbers1)
  }
  # create site batches
  batch1 = NIFTI_df$location
  
  mod_tog <- model.matrix(~NIFTI_df$day+NIFTI_df$treatment+NIFTI_df$sex)
  print(dim(numbers1))
  h_tog_NIFTI <- neuroCombat(numbers1, batch = batch1, mod = mod_tog)
  if(log == TRUE){
    h_tog_NIFTI$dat.combat <- exp(h_tog_NIFTI$dat.combat)
    h_tog_NIFTI$dat.combat <- h_tog_NIFTI$dat.combat - 1
  }
  h_tog_NIFTI$zeroes <- numbers[, which_are_constant(numbers, verbose = FALSE)]
  return(h_tog_NIFTI)
}

#numbers1 <- t(numbers[,-which_are_constant(numbers, verbose = FALSE)])
output <- harmonize_tog(numbers, NIFTI_df, log = FALSE)
h_data <- as.data.frame(output$dat.combat)

new_NIFTI <- matrix(ncol = ncol(output$dat.combat), nrow = ncol(numbers))
new_NIFTI[output$zeroes,] <- c(rep(0, ncol(output$dat.combat) ))
new_NIFTI[setdiff(c(1:nrow(new_NIFTI)),output$zeroes),] <- output$dat.combat

# recreate NIFTI files
make_files <- function(inp, method){
  z_df <- matrix(0, nrow = ncol(inp), ncol = 338688)
  row_vector <- row.names(inp)
  #str_sub(row_vector, 1, 1) <- ""
  row_vector <- as.numeric(row_vector)
  print(row_vector[1])
  j = 1
  for (i in row_vector){
    z_df[, i] <- inp[j, 1:337]
    j <- j + 1
  }
  for (i in 1:nrow(z_df)){
    x = array(z_df[i,], dim = dim(readnii(files[i])))
    writeNIfTI(
      x,
      paste0(str_sub(files[i], end = -7), method)
    )
  }
}

make_files(output$dat.combat, "_Harmonizednlog")

# read harmonized data
h_files <- list.files(path = "", pattern="_Harmonizednlog", full.names = FALSE, no.. = TRUE)
harmonized_NIFTI_df <- data.frame()
#for(i in h_files){
#print(i)
#x= c(readnii(i))
#assign(paste0("X_",i), x , .GlobalEnv)
#harmonized_NIFTI_df <- rbind(harmonized_NIFTI_df, x)
#}

cl <- makeCluster(detectCores()-3)
clusterEvalQ(cl, library(oro.nifti))
clusterExport(cl = cl, varlist = c('readnii', 'h_files'))
doParallel::registerDoParallel(cl = cl)

h_df <- foreach(i = 1:length(h_files), .combine = rbind, .multicombine = TRUE) %dopar% {
  readnii(h_files[i])
}
h_df
h_df <- data.frame(h_files, h_df)
stopCluster(cl)
harmonized_NIFTI_df <- h_df

# Compute mean FA (for comparison to harmonized means) and extract only means (faster processing)
attr <- NIFTI_df[, c('treatment', 'number','day','sex','ID')]
meanFA <- rowMeans(harmonized_NIFTI_df[, 2:338688])
avg_NIFTI_df <- cbind(meanFA, attributes) # harmonized data
colnames(avg_NIFTI_df)[1] <- "meanFA"
orig_NIFTI_df <- cbind(rowMeans(numbers[, 2:338688]), attributes) # original data
colnames(orig_NIFTI_df)[1] <- "meanFA"

# effect size
effect_size_original <- (mean(orig_NIFTI_df[,1][which(orig_NIFTI_df$treatment == "CCI")]) - mean(orig_NIFTI_df[,1][which(orig_NIFTI_df$treatment == "SHM")]))/sd_pooled(orig_NIFTI_df[,1])
effect_size <- (mean(avg_NIFTI_df[,1][which(avg_NIFTI_df$treatment == "CCI")]) - mean(avg_NIFTI_df[,1][which(avg_NIFTI_df$treatment == "SHM")]))/sd_pooled(avg_NIFTI_df[,1])

all_NIFTI <- rbind(orig_NIFTI_df, avg_NIFTI_df)
all_NIFTI$harmonized <- rep(c("Original", "Harmonized"), each=nrow(all_NIFTI)/2)

# add locations
all_NIFTI$location <- rep(c("UCLA", "UJHU", "USUH", "UNFL"), times=c(99, 89, 96, 53))

# show plot of comparison between voxel and mean
modify_df <- function(ldf){
  ldf$treatLoc[ldf$treatment=="CCI" & ldf$location=="UCLA"] <- "CCI.Site1"
  ldf$treatLoc[ldf$treatment=="CCI" & ldf$location=="UJHU"] <- "CCI.Site2"
  ldf$treatLoc[ldf$treatment=="CCI" & ldf$location=="USUH"] <- "CCI.Site3"
  ldf$treatLoc[ldf$treatment=="CCI" & ldf$location=="UNFL"] <- "CCI.Site4"
  ldf$treatLoc[ldf$treatment=="SHM" & ldf$location=="UCLA"] <- "SHM.Site1"
  ldf$treatLoc[ldf$treatment=="SHM" & ldf$location=="UJHU"] <- "SHM.Site2"
  ldf$treatLoc[ldf$treatment=="SHM" & ldf$location=="USUH"] <- "SHM.Site3"
  ldf$treatLoc[ldf$treatment=="SHM" & ldf$location=="UNFL"] <- "SHM.Site4"
  return(ldf)
}

modify_df_means <- function(ldf){
  ldf$treatLoc[ldf$treatment=="CCI" & ldf$location==1] <- "CCI.Site1"
  ldf$treatLoc[ldf$treatment=="CCI" & ldf$location==2] <- "CCI.Site2"
  ldf$treatLoc[ldf$treatment=="CCI" & ldf$location==3] <- "CCI.Site3"
  ldf$treatLoc[ldf$treatment=="SHM" & ldf$location==1] <- "SHM.Site1"
  ldf$treatLoc[ldf$treatment=="SHM" & ldf$location==2] <- "SHM.Site2"
  ldf$treatLoc[ldf$treatment=="SHM" & ldf$location==3] <- "SHM.Site3"
  return(ldf)
}

all_NIFTI_mod <- modify_df(all_NIFTI)


#just imaging harmonization
violinplot1 <- ggplot(all_NIFTI_mod, aes(x=harmonized, y=meanFA, fill=treatLoc)) + geom_violin(position=position_dodge(0.8), width=1.1, alpha=0.75) + 
  geom_boxplot(position=position_dodge(0.8), width=0.2) + 
  stat_summary(fun=median, show.legend = FALSE, geom="crossbar", position=position_dodge(0.8), width=0.5) +
  labs(x="", y="FA", title="Fractional Anisotropy") +
  theme_classic() + scale_x_discrete(limits = c("Original","Harmonized")) + 
  scale_fill_manual(values=c("red", "orange", "yellow", "green", "blue", "violet", "brown", "black"), name="Treatment and Location", labels=c("TBI Site 1", "TBI Site 2", "TBI Site 3", "TBI Site 4", "Sham Site 1", "Sham Site 2", "Sham Site 3", "Sham Site 4"))


# B-A plot --> want each voxel to be a point and find the mean and difference of each voxel
NIFTI_df$location <- rep(c("UCLA", "UJHU", "USUH", "UNFL"), times=c(99, 89, 96, 53))
UCLA <- as.data.frame(colMeans(NIFTI_df[NIFTI_df$location=="UCLA", 8:338695]))
UJHU <- as.data.frame(colMeans(NIFTI_df[NIFTI_df$location=="UJHU", 8:338695]))
USUH <- as.data.frame(colMeans(NIFTI_df[NIFTI_df$location=="USUH", 8:338695]))
UNFL <- as.data.frame(colMeans(NIFTI_df[NIFTI_df$location=="UNFL", 8:338695]))

gen_BA <- function(location_list, df){
  combos <- combn(location_list, 2)
  color_options <- c("red", "orange", "yellow", "green", "blue", "violet")
  list_plots <- list()
  df$location <- rep(c("UCLA", "UJHU", "USUH", "UNFL"), times=c(99, 89, 96, 53))
  for(i in 1:ncol(combos)){
    print(i)
    BA_df_pair <- cbind(as.data.frame(colMeans(df[df$location==combos[1,i],1:338687])), as.data.frame(colMeans(df[df$location==combos[2,i],1:338687])))
    colnames(BA_df_pair) <- c(combos[1,i], combos[2,i])
    BA_df_pair$avg <- rowMeans(BA_df_pair)
    BA_df_pair$diff <- BA_df_pair[,1] - BA_df_pair[,2]
    plot1 <- ggplot(BA_df_pair, aes(x = avg, y = diff)) + geom_point(size=2, color=color_options[i]) + geom_hline(yintercept = 0) + xlim(0, 1) + ylim(-1, 1) + ggtitle(paste(combos[1,i], combos[2,i])) + ylab("Difference Across Sites") + xlab("Average Measurement") + theme_pubr()
    list_plots[[i]] <- plot1
    
  }
  return(list_plots)
}

original_BA <- gen_BA(c("UCLA", "UJHU", "USUH", "UNFL"), NIFTI_df[,8:338695])
#multiplot(plotlist = original_BA, cols = 1)
harmonized_BA <- gen_BA(c("UCLA", "UJHU", "USUH", "UNFL"), harmonized_NIFTI_df[,2:338688])

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#to generate Bland-Altman plots
multiplot(plotlist = original_BA, layout = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2))
multiplot(plotlist = harmonized_BA, layout = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2))


# inputs: harmonized combat data directly from neuroCombat output, template (usually harmonized data
# typically 338688 elements for FA data), make sure it is the correct size)
# and treatment_vector, must be same length as number of columns in template

#output: array of significant voxels based on two way ANOVA results. Useful for visualization with overlay function
gen_Array <- function(combat_data, template, treatment_vector){
  zero_matrix <- matrix(0, nrow(template), ncol(template))
  df1 <- data.frame(t(combat_data))
  stats <- t(sapply(df1, function(y) 
    unlist(summary(aov(y~attributes$treatment+locs)))[c("Pr(>F)1","Pr(>F)2")]))
  
  stats <- apply(stats, 2, FUN = p.adjust, method = 'fdr')
  t_stats_trt <- subset(stats, stats[,1] < 0.01)
  
  s_index <- rownames(t_stats_trt)
  print(nrow(t_stats_trt))
  
  str_sub(s_index, 1, 1) <- ""
  s_index <- as.numeric(s_index)
  template_adjusted <- template
  template_adjusted[, s_index] <- 1
  template_adjusted[, !(seq.int(1, ncol(template)) %in% s_index)] <- 0
  array_adjusted <- array(colMeans(template_adjusted), dim = dim(readnii(files[1])))
  return(array_adjusted)
}

#for outlier counting
detect_outlier <- function(x) {
  
  # calculate first quantile
  Quantile1 <- quantile(x, probs=.25)
  
  # calculate third quantile
  Quantile3 <- quantile(x, probs=.75)
  
  # calculate inter quartile range
  IQR = Quantile3-Quantile1
  
  # return true or false
  x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}
gen_outlier_score <- function(input_df){
  Q1 <- sapply(input_df, stats::quantile, probs = 0.25, na.rm = TRUE)
  print(Q1[1:5])
  Q3 <- sapply(input_df, stats::quantile, probs = 0.75, na.rm = TRUE)
  print(Q3[1:5])
  InterQR <- sapply(input_df, IQR, na.rm = FALSE)
  print(class(InterQR))
  print(InterQR[1:5])
  os_1 <- abs(sweep(input_df, 2, Q3, "-")) # subtract Q3 to get smallest value, since all values should be positive at this point
  outlier_score <- sweep(input_df, 2, InterQR, "/")
  print("done")
  return(outlier_score)
}


#checking site effects - harmonized data
locs <- rep(c("UCLA", "UJHU", "USUH", "UNFL"), times=c(99, 89, 96, 53))
sig_site <- function(df){
  site_stuff <- t(sapply(data.frame(t(df)), function(y) 
    unlist(summary(aov(y~attributes$treatment+locs)))[c("Pr(>F)1","Pr(>F)2")]))
}

#check site effects  
test_site <- sig_site(output$dat.combat)
test_site <- apply(test_site, 2, FUN = p.adjust, method = 'fdr')
test_site_filtered <- subset(test_site, test_site[,2] < 0.01)
test_loc_filtered <- subset(test_site, test_site[,1] < 0.01)

#checking site effects - original data
o_site <- sig_site(t(numbers))
o_site <- na.omit(o_site)
o_site <- apply(o_site, 2, FUN = p.adjust, method = 'fdr')
o_site_filtered <- subset(o_site, o_site[,2] < 0.01)
o_loc_filtered <- subset(o_site, o_site[,1] < 0.01)


#make bar chart
list_values <- c(nrow(o_site_filtered), nrow(test_site_filtered))
list_names <- c("Raw", "Harmonized")

bplot <- barplot(list_values, names.arg = list_names, xlab = "Method", ylab = "Number of Voxels", ylim = c(0, 200000), main = "Voxels Associated with Difference in Location", col = c("cyan", "blue", "purple", "violet"))
text(bplot, list_values, list_values, pos = 3)

treatment_values <- c(nrow(o_loc_filtered), nrow(test_loc_filtered))
treatment_plot <- barplot(treatment_values, names.arg = list_names, xlab = "Method", ylab = "Number of Voxels", ylim = c(0, 4000), main = "Voxels Associated with Difference in Group", col = c("cyan", "blue", "purple", "violet"))
text(treatment_plot, treatment_values, treatment_values, pos = 3)

bottom_array <- array(colMeans(numbers), dim = dim(readnii(files[1])))
original_array <- gen_Array(t(numbers), numbers, NIFTI_df$treatment)
basic_array <- gen_Array(output$dat.combat, numbers, NIFTI_df$treatment)


make_files(sig_output$dat.combat, "significantlog")
make_files(outliers_output$dat.combat, "outlierlog ")

sig_files <- list.files(path = "", pattern="significant", full.names = FALSE, no.. = TRUE)
#for(i in h_files){
#print(i)
#x= c(readnii(i))
#assign(paste0("X_",i), x , .GlobalEnv)
#harmonized_NIFTI_df <- rbind(harmonized_NIFTI_df, x)
#}

my.cluster <- makeCluster(detectCores() - 3)
clusterEvalQ(my.cluster, library(oro.nifti))
clusterExport(cl = cl, varlist = c('readnii', 'sig_files'))
doParallel::registerDoParallel(cl = cl)

check_sig_df <- foreach(i = 1:length(sig_files), .combine = rbind, .multicombine = TRUE) %dopar% {
  readnii(sig_files[i])
}
check_sig_df
check_sig_df <- data.frame(sig_files, check_sig_df)
stopCluster(cl)

overlay(bottom_array, original_array, z = 45, plane = "coronal", plot.type = "single" , col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)
overlay(x = bottom_array, y = basic_array, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)

#density plots for voxel data
movement_data <- cbind(test_site[,1], o_site[,1])
dh <- density(movement_data[,1])
do <- density (movement_data[,2])

plot(dh, lwd = 2, col = "blue", main = "FDR adjusted p-values by Group", ylim = c(0, 5))
lines(do, lwd = 2, col = "green")
legend(0, 5, legend=c("Harmonized", "Original"), fill = c("blue", "green"))

#frequencies of significant voxels 
#input_things is some kind of neuroCombat output - df$dat.combat
freq_sig <- function(input_things){
  shm_only <- input_things[attr$treatment == "SHM",]
  m_shm <- colMeans(shm_only)
  print(paste("Head of m_shm:", head(m_shm)))
  sd_shm <- colSds(shm_only)
  print(paste("Head of sd_shm:", head(sd_shm)))
  m_df <- sweep(input_things, 2, m_shm, '-')
  print(paste("Dimensions of m_df:", dim(m_df)))
  z_df <- sweep(m_df, 2, sd_shm, '/')
  print(paste("Dimensions of z_df:", dim(z_df)))
  p_df <- pnorm(z_df, mean = 0, sd = 1)
  p_df <- apply(p_df, 2, FUN = p.adjust, method = 'fdr')
  print(dim(p_df))
  p_df[p_df > 0.01] <- 5
  p_df[p_df < 0.01] <- 1
  p_df[p_df == 5] <- 0
  csums_p <- colSums(p_df)
  csums_p <- csums_p/nrow(input_things)
  tmplt <- replicate(338688, 0)
  replacement_cols <- as.numeric(colnames(input_things))
  tmplt[replacement_cols] <- csums_p
  print(mean(tmplt))
  op_array <- array(tmplt, dim = dim(readnii(files[1])))
  return(op_array)
}     
harm_pic <- freq_sig(t(output$dat.combat))
orig_pic <- freq_sig(t(output$dat.original))
overlay(x = bottom_array, y = harm_pic, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)
overlay(x = bottom_array, y = orig_pic, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)

# effect size by voxel 
gen_esize <- function(input_df, location){
  shms <- input_df[attr$treatment == "SHM" & attr$location == location,]
  # if you want site-specific effect sizes, then keep the & attr$location == location, remove if you want cohen's d relative to all shams
  ccis <- input_df[attr$treatment == "CCI" & attr$location == location, ]
  # generate means & sd for cohen's d calculation
  shms_mean <- colMeans(shms)
  ccis_mean <- colMeans(ccis)
  pooled_stdev <- colSds(input_df)

  # calculate effect sizes 
  cdvector <- (ccis_mean-shms_mean)/pooled_stdev

  # generate array
  tmplt1 <- replicate(338688, 0)
  replacement_cols1 <- as.numeric(colnames(input_df))
  tmplt1[replacement_cols1] <- cdvector
  #tmplt1[abs(tmplt1) < 0.2] <- 0
  d_array <- array(tmplt1, dim = dim(readnii(files[1])))
  return(d_array)
}

orig_d <- gen_esize(t(output$dat.original), location = "UNFL")
harm_d <- gen_esize(t(output$dat.combat), location = "UNFL")
overlay(x = bottom_array, y = orig_d, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)
overlay(x = bottom_array, y = harm_d, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)

# write arrays to NIFTI files
writeNIfTI(orig_d, "original_cohensd_specificUNFL")
writeNIfTI(harm_d, "harmonized_cohensd_specificUNFL")

test_image <- readnii("original_cohensd.nii.gz")
test_image2 <- readnii("harmonized_cohensd.nii.gz")

# density and violin plots to show effect size movement
harm_vector <- as.vector(harm_d[abs(harm_d) > 0.2])
orig_vector <- as.vector(orig_d[abs(orig_d) > 0.2])

hd <- density(harm_vector)
od <- density(orig_vector)

#in case you want to filter for certain effect sizes
plot(hd, lwd = 2, col = "blue", main = ">0.2 Cohen's d Values", ylim = c(0, 7))
lines(od, lwd = 2, col = "green")
legend(0, 5, legend=c("Harmonized", "Original"), fill = c("blue", "green"), cex = 0.5)

#comparing repeated measures anova and mixed methods design
sd_numbers <- apply(numbers, 2, sd)
no_zero_NIFTI <- cbind(NIFTI_df[,1:7], numbers[,sd_numbers != 0])
no_zero_NIFTI$treatment <- as.factor(no_zero_NIFTI$treatment)
no_zero_NIFTI$files <- as.factor(no_zero_NIFTI$files)
no_zero_NIFTI$location <- as.factor(no_zero_NIFTI$location)
#no_zero_NIFTI$files <- substr(no_zero_NIFTI$files, 1, nchar(no_zero_NIFTI$files) - 13)

  #run repeated measures
no_zero_colnames <- colnames(no_zero_NIFTI[,8:110257])
no_zero_NIFTI_g <- no_zero_NIFTI %>%
  pivot_longer(cols =  8:110257, names_to = "voxel index", values_to = "FA")
results_set <- list()
# original data
aov_res <- sapply(no_zero_NIFTI[8:110257], function(y) 
  unlist(summary(aov(y~((no_zero_NIFTI$treatment*no_zero_NIFTI$location*no_zero_NIFTI$day) 
                        + Error(no_zero_NIFTI$files/no_zero_NIFTI$day)))))[c("Error: no_zero_NIFTI$files.Pr(>F)1","Error: no_zero_NIFTI$files.Pr(>F)2", "Error: no_zero_NIFTI$files.Pr(>F)3", "Error: no_zero_NIFTI$files.Sum Sq1", "Error: no_zero_NIFTI$files.Sum Sq2", "Error: no_zero_NIFTI$files.Sum Sq3", "Error: no_zero_NIFTI$files.Sum Sq4", "Error: no_zero_NIFTI$files.Sum Sq5", "Error: no_zero_NIFTI$files.Sum Sq6", "Error: no_zero_NIFTI$files.Sum Sq7", "Error: no_zero_NIFTI$files.Sum Sq8")])
aov_res <- t(aov_res)
colnames(aov_res) <- c("treatment", "location", "day", "SS_treatment", "SS_location", "SS_day", "SS_treatment:location", "SS_treatment_day", "SS_location:day", "SS_treatment:day:location", "SS_Res")
aov_ss_sums <- rowSums(aov_res[,4:11])
aov_eta_test <- aov_res[,4:11]/aov_ss_sums
aov_eta_test_means_o <- colMeans(aov_eta_test)

#harmonized data
h_output <- data.frame(t(output$dat.combat))
aov_res_h <- sapply(h_output, function(y) 
  unlist(summary(aov(y~((no_zero_NIFTI$treatment*no_zero_NIFTI$location*no_zero_NIFTI$day) 
                        + Error(no_zero_NIFTI$files/no_zero_NIFTI$day)))))[c("Error: no_zero_NIFTI$files.Pr(>F)1","Error: no_zero_NIFTI$files.Pr(>F)2", "Error: no_zero_NIFTI$files.Pr(>F)3", "Error: no_zero_NIFTI$files.Sum Sq1", "Error: no_zero_NIFTI$files.Sum Sq2", "Error: no_zero_NIFTI$files.Sum Sq3", "Error: no_zero_NIFTI$files.Sum Sq4", "Error: no_zero_NIFTI$files.Sum Sq5", "Error: no_zero_NIFTI$files.Sum Sq6", "Error: no_zero_NIFTI$files.Sum Sq7", "Error: no_zero_NIFTI$files.Sum Sq8")])
aov_res_h <- t(aov_res_h)
colnames(aov_res_h) <- c("treatment", "location", "day", "SS_treatment", "SS_location", "SS_day", "SS_treatment:location", "SS_treatment_day", "SS_location:day", "SS_treatment:day:location", "SS_Res")
aov_ss_sums_h <- rowSums(aov_res_h[,4:11])
aov_eta_test_h <- aov_res_h[,4:11]/aov_ss_sums_h
aov_eta_test_means_h <- colMeans(aov_eta_test_h)

# look at AIC and BIC
model_res <- sapply(no_zero_NIFTI[8:110257], function(y) 
  unlist(glance(lm(y~(no_zero_NIFTI$treatment*no_zero_NIFTI$location*no_zero_NIFTI$day 
                        ))))[c("AIC", "BIC")])
model_res1 <- sapply(no_zero_NIFTI[8:110257], function(y) 
  unlist(glance(lm(y~(no_zero_NIFTI$treatment*no_zero_NIFTI$location
  ))))[c("AIC", "BIC")])

#calculate power
h_pwr <- pwr.t2n.test(h=harm_d, n1 = 233, n2 = 104, sig.level = 0.01, power = NULL)
o_pwr <- pwr.t2n.test(h=orig_d, n1 = 233, n2 = 104, sig.level = 0.01, power = NULL)

#rough estimates of power before/after harmonization - not accurate given that original data should be separated by site
# this can be done with the gen_esize function and setting location to get site specific effect sizes, which can then be 
# used here. 
h_pwr_map <- overlay(x = bottom_array, y = h_pwr$power, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)
o_pwr_map <- overlay(x = bottom_array, y = o_pwr$power, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)

#outlier counting
count_outliers <- 0
df_outliers <- numbers
h_df_outliers <- harmonized_NIFTI_df[,-1]
h_df_outliers <- data.frame(h_df_outliers[,colSds(as.matrix(h_df_outliers)) != 0])
df_outliers <- data.frame(numbers[,colSds(as.matrix(numbers)) != 0])
outlier_scores <- gen_outlier_score(df_outliers)
detected_outliers <- sapply(df_outliers, detect_outlier)
h_detected_outliers <- sapply(h_df_outliers, detect_outlier)

UUUF_os <- outlier_scores[attr$location == "UUUF",]
UUUF_outlier_index <- mean(detected_outliers[attr$location == "UUUF",])
h_UUUF_outlier_index <- mean(h_detected_outliers[attr$location == "UUUF",])

UUUF_os_true <- UUUF_os[UUUF_outlier_index == TRUE]
UUUF_os_true[sapply(UUUF_os_true, is.infinite)] <-NA
UUUF_os_true <- na.omit(UUUF_os_true)

USUH_os <- outlier_scores[attr$location == "USUH",]
USUH_outlier_index <- mean(detected_outliers[attr$location == "USUH",])
h_USUH_outlier_index <- mean(h_detected_outliers[attr$location == "USUH",])

USUH_os_true <- USUH_os[USUH_outlier_index == TRUE]
USUH_os_true[sapply(USUH_os_true, is.infinite)] <-NA
USUH_os_true <- na.omit(USUH_os_true)

UJHU_os <- outlier_scores[attr$location == "UJHU",]
UJHU_outlier_index <- mean(detected_outliers[attr$location == "UJHU",])
h_UJHU_outlier_index <- mean(h_detected_outliers[attr$location == "UJHU",])

UJHU_os_true <- UJHU_os[UJHU_outlier_index == TRUE]
UJHU_os_true[sapply(UJHU_os_true, is.infinite)] <-NA
UJHU_os_true <- na.omit(UJHU_os_true)

UCLA_os <- outlier_scores[attr$location == "UCLA",]
UCLA_outlier_index <-  mean(detected_outliers[attr$location == "UCLA",])
h_UCLA_outlier_index <- mean(h_detected_outliers[attr$location == "UCLA",])

UCLA_os_true <- UCLA_os[UCLA_outlier_index == TRUE]
UCLA_os_true[sapply(UCLA_os_true, is.infinite)] <-NA
UCLA_os_true <- na.omit(UCLA_os_true)