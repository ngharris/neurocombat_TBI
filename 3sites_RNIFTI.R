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
#library(ggthemes)
#library(ggpubr)
#library(ggfortify)
library(lme4)
library(grid)
library(matrixStats)
library(stats)
library(tidyverse)
#library(lsr)
library(pwr)
#library(mosaic)
#library(fMRItools)
library(ggstatsplot)
setwd("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\MD_new")


#   list of files * excludes SHM mean/SD files from FA folder
files <- list.files(path = "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\MD_new", pattern=".nii", full.names = FALSE, no.. = TRUE)
#   read NIFTI files, turn them into vectors (72 x 49 x 96 voxels = vector of length 338688)
NIFTI_df <- data.frame()
# read in NIFTI files
files <- files[!str_detect(files, "Harmonized")] #if harmonized files exist
files <- files[!str_detect(files, "significant")]
files <- files[!str_detect(files, "outlier")]
my.cluster <- makeCluster(detectCores() - 3)
clusterEvalQ(my.cluster, library(oro.nifti))
clusterExport(cl = my.cluster, varlist = c('readnii', 'files'))
doParallel::registerDoParallel(cl = my.cluster)

df <- foreach(i = 1:length(files), .combine = rbind, .multicombine = TRUE) %dopar% {
  readnii(files[i])
}

df <- as.data.frame(df)

# read in atrophy data 
setwd("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\FA_Maps_UHexpt\\FA\\Atrophy")
files_atrophy <- list.files(path = "C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\FA_Maps_UHexpt\\FA\\Atrophy", pattern = '_', full.names = FALSE, no.. = TRUE)
my.cluster <- makeCluster(detectCores() - 3)
clusterEvalQ(my.cluster, library(oro.nifti))
clusterExport(cl = my.cluster, varlist = c('readnii', 'files'))
doParallel::registerDoParallel(cl = my.cluster)

atro_df <- foreach(i = 1:length(files_atrophy), .combine = rbind, .multicombine = TRUE) %dopar% {
  read.csv(files_atrophy[i], header = FALSE, sep = " ")
}
atro_df[c('treatment', 'number', 'day', 'sex', 'random')] <- str_split_fixed(atro_df$V1, '_', 5)
atro_df_trimmed <- atro_df %>% dplyr::select('treatment', 'number', 'day', 'sex', 'V3')
atro_df_trimmed$V3 <- atro_df_trimmed$V3/1000
names(atro_df_trimmed)[names(atro_df_trimmed) == 'V3'] <- 'atrophy'
stopCluster(my.cluster)
setwd("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\MD_new")

#read in expansion data from Rachel
expansion <- read.csv("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\FA_Maps_UHexpt\\FA\\EXPANSION_all.txt", sep = " ", header = FALSE)
expansion[c('treatment', 'number', 'day', 'sex', 'random')] <- str_split_fixed(expansion$V1, '_', 5)
expansion_trimmed <- expansion %>% dplyr::select('treatment', 'number', 'day', 'sex', 'V3')
expansion_trimmed$V3 <- expansion_trimmed$V3/1000
names(expansion_trimmed)[names(expansion_trimmed) == 'V3'] <- 'expansion'

#add identifiers for harmonization: treatment, day, sex, ID (number is not super relevant here due to duplicates)
NIFTI_df <- data.frame(files, df)
NIFTI_df[c('location', 'ID')] <- str_split_fixed(NIFTI_df$files, '_', n=2)
NIFTI_df <- NIFTI_df %>% dplyr::select('location', 'ID', everything())
NIFTI_df$ID <- substring(NIFTI_df$ID, 18)
#substring(NIFTI_df$ID[NIFTI_df$location == 'USUH'], 11, 11) = '-'
NIFTI_df[c('treatment', 'number', 'day', 'sex')] <- str_split_fixed(substring(NIFTI_df$ID, 1, nchar(NIFTI_df$ID)-7), '_', 4)
NIFTI_df <- NIFTI_df %>% dplyr::select('treatment', 'number', 'day', 'sex', everything())
NIFTI_df <- NIFTI_df[NIFTI_df$location != "UUUF",]


NIFTI_df <- merge(NIFTI_df, expansion_trimmed, by = c("treatment", "number", "day"))
NIFTI_df <- NIFTI_df[, -which(names(NIFTI_df) %in% c('sex.y'))]
names(NIFTI_df)[names(NIFTI_df) == 'sex.x'] <- 'sex'

NIFTI_df$number[NIFTI_df$location == "UCLA"] <- paste0("0",NIFTI_df$number[NIFTI_df$location == "UCLA"])
NIFTI_df <- merge(NIFTI_df, atro_df_trimmed, by = c("treatment", "number", "day"))
NIFTI_df <- NIFTI_df[, -which(names(NIFTI_df) %in% c('sex.y'))]
names(NIFTI_df)[names(NIFTI_df) == 'sex.x'] <- 'sex'
#NIFTI_df <- NIFTI_df[order(NIFTI_df$files),]
NIFTI_df <- NIFTI_df %>% select(treatment, number, day, sex, atrophy, expansion, everything())
numbers <- NIFTI_df[, -1:-9]
numbers[numbers<0]<-0
attr <- NIFTI_df[, c('treatment', 'number','day','sex','ID', 'location','files')]



harmonize_tog <- function(numbers, NIFTI_df, log = TRUE){
  colnames(numbers) <- (1:ncol(numbers))
  sd_numbers <- sapply(numbers, sd)
  # remove columns that are 0 and use for harmonization
  numbers1 <- t(numbers[,sd_numbers != 0])
  #regress out atrophy (confounded with treatment and FA measurements; CCI will affect both)
  #numbers1 <- nuisance_regression(numbers1, NIFTI_df$atrophy)
  if(log == TRUE){
    numbers1 <- numbers1 + 1
    numbers1 <- log(numbers1)
  }
  # create site batches
  batch1 = NIFTI_df$location
  
  mod_tog <- model.matrix(~NIFTI_df$day+NIFTI_df$sex+NIFTI_df$treatment+NIFTI_df$atrophy)
  print(dim(numbers1))
  h_tog_NIFTI <- neuroCombat(numbers1, batch = batch1, mod = mod_tog)
  if(log == TRUE){
    h_tog_NIFTI$dat.combat <- exp(h_tog_NIFTI$dat.combat)
    h_tog_NIFTI$dat.combat <- h_tog_NIFTI$dat.combat - 1
  }
  h_tog_NIFTI$zeroes <- numbers[, which_are_constant(numbers, verbose = FALSE)]
  return(h_tog_NIFTI)
}

ks <- function(location_list, list_locations, plugin){
  combos <- combn(location_list, 2)
  plugin$location <- list_locations
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

gen_BA <- function(location_list, df, location_specific){
  combos <- combn(location_list, 2)
  color_options <- c("red", "orange", "yellow", "green", "blue", "violet")
  list_plots <- list()
  df$location <- location_specific
  print(df$location)
  for(i in 1:ncol(combos)){
    print(i)
    BA_df_pair <- cbind(as.data.frame(colMeans(df[df$location==combos[1,i],1:(ncol(df)-1)])), as.data.frame(colMeans(df[df$location==combos[2,i],1:(ncol(df)-1)])))
    #BA_df_pair <- BA_df_pair/max(abs(BA_df_pair))
    colnames(BA_df_pair) <- c(combos[1,i], combos[2,i])
    BA_df_pair$avg <- rowMeans(BA_df_pair)
    BA_df_pair$diff <- BA_df_pair[,1] - BA_df_pair[,2]
    plot1 <- ggplot(BA_df_pair, aes(x = avg, y = diff)) + geom_point(size=2, color=color_options[i]) + geom_hline(yintercept = 0) + xlim(0, .003) + ylim(-.003, .003) + ggtitle(paste(combos[1,i], combos[2,i])) + ylab("Difference Across Sites") + xlab("Average Measurement") + theme_pubr()
    list_plots[[i]] <- plot1
    
  }
  return(list_plots)
}

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

gen_Array <- function(combat_data, template, treatment_vector, location_vector){
  zero_matrix <- matrix(0, nrow(template), ncol(template))
  df1 <- data.frame(combat_data)
  stats <- t(sapply(df1, function(y) 
    unlist(summary(aov(y~treatment_vector+location_vector)))[c("Pr(>F)1","Pr(>F)2")]))
  
  stats <- apply(stats, 2, FUN = p.adjust, method = 'fdr')
  t_stats_trt <- subset(stats, stats[,1] < 0.01)
  
  s_index <- rownames(t_stats_trt)
  print(nrow(t_stats_trt))
  
  str_sub(s_index, 1, 5) <- ""
  s_index <- as.numeric(s_index)
  print(s_index)
  #print(class(s_index[1]))
  template_adjusted <- data.frame(matrix(0,nrow = nrow(template), ncol = ncol(template)))
  print(dim(template_adjusted))
  print("done1")
  template_adjusted[, c(s_index)] <- 1
  print("done2")
  #template_adjusted[, c(rownames(template_adjusted) != s_index)] <- 0
  array_adjusted <- array(colMeans(template_adjusted), dim = dim(readnii(files[1])))
  return(array_adjusted)
}

freq_sig <- function(input_things, treatment_vector,fname){
  shm_only <- input_things[treatment_vector == "SHM",]
  m_shm <- colMeans(shm_only)
  print(paste("Head of m_shm:", head(m_shm)))
  sd_shm <- apply(shm_only, MARGIN = 2, FUN = sd)
  print(paste("Head of sd_shm:", head(sd_shm)))
  m_df <- sweep(input_things, 2, m_shm, '-')
  print(paste("Dimensions of m_df:", dim(m_df)))
  z_df <- sweep(m_df, 2, sd_shm, '/') 
  print(paste("Dimensions of z_df:", dim(z_df)))
  #print(head(z_df))
  p_df <- 2 * data.frame(lapply(z_df, pnorm, mean = 0, sd = 1))
  p_df <- apply(p_df, 2, FUN = p.adjust, method = 'fdr')
  p_df[is.na(p_df) == TRUE] <- 0 #na's come from SD of shams being 0, so replacing them with a p of 0 is valid, just means CCI was much different
  p_df[p_df > 0.01] <- 5
  p_df[p_df < 0.01] <- 1
  p_df[p_df == 5] <- 0
  p_df_shm <- p_df[treatment_vector=='SHM',]
  p_df_inj <- p_df[treatment_vector!='SHM',]
  print(dim(p_df_inj))
  assign(paste0("p_df_inj_",fname), p_df_inj, .GlobalEnv)
  write.table(p_df_inj,file=paste(fname,"MD_inj_pdf.txt"),quote = FALSE,row.names = FALSE,col.names = FALSE)
  csums_p_shm <- colSums(p_df_shm)
  csums_p_shm <- csums_p_shm/nrow(shm_only)
  
  csums_p_inj <- colSums(p_df_inj)
  csums_p_inj <- csums_p_inj/nrow(input_things[treatment_vector!='SHM',])
  
  tmplt_shm <- replicate(338688, 0)
  tmplt_inj <- replicate(338688, 0)
  
  replacement_cols <- colnames(input_things)

  str_sub(replacement_cols, 1, 5) <- "" #if doing global shms change to five
  replacement_cols <- as.numeric(replacement_cols)
  tmplt_shm[replacement_cols] <- csums_p_shm
  tmplt_inj[replacement_cols] <- csums_p_inj
  
  op_array_shm <- array(tmplt_shm, dim = dim(readnii(files[1])))
  op_array_inj <- array(tmplt_inj, dim = dim(readnii(files[1])))
  
  assign(paste0("shm_",fname), op_array_shm, .GlobalEnv)
  assign(paste0("inj_",fname), op_array_inj, .GlobalEnv)
  
  writeNIfTI(op_array_shm, paste0(fname,'_shm'))
  writeNIfTI(op_array_inj, paste0(fname,'_inj'))
}    

check_change <- function(orig_array,harm_array,fname){
  sig_orig <- which(orig_array!=0)
  sig_harm <- which(harm_array!=0)
  
  new_sig <- sig_harm[!(sig_harm %in% sig_orig)]

    tmplt <- replicate(338688, 0)
  tmplt[new_sig] <- 1
  newsig_tmplt <- array(tmplt, dim = dim(readnii(files[1])))
  
  assign(paste0("newSig",fname), newsig_tmplt, .GlobalEnv)
  writeNIfTI(newsig_tmplt, paste0(fname,'_newsig'))
  
}

gen_esize <- function(input_df, attributes_vector, location_vector, locations, global_shms){
  if(global_shms == TRUE){
    shms <- input_df[attributes_vector == "SHM",] 
  } else if(global_shms == FALSE){
    shms <- input_df[attributes_vector == "SHM" & location_vector %in% locations > 0,] 
  }
  ccis <- input_df[attributes_vector == "CCI" & location_vector %in% locations > 0 , ] # dimensions check indicates this works - when using pooled locations 
  print(dim(ccis))
  print(dim(shms))
  # generate means & sd for cohen's d calculation
  shms_mean <- colMeans(shms)
  ccis_mean <- colMeans(ccis)
  pooled_stdev <- apply(input_df, 2, sd)
  
  # calculate effect sizes 
  cdvector <- (ccis_mean-shms_mean)/pooled_stdev
  
  # generate array
  tmplt1 <- replicate(338688, 0)
  replacement_cols1 <- colnames(input_df)
  print(replacement_cols1[1:5])
  str_sub(replacement_cols1, 1, 1) <- ""
  replacement_cols1 <- as.numeric(replacement_cols1)
  tmplt1[replacement_cols1] <- cdvector
  #tmplt1[abs(tmplt1) < 0.2] <- 0
  d_array <- array(tmplt1, dim = dim(readnii(files[1])))
  return(d_array)
}

gen_sd <- function(input_df){
  pooled_stdev <- apply(input_df, 2, sd)
  
  tmplt1 <- replicate(338688, 0)
  replacement_cols1 <- colnames(input_df)
  print(replacement_cols1[1:5])
  str_sub(replacement_cols1, 1, 10) <- "" # change between 1 and 10 depending on pooled vs location/day specific
  replacement_cols1 <- as.numeric(replacement_cols1)
  tmplt1[replacement_cols1] <- pooled_stdev
  #tmplt1[abs(tmplt1) < 0.2] <- 0
  d_array <- array(tmplt1, dim = dim(readnii(files[1])))
  return(d_array)
}

make_files <- function(inp, method){
  z_df <- matrix(0, nrow = ncol(inp), ncol = 338688)
  row_vector <- row.names(inp)
  #str_sub(row_vector, 1, 1) <- ""
  row_vector <- as.numeric(row_vector)
  print(dim(z_df))
  j = 1
  for (i in row_vector){
    z_df[, i] <- inp[j, 1:ncol(inp)]
    j <- j + 1
  }
  for (i in 1:nrow(z_df)){
    x = array(z_df[i,], dim = dim(readnii(files[i])))
    writeNIfTI(
      x,
      paste0(str_sub(attr$files[i],end=-7), method)
    )
  }
}


#harmonize
setwd("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\MD_new")
output <- harmonize_tog(numbers, NIFTI_df, log = FALSE)
make_files(output$dat.combat, "_Harmonized")

o_data <- data.frame(t(output$dat.original))
h_data <- data.frame(t(output$dat.combat))
o_split_by_day <- split(o_data, attr$day)
h_split_by_day <- split(h_data, attr$day)

day_3_locations <- c(unlist(split(attr$location, attr$day)[1]))
day_3_treatments <- c(unlist(split(attr$treatment, attr$day)[1]))
day_30_locations <- c(unlist(split(attr$location, attr$day)[2]))
day_30_treatments <- c(unlist(split(attr$treatment, attr$day)[2]))

#KS test
ks_test <- ks(c("UCLA", "UJHU", "USUH"), attr$location, numbers)
ks_test_harmonized <- ks(c("UCLA", "UJHU", "USUH"), attr$location, data.frame(t(output$dat.combat)))
#BA plots - for Figure 5B
original_BA_d03 <- gen_BA(c("UCLA", "UJHU", "USUH"), o_split_by_day$d03, day_3_locations)
harmonized_BA_d03 <- gen_BA(c("UCLA", "UJHU", "USUH"), h_split_by_day$d03, day_3_locations)

original_BA_d30 <- gen_BA(c("UCLA", "UJHU", "USUH"), o_split_by_day$d30, day_30_locations)
harmonized_BA_d30 <- gen_BA(c("UCLA", "UJHU", "USUH"), h_split_by_day$d30, day_30_locations)

multiplot(plotlist = original_BA_d03, layout = matrix(c(1, 2, 3), nrow = 1))
multiplot(plotlist = harmonized_BA_d03, layout = matrix(c(1, 2, 3), nrow = 1))

multiplot(plotlist = original_BA_d30, layout = matrix(c(1, 2, 3), nrow = 1))
multiplot(plotlist = harmonized_BA_d30, layout = matrix(c(1, 2, 3), nrow = 1))

#stats
bottom_array <- array(colMeans(numbers), dim = dim(readnii(files[1])))
writeNIfTI(bottom_array, "bottom_array")

original_array_d03 <- gen_Array(data.frame(o_split_by_day[1]),numbers, day_3_treatments, day_3_locations)
original_array_d30 <- gen_Array(data.frame(o_split_by_day[2]),numbers, day_30_treatments, day_30_locations)
harmonized_array_d03 <- gen_Array(data.frame(h_split_by_day[1]),numbers, day_3_treatments, day_3_locations)
harmonized_array_d30 <- gen_Array(data.frame(h_split_by_day[2]),numbers, day_30_treatments, day_30_locations)

writeNIfTI(original_array_d03, "original_array_d03")
writeNIfTI(original_array_d30, "original_array_d30")
writeNIfTI(harmonized_array_d03, "harmonized_array_d03")
writeNIfTI(harmonized_array_d30, "harmonized_array_d30")
 
#frequency data - for Figure 8
orig_freq_d03 <- freq_sig(data.frame(o_split_by_day[1]), day_3_treatments,'orig_freq_d03')
orig_freq_d30 <- freq_sig(data.frame(o_split_by_day[2]), day_30_treatments, "orig_freq_d30")
harm_freq_d03 <- freq_sig(data.frame(h_split_by_day[1]), day_3_treatments, "harm_freq_d03")
harm_freq_d30 <- freq_sig(data.frame(h_split_by_day[2]), day_30_treatments, "harm_freq_d30")


#check change
check_change(inj_orig_freq_d03,inj_harm_freq_d03,'INJchangeD3')

#e sizes - for Figure 7
orig_esize_d03_UCLA <- gen_esize(data.frame(o_split_by_day$d03), attributes_vector = day_3_treatments, location_vector = day_3_locations, locations = "UCLA", global_shms = TRUE)
orig_esize_d03_UJHU <- gen_esize(data.frame(o_split_by_day$d03), attributes_vector = day_3_treatments, location_vector = day_3_locations, locations = "UJHU", global_shms = TRUE)
orig_esize_d03_USUH <- gen_esize(data.frame(o_split_by_day$d03), attributes_vector = day_3_treatments, location_vector = day_3_locations, locations = "USUH", global_shms = TRUE)
orig_esize_d03_pooled <- gen_esize(data.frame(o_split_by_day$d03), attributes_vector = day_3_treatments, location_vector = day_3_locations, locations = c("UCLA", "UJHU", "USUH", "UNFL"), global_shms = TRUE)

harm_esize_d03_UCLA <- gen_esize(data.frame(h_split_by_day$d03), attributes_vector = day_3_treatments, location_vector = day_3_locations, locations = "UCLA", global_shms = TRUE)
harm_esize_d03_UJHU <- gen_esize(data.frame(h_split_by_day$d03), attributes_vector = day_3_treatments, location_vector = day_3_locations, locations = "UJHU", global_shms = TRUE)
harm_esize_d03_USUH <- gen_esize(data.frame(h_split_by_day$d03), attributes_vector = day_3_treatments, location_vector = day_3_locations, locations = "USUH", global_shms = TRUE)
harm_esize_d03_pooled <- gen_esize(data.frame(h_split_by_day$d03), attributes_vector = day_3_treatments, location_vector = day_3_locations, locations = c("UCLA", "UJHU", "USUH", "UNFL"), global_shms = TRUE)

orig_esize_d30_UCLA <- gen_esize(data.frame(o_split_by_day$d30), attributes_vector = day_30_treatments, location_vector = day_30_locations, locations = "UCLA", global_shms = TRUE)
orig_esize_d30_UJHU <- gen_esize(data.frame(o_split_by_day$d30), attributes_vector = day_30_treatments, location_vector = day_30_locations, locations = "UJHU", global_shms = TRUE)
orig_esize_d30_USUH <- gen_esize(data.frame(o_split_by_day$d30), attributes_vector = day_30_treatments, location_vector = day_30_locations, locations = "USUH", global_shms = TRUE)
orig_esize_d30_pooled <- gen_esize(data.frame(o_split_by_day$d30), attributes_vector = day_30_treatments, location_vector = day_30_locations, locations = c("UCLA", "UJHU", "USUH", "UNFL"), global_shms = TRUE)

harm_esize_d30_UCLA <- gen_esize(data.frame(h_split_by_day$d30), attributes_vector = day_30_treatments, location_vector = day_30_locations, locations = "UCLA", global_shms = TRUE)
harm_esize_d30_UJHU <- gen_esize(data.frame(h_split_by_day$d30), attributes_vector = day_30_treatments, location_vector = day_30_locations, locations = "UJHU", global_shms = TRUE)
harm_esize_d30_USUH <- gen_esize(data.frame(h_split_by_day$d30), attributes_vector = day_30_treatments, location_vector = day_30_locations, locations = "USUH", global_shms = TRUE)
harm_esize_d30_pooled <- gen_esize(data.frame(h_split_by_day$d30), attributes_vector = day_30_treatments, location_vector = day_30_locations, locations = c("UCLA", "UJHU", "USUH", "UNFL"), global_shms = TRUE)

esize_diff_d03_UCLA <- harm_esize_d03_UCLA-orig_esize_d03_UCLA
esize_diff_d03_UJHU <- harm_esize_d03_UJHU - orig_esize_d03_UJHU
esize_diff_d03_USUH <- harm_esize_d03_USUH - orig_esize_d03_USUH
esize_diff_d03_pooled <- harm_esize_d03_pooled - orig_esize_d03_pooled

esize_diff_d30_UCLA <- harm_esize_d30_UCLA-orig_esize_d30_UCLA
esize_diff_d30_UJHU <- harm_esize_d30_UJHU - orig_esize_d30_UJHU
esize_diff_d30_USUH <- harm_esize_d30_USUH - orig_esize_d30_USUH
esize_diff_d30_pooled <- harm_esize_d30_pooled - orig_esize_d30_pooled

writeNIfTI(orig_esize_d03_pooled, "orig_esize_d03_pooled")
writeNIfTI(orig_esize_d03_UCLA, "orig_esize_d03_UCLA")
writeNIfTI(orig_esize_d03_UJHU, "orig_esize_d03_UJHU")
writeNIfTI(orig_esize_d03_USUH, "orig_esize_d03_USUH")

writeNIfTI(orig_esize_d30_pooled, "orig_esize_d30_pooled")
writeNIfTI(orig_esize_d30_UCLA, "orig_esize_d30_UCLA")
writeNIfTI(orig_esize_d30_UJHU, "orig_esize_d30_UJHU")
writeNIfTI(orig_esize_d30_USUH, "orig_esize_d30_USUH")

writeNIfTI(harm_esize_d03_pooled, "harm_esize_d03_pooled")
writeNIfTI(harm_esize_d03_UCLA, "harm_esize_d03_UCLA")
writeNIfTI(harm_esize_d03_UJHU, "harm_esize_d03_UJHU")
writeNIfTI(harm_esize_d03_USUH, "harm_esize_d03_USUH")

writeNIfTI(harm_esize_d30_pooled, "harm_esize_d30_pooled")
writeNIfTI(harm_esize_d30_UCLA, "harm_esize_d30_UCLA")
writeNIfTI(harm_esize_d30_UJHU, "harm_esize_d30_UJHU")
writeNIfTI(harm_esize_d30_USUH, "harm_esize_d30_USUH")

writeNIfTI(esize_diff_d03_UCLA, "esize_diff_d03_UCLA")
writeNIfTI(esize_diff_d03_UJHU, "esize_diff_d03_UJHU")
writeNIfTI(esize_diff_d03_USUH, "esize_diff_d03_USUH")
writeNIfTI(esize_diff_d03_pooled, "esize_diff_d03_pooled")

writeNIfTI(esize_diff_d30_UCLA, "esize_diff_d30_UCLA")
writeNIfTI(esize_diff_d30_UJHU, "esize_diff_d30_UJHU")
writeNIfTI(esize_diff_d30_USUH, "esize_diff_d30_USUH")
writeNIfTI(esize_diff_d30_pooled, "esize_diff_d30_pooled")
#sd 
o_sd_d03 <- gen_sd(data.frame(o_split_by_day$d03))
o_sd_d30 <- gen_sd(data.frame(o_split_by_day$d30))
h_sd_d03 <- gen_sd(data.frame(h_split_by_day$d03))
h_sd_d30 <- gen_sd(data.frame(h_split_by_day$d30))

sd_diff_d03 <- h_sd_d03 - o_sd_d03
sd_diff_d30 <- h_sd_d30 - o_sd_d30
overlay(x = bottom_array, y = sd_diff_d03, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)
overlay(x = bottom_array, y = sd_diff_d30, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)


#SD differences, 3 sites
writeNIfTI(o_sd_d03, "o_sd_d03")
writeNIfTI(o_sd_d30, "o_sd_d30")
writeNIfTI(h_sd_d03, "h_sd_d03")
writeNIfTI(h_sd_d30, "h_sd_d30")
writeNIfTI(sd_diff_d03, "sd_diff_d03")
writeNIfTI(sd_diff_d30, "sd_diff_d30")

o_split_by_day_location <- split(data.frame(t(output$dat.original)), list(attr$location, attr$day))
h_split_by_day_location <- split(data.frame(t(output$dat.combat)), list(attr$location, attr$day))
treatments_by_day_site <- split(attr$treatment, list(attr$location, attr$day))
for(i in 1:length(o_split_by_day_location)){
  sd_df <- data.frame(o_split_by_day_location[i])
  sd_array <- gen_sd(sd_df)
  assign(paste0("o_",names(o_split_by_day_location)[i]), sd_array, .GlobalEnv)
}

for(i in 1:length(h_split_by_day_location)){
  sd_df <- data.frame(h_split_by_day_location[i])
  sd_array <- gen_sd(sd_df)
  assign(paste0("h_",names(o_split_by_day_location)[i]), sd_array, .GlobalEnv)
}

sd_diff_UCLA_d03 <- h_UCLA.d03 - o_UCLA.d03
writeNIfTI(sd_diff_UCLA_d03, "sd_diff_UCLA_d03")
sd_diff_UCLA_d30 <- h_UCLA.d30 - o_UCLA.d30
writeNIfTI(sd_diff_UCLA_d30, "sd_diff_UCLA_d30")
sd_diff_UJHU_d03 <- h_UJHU.d03 - o_UJHU.d03
writeNIfTI(sd_diff_UJHU_d03, "sd_diff_UJHU_d03")
sd_diff_UJHU_d30 <- h_UJHU.d30 - o_UJHU.d30
writeNIfTI(sd_diff_UJHU_d30, "sd_diff_UJHU_d30")
sd_diff_USUH_d03 <- h_USUH.d03 - o_USUH.d03
writeNIfTI(sd_diff_USUH_d03, "sd_diff_USUH_d03")
sd_diff_USUH_d30 <- h_USUH.d30 - o_USUH.d30
writeNIfTI(sd_diff_USUH_d30, "sd_diff_USUH_d30")

avg_sd_diff_d03 <- (sd_diff_UCLA_d03+sd_diff_UJHU_d03+sd_diff_USUH_d03)/3
writeNIfTI(avg_sd_diff_d03 , "avg_sd_diff_d03 ")
avg_sd_diff_d30 <- (sd_diff_UCLA_d30+sd_diff_UJHU_d30+sd_diff_USUH_d30)/3
writeNIfTI(avg_sd_diff_d30, "avg_sd_diff_d30")

#power - for Figure 6
split_by_loc <- split(numbers, NIFTI_df$location)
h_pwr_d03 <- pwr.t2n.test(d=harm_esize_d03_pooled, n1 = length(day_3_treatments == "SHM"), n2 = length(day_3_treatments == "CCI"), sig.level = 0.01, power = NULL)
o_pwr_d03 <- pwr.t2n.test(d=orig_esize_d03_pooled, n1 = length(day_3_treatments == "SHM"), n2 = length(day_3_treatments == "CCI"), sig.level = 0.01, power = NULL)
h_pwr_d30 <- pwr.t2n.test(d=harm_esize_d30_pooled, n1 = length(day_30_treatments == "SHM"), n2 = length(day_30_treatments == "CCI"), sig.level = 0.01, power = NULL)
o_pwr_d30 <- pwr.t2n.test(d=orig_esize_d30_pooled, n1 = length(day_30_treatments == "SHM"), n2 = length(day_30_treatments == "CCI"), sig.level = 0.01, power = NULL)

d03_pwr_diff <- h_pwr_d03$power - o_pwr_d03$power
overlay(x = bottom_array, y = d03_pwr_diff, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)

d30_pwr_diff <- h_pwr_d30$power - o_pwr_d30$power
overlay(x = bottom_array, y = d30_pwr_diff, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)

writeNIfTI(d03_pwr_diff, "d03_pwr_diff")
writeNIfTI(d30_pwr_diff, "d30_pwr_diff")

table(attr[c("treatment", "location", "day")])

UCLA_h_pwr_d03 <- pwr.t2n.test(d=harm_esize_d03_UCLA, n1 = 33, n2 = 40, sig.level = 0.01, power = NULL)
UCLA_h_pwr_d30 <- pwr.t2n.test(d=harm_esize_d30_UCLA, n1 = 33, n2 = 45, sig.level = 0.01, power = NULL)
UCLA_o_pwr_d03 <- pwr.t2n.test(d=orig_esize_d03_UCLA, n1 = 33, n2 = 40, sig.level = 0.01, power = NULL)
UCLA_o_pwr_d30 <- pwr.t2n.test(d=orig_esize_d30_UCLA, n1 = 33, n2 = 45, sig.level = 0.01, power = NULL)
diff_pwr_UCLA_d03 <- UCLA_h_pwr_d03$power - UCLA_o_pwr_d03$power
diff_pwr_UCLA_d30 <- UCLA_h_pwr_d30$power - UCLA_o_pwr_d30$power
overlay(x = bottom_array, y = diff_pwr_UCLA_d03, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)

writeNIfTI(diff_pwr_UCLA_d03, "diff_pwr_UCLA_d03")
writeNIfTI(diff_pwr_UCLA_d30, "diff_pwr_UCLA_d30")

UJHU_h_pwr_d03 <- pwr.t2n.test(d=harm_esize_d03_UJHU, n1 = 33, n2 = 40, sig.level = 0.01, power = NULL)
UJHU_h_pwr_d30 <- pwr.t2n.test(d=harm_esize_d30_UJHU, n1 = 31, n2 = 45, sig.level = 0.01, power = NULL)
UJHU_o_pwr_d03 <- pwr.t2n.test(d=orig_esize_d03_UJHU, n1 = 33, n2 = 40, sig.level = 0.01, power = NULL)
UJHU_o_pwr_d30 <- pwr.t2n.test(d=orig_esize_d30_UJHU, n1 = 31, n2 = 45, sig.level = 0.01, power = NULL)
diff_pwr_UJHU_d03 <- UJHU_h_pwr_d03$power - UJHU_o_pwr_d03$power
diff_pwr_UJHU_d30 <- UJHU_h_pwr_d30$power - UJHU_o_pwr_d30$power

writeNIfTI(diff_pwr_UJHU_d03, "diff_pwr_UJHU_d03")
writeNIfTI(diff_pwr_UJHU_d30, "diff_pwr_UJHU_d30")

overlay(x = bottom_array, y = diff_pwr_UJHU_d03, z = 45, plane = "coronal", plot.type = "single", col.x = grey(0:64/64), col.y = tim.colors(), NA.y = TRUE)

USUH_h_pwr_d03 <- pwr.t2n.test(d=harm_esize_d03_USUH, n1 = 33, n2 = 40, sig.level = 0.01, power = NULL)
USUH_h_pwr_d30 <- pwr.t2n.test(d=harm_esize_d30_USUH, n1 = 33, n2 = 45, sig.level = 0.01, power = NULL)
USUH_o_pwr_d03 <- pwr.t2n.test(d=orig_esize_d03_USUH, n1 = 33, n2 = 40, sig.level = 0.01, power = NULL)
USUH_o_pwr_d30 <- pwr.t2n.test(d=orig_esize_d30_USUH, n1 = 33, n2 = 45, sig.level = 0.01, power = NULL)
diff_pwr_USUH_d03 <- USUH_h_pwr_d03$power - USUH_o_pwr_d03$power
diff_pwr_USUH_d30 <- USUH_h_pwr_d30$power - USUH_o_pwr_d30$power

writeNIfTI(diff_pwr_USUH_d03, "diff_pwr_USUH_d03")
writeNIfTI(diff_pwr_USUH_d30, "diff_pwr_USUH_d30")

avg_pwr_d03 <- (diff_pwr_UCLA_d03 + diff_pwr_UJHU_d03 + diff_pwr_USUH_d03)/3
avg_pwr_d30 <- (diff_pwr_UCLA_d30 + diff_pwr_UJHU_d30 + diff_pwr_USUH_d30)/3

writeNIfTI(avg_pwr_d03, "avg_pwr_d03")
writeNIfTI(avg_pwr_d30, "avg_pwr_d30")

#correlate power with frequency
d3_pwr_freq <- data.frame(cbind(avg_pwr_d03,(orig_freq_d03)))
d3_pwr_freq <- data.frame(d3_pwr_freq[!(d3_pwr_freq[,1] == 0 & d3_pwr_freq[,2] ==0),])
colnames(d3_pwr_freq) <- c('pwr','freq')
ggscatterstats(data=d3_pwr_freq,x=freq,y=pwr,xlab='Original Frequency day 3',
               ylab = 'Change in Power',type='np')+theme(axis.text=element_text(size=30,color='black'),
                                plot.subtitle = element_text(size=40,color = 'black'),
                                axis.title = element_text(size=40))+ylim(-1,1)
ggsave('Day3pwrfreq_marge.png',width = 20,height =13)

d30_pwr_freq <- data.frame(cbind(avg_pwr_d30,(orig_freq_d30)))
d30_pwr_freq <- data.frame(d30_pwr_freq[!(d30_pwr_freq[,1] == 0 & d30_pwr_freq[,2] ==0),])
colnames(d30_pwr_freq) <- c('pwr','freq')
ggscatterstats(data=d30_pwr_freq,x=freq,y=pwr,xlab='Original Frequency day 30',ylab = 'Change in Power',
               type = 'np')+theme(axis.text=element_text(size=30,color='black'),
                                  plot.subtitle = element_text(size=40,color = 'black'),
                                  axis.title = element_text(size=40))+ylim(-1,1)
ggsave('Day30pwrfreq_marge.png',width = 20,height =13)
#Voxel Counts
#global
freq_counts <- read.csv("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\FA\\Atrophy_Metrics\\Freq_Voxel_Count.csv")
freq_counts <- freq_counts[-c(1),]
od03 <- ggplot(freq_counts, aes(x = Frequency.Bin, y = Original_d03)) + geom_col() + ylim(0,500) + xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Original d03")
od30 <- ggplot(freq_counts, aes(x = Frequency.Bin, y = Original_d30)) + geom_col()+ ylim(0,500)+ xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Original d30")
hd03 <- ggplot(freq_counts, aes(x = Frequency.Bin, y = Harmonized_d03)) + geom_col()+ ylim(0,500)+ xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Harmonized d03")
hd30 <- ggplot(freq_counts, aes(x = Frequency.Bin, y = Harmonized_d30)) + geom_col()+ ylim(0,500)+ xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Harmonized d30")

freqs <- ggarrange(od03, od30, hd03, hd30, labels = c("i", "ii", "iii", "iv"))
freqs

#site specific
for(i in 1:length(o_split_by_day_location)){
  frq_df <- data.frame(o_split_by_day_location[i])
  frq_array <- freq_sig(frq_df, c(unlist(treatments_by_day_site[i])))
  frq_array <- frq_array * nrow(frq_df)
  assign(paste0("o_",names(o_split_by_day_location)[i]), frq_array, .GlobalEnv)
}

for(i in 1:length(h_split_by_day_location)){
  frq_df <- data.frame(h_split_by_day_location[i])
  frq_array <- freq_sig(frq_df, c(unlist(treatments_by_day_site[i])))
  frq_array <- frq_array * nrow(frq_df)
  assign(paste0("h_",names(o_split_by_day_location)[i]), frq_array, .GlobalEnv)
}

o_d03_freq_site <- (o_UCLA.d03 + o_UJHU.d03 + o_USUH.d03)/281
o_d30_freq_site <- (o_UCLA.d30 + o_UJHU.d30 + o_USUH.d30)/281
h_d03_freq_site <- (h_UCLA.d03 + h_UJHU.d03 + h_USUH.d03)/281
h_d30_freq_site <- (h_UCLA.d30 + h_UJHU.d30 + h_USUH.d30)/281

writeNIfTI(o_d03_freq_site, "o_d03_freq_site")
writeNIfTI(o_d30_freq_site, "o_d30_freq_site")
writeNIfTI(h_d03_freq_site, "h_d03_freq_site")
writeNIfTI(h_d30_freq_site, "h_d30_freq_site")

freq_counts_s <- read.csv("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\FA_Maps_UHexpt\\FA\\Atrophy_Metrics\\Counts_Site.csv")
freq_counts_s <- freq_counts_s[-c(1),]
od03_s <- ggplot(freq_counts_s, aes(x = Frequency.Bin, y = Original_d03)) + geom_col() + ylim(0,650) + xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Original d03")
od30_s <- ggplot(freq_counts_s, aes(x = Frequency.Bin, y = Original_d30)) + geom_col()+ ylim(0,650)+ xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Original d30")
hd03_s <- ggplot(freq_counts_s, aes(x = Frequency.Bin, y = Harmonized_d03)) + geom_col()+ ylim(0,650)+ xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Harmonized d03")
hd30_s <- ggplot(freq_counts_s, aes(x = Frequency.Bin, y = Harmonized_d30)) + geom_col()+ ylim(0,650)+ xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Harmonized d30")


freqs_s <- ggarrange(od03_s, od30_s, hd03_s, hd30_s, labels = c("i", "ii", "iii", "iv"))
freqs_s

freq_counts <- read.csv("C:\\Users\\gkisl\\Downloads\\dataharmonizationuclajhopkins\\FA_Maps_UHexpt\\FA\\Atrophy_Metrics\\Freq_Voxel_Count.csv")
freq_counts <- freq_counts[-c(1),]
od03 <- ggplot(freq_counts, aes(x = Frequency.Bin, y = Original_d03)) + geom_col() + ylim(0,650) + xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Original d03")
od30 <- ggplot(freq_counts, aes(x = Frequency.Bin, y = Original_d30)) + geom_col()+ ylim(0,650)+ xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Original d30")
hd03 <- ggplot(freq_counts, aes(x = Frequency.Bin, y = Harmonized_d03)) + geom_col()+ ylim(0,650)+ xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Harmonized d03")
hd30 <- ggplot(freq_counts, aes(x = Frequency.Bin, y = Harmonized_d30)) + geom_col()+ ylim(0,650)+ xlab("Proportion") + ylab("Number of Voxels") + ggtitle("Harmonized d30")

freqs <- ggarrange(od03, od30, hd03, hd30, labels = c("i", "ii", "iii", "iv"))
freqs

short_file = list()
for(file in files){
  file_new <- substr(file, 1, nchar(file) - 13)
  short_file <- append(short_file, file_new)
}
u_sf <- unique(short_file)
tfCCI <- c()
for(file in u_sf){
  val <- grepl("CCI", file)
  tfCCI <- append(tfCCI, val)
}
