rm(list=ls())
library(tidyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

ROOT_FOLDER <- "/Users/viana/projects/qcb/data-raw/"

DATASET_NAME <- "NUCLEUS_CARDIO"

setwd(ROOT_FOLDER)

#
# Loading data
#

Table_meta <- data.frame(read.table(paste(DATASET_NAME,"_meta.csv", sep=""), header=T, sep=","))
dim(Table_meta)
levels(Table_meta$condition)
names(Table_meta)

Table_features <- data.frame(read.table(paste(DATASET_NAME,"_feature.csv", sep=""), header=T, sep=","))
dim(Table_features)

Table_raw <- cbind(Table_meta,Table_features[,-ncol(Table_features)])

#
# Creating new features
#

Table_raw <- within(Table_raw, group <- paste(structure,condition))
Table_raw <- within(Table_raw, dna_io_intensity_ratio <- dna_io_intensity_outer_mean/dna_io_intensity_mid_mean)
Table_raw <- within(Table_raw, dna_io_intensity_ratio_slice <- dna_io_intensity_outer_slice_mean/dna_io_intensity_mid_slice_mean)
Table_raw <- within(Table_raw, dna_intensity_sum <- dna_volume*dna_intensity_mean)
Table_raw <- within(Table_raw, dna_bright_spots_intensity_mean_norm <- dna_bright_spots_intensity_mean/dna_intensity_mean)
Table_raw <- within(Table_raw, dna_bright_spots_xy_rad <- sqrt(dna_bright_spots_xy_cross_sec_area_mean/pi))

#
# Outliers detection
#

detect_outliers <- function(Table_or, fac_var, num_var, nsigma) {
  
  # Find std_min accross groups
  std_min <- c()
  for (group in unique(Table_or[,fac_var])) {
    ids <- which(Table_or[,fac_var]==group)
    std_min <- c(std_min,sd(Table_or[ids,num_var]))
  }
  std_min <- min(std_min)
  
  # Find outliers accross groups
  Table_or$tmp_outlier <- FALSE
  for (group in unique(Table_or[,fac_var])) {
    ids <- which(Table_or[,fac_var]==group)
    avg <- mean(Table_or[ids,num_var])
    
    Table_or$tmp_outlier[ids] <- ifelse( abs(Table_or[ids,num_var]-avg)>nsigma*std_min, TRUE, FALSE)
  }
  std_min <- min(std_min)
  
  return(Table_or$tmp_outlier)
}

# out1 <- detect_outliers(Table_or=Table_raw, fac_var="condition", num_var="dna_volume", nsigma=2)
# out2 <- detect_outliers(Table_or=Table_raw, fac_var="condition", num_var="dna_intensity_max", nsigma=2)
# out2 <- detect_outliers(Table_or=Table_raw, fac_var="condition", num_var="dna_roundness", nsigma=2)

Table_raw$outlier <- out1 + out2

Table <- subset(Table_raw, outlier==0)

#
# Save processed table
#

write.table(x=Table, file=paste(ROOT_FOLDER,"../engine/data-processed/",DATASET_NAME,".csv",sep=""), row.names=F, sep=";")

ggplot(Table_raw) + geom_point(aes(dna_volume,dna_roundness_roughness_xy,col=outlier))

names(Table_raw)
