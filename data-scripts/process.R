rm(list=ls())
library(tidyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

ROOT_FOLDER <- "/Users/viana/projects/qcb/data-raw/"

DATASET_NAME <- "NUCLEUS_HIPS"

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

Table_raw$date <- unlist(lapply(Table_raw$czi, function(x) strsplit(as.character(x),"_")[[1]][1]))
Table_raw$date <- paste("D", Table_raw$date, sep="")

should_not_be_used <- c("20180924_R07.czi","20180927_R04.czi","20180927_R05.czi","20180928_R03.czi","20180928_R04.czi",
                        "20181106_J01_020_Out.czi", "20181106_J02_015_Out.czi", "20181106_J03_015_Out.czi")

Table_raw <- Table_raw[is.na(match(Table_raw$czi, should_not_be_used)),]

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

Table <- Table_raw
Table$outlier <- FALSE

out1 <- detect_outliers(Table_or=Table, fac_var="condition", num_var="dna_volume", nsigma=2)
out2 <- detect_outliers(Table_or=Table, fac_var="group", num_var="dna_intensity_max", nsigma=2)
out3 <- detect_outliers(Table_or=Table, fac_var="condition", num_var="dna_roundness_roughness_xy", nsigma=3)
Table$outlier <- ifelse(out1+out2+out3>0, "yes", "no")
#Table$outlier <- ifelse(out2>0, "yes", "no")

Table <- subset(Table, outlier=="no")

#Table <- Table[-as.numeric(which(apply(is.na(Table),MARGIN = 1,FUN = sum)>0)),]

#
# Save processed table
#

write.table(x=Table, file=paste(ROOT_FOLDER,"../engine/data-processed/",DATASET_NAME,".csv",sep=""), row.names=F, sep=";")

ggplot(Table) + geom_point(aes(dna_volume,dna_intensity_mean,col=outlier))
