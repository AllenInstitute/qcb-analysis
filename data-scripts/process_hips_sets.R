rm(list=ls())
library(tidyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

ROOT_FOLDER <- "/Users/viana/projects/qcb/data-raw/"

DATASET_NAME <- "NUCLEUS_HIPS"

CONFIG_PLOT <- "set2"

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
# Setting specific set configuration
#

config <- list(
            set1=
              list(
                labels=c("Control","DMSO","TSA","TSA","TSA","Control","DMSO","TSA"),
                experiment_id=c("20181005_R02_002","20181005_R02_001","20181005_R01_001","20181005_R01_002","20181005_R01_003","20181016_R02","20181016_R03","20181016_R01")
              ),
            set2=
              list(
                labels=c("Control","Control","Control","TSA","TSA","TSA","TSA","TSA","Control","Control","TSA","TSA"),
                experiment_id=c("20181005_R02_002","20181012_R02","20181026_R02","20181005_R01_001","20181005_R01_002","20181005_R01_003","20181012_R01","20181026_R01","20181016_R02","20181019_R04","20181016_R01","20181019_R04")
              )
)

Table_raw$valid <- FALSE
Table_raw$condition <- as.character(Table_raw$condition)
for (exp_id in unique(Table_raw$czi)) {
  exp_id_char = as.character(exp_id)
  for (n in seq(1,length(config[[CONFIG_PLOT]]$labels),1)) {
    if (grepl(config[[CONFIG_PLOT]]$experiment_id[n], exp_id_char) == TRUE) {
      valids <- which(Table_raw$czi==exp_id)
      Table_raw$valid[valids] <- TRUE
      Table_raw$condition[valids] <- config[[CONFIG_PLOT]]$label[n]
    }
  }
}
Table_raw <- subset(Table_raw, valid==TRUE)
Table_raw$condition <- as.factor(Table_raw$condition)

#
# Creating new features
#

Table_raw <- within(Table_raw, dna_io_intensity_ratio <- dna_io_intensity_outer_mean/dna_io_intensity_mid_mean)
Table_raw <- within(Table_raw, dna_io_intensity_ratio_slice <- dna_io_intensity_outer_slice_mean/dna_io_intensity_mid_slice_mean)
Table_raw <- within(Table_raw, dna_intensity_sum <- dna_volume*dna_intensity_mean)
Table_raw <- within(Table_raw, dna_bright_spots_intensity_mean_norm <- dna_bright_spots_intensity_mean/dna_intensity_mean)
Table_raw <- within(Table_raw, dna_bright_spots_xy_rad <- sqrt(dna_bright_spots_xy_cross_sec_area_mean/pi))

Table_raw$date <- unlist(lapply(Table_raw$czi, function(x) strsplit(as.character(x),"_")[[1]][1]))
Table_raw$date <- paste("D", Table_raw$date, sep="")
Table_raw <- within(Table_raw, group <- paste(structure,condition))

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
Table$outlier <- ifelse(out1 + out2 + out3 > 0, "yes", "no")
Table <- subset(Table, outlier=="no")

#
# Save processed table
#

write.table(x=Table, file=paste(ROOT_FOLDER,"../engine/data-processed/",DATASET_NAME,"_set2.csv",sep=""), row.names=F, sep=";")

ggplot(Table) + geom_point(aes(dna_volume,dna_intensity_mean,col=group))
