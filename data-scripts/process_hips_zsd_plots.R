rm(list=ls())
library(tidyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(gridExtra)

PIXEL_SIZE <- 0.108

CONFIG_PLOT <- "set2"

SIZE_CUTOFF = 230000

COLOR_PALLETE <- "Set3" #Dark2

ROOT_FOLDER <- "/Users/viana/projects/qcb/data-raw/"

setwd(ROOT_FOLDER)

#
# Functions
#

get_single_texture_feature_table <- function(table_full, feature_name, pixel_size) {
  
  table_feature <- table_full[,grep(pattern = "texture", x = names(table_full))]
  
  npixels = as.character(table_feature$dna_texture_distance_used[1])
  npixels <- gsub(pattern = "\\[", replacement = "", x = npixels)
  npixels <- gsub(pattern = "\\]", replacement = "", x = npixels)
  npixels <- as.numeric(unlist(strsplit(npixels, split=",")))
  
  subTable <- table_feature[,which(names(table_feature)==feature_name)]
  subTable <- as.character(subTable)
  subTable <- gsub(pattern = "\\[", replacement = "", x = subTable)
  subTable <- gsub(pattern = "\\]", replacement = "", x = subTable)
  subTable <- data.frame(str_split_fixed(string = subTable, pattern = ",", length(npixels)), stringsAsFactors=F)
  subTable <- as.data.frame(apply(subTable, MARGIN = 2, as.numeric))
  
  names(subTable) <- npixels
  
  aggTable <- NULL
  for (g in unique(table_full$group)) {
    
    subTable_long <- melt(subTable[which(table_full$group==g),])
    
    names(subTable_long) <- c("npixels", "feature")
    
    aggTable <- rbind(aggTable, data.frame(aggregate(feature~npixels, data=subTable_long, FUN=mean),
                                           sd=aggregate(.~npixels, data=subTable_long, FUN=sd)[,2],
                                           count=aggregate(.~npixels, data=subTable_long, FUN=length)[,2],
                                           structure=strsplit(g," ")[[1]][1],
                                           condition=strsplit(g," ")[[1]][2]))
  }
  
  aggTable$npixels <- as.numeric(as.character(aggTable$npixels))
  
  aggTable <- within(aggTable, distance<-pixel_size*npixels)
  
  return (aggTable)
}

is_different <- function(Table, num_var,fac_var) {
  n <- 1
  comp_names <- c()
  comp_pvals <- c()
  fac_levels <- levels(Table[,fac_var])
  for (i in seq(from=1, to=length(fac_levels)-1, by=1)) {
    ids_i <- which(Table[,fac_var]==fac_levels[i])
    for (j in seq(from=i+1, to=length(fac_levels), by=1)) {
      ids_j <- which(Table[,fac_var]==fac_levels[j])
      result <- t.test(
        x=Table[ids_i,num_var],
        y=Table[ids_j,num_var])
      comp_names <- c(comp_names,paste(fac_levels[i],fac_levels[j]))
      comp_pvals <- c(comp_pvals,result$p.value)
    }
  }
  print(num_var)
  print(paste(comp_names,p.adjust(comp_pvals)<0.05))
}

#
# Loading data
#

conditions <- c("H2B-TSA", "H2B+TSA", "CBX1-TSA", "CBX1+TSA")

table_full <- NULL
for (cond in conditions) {
  table_meta <- data.frame(read.table(paste("NUCLEUS_HIPS_",cond,"_meta.csv", sep=""), header=T, sep=","))
  table_feat <- data.frame(read.table(paste("NUCLEUS_HIPS_",cond,"_feature.csv", sep=""), header=T, sep=","))
  table_feat <- table_feat[,-which(names(table_feat)=="cell_id")]
  table_full <- rbind(table_full, cbind(table_meta, table_feat))
}

#
# Setting specific configuration set
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
      experiment_id=c("20181005_R02_002","20181012_R02","20181026_R02","20181005_R01_001","20181005_R01_002","20181005_R01_003","20181012_R01","20181026_R01","20181016_R02","20181019_R04","20181016_R01","20181019_R03")
    )
)

table_full$valid <- FALSE
table_full$condition <- as.character(table_full$condition)
for (exp_id in unique(table_full$czi)) {
  exp_id_char = as.character(exp_id)
  for (n in seq(1,length(config[[CONFIG_PLOT]]$labels),1)) {
    if (grepl(config[[CONFIG_PLOT]]$experiment_id[n], exp_id_char) == TRUE) {
      valids <- which(table_full$czi==exp_id)
      table_full$valid[valids] <- TRUE
      table_full$condition[valids] <- config[[CONFIG_PLOT]]$label[n]
    }
  }
}
table_full <- subset(table_full, valid==TRUE)
table_full$condition <- as.factor(table_full$condition)

#
# Pre process the table and create variable
#

table_full <- within(table_full, group <- paste(structure, condition))
table_full <- within(table_full, dna_io_intensity_ratio <- dna_io_intensity_outer_mean/dna_io_intensity_mid_mean)
table_full <- within(table_full, dna_io_intensity_ratio_slice <- dna_io_intensity_outer_slice_mean/dna_io_intensity_mid_slice_mean)
table_full <- within(table_full, dna_bright_spots_intensity_mean_norm <- dna_bright_spots_intensity_mean/dna_intensity_mean)

#
# Size cutoff (230000)
#

table_full <- subset(table_full, dna_volume>SIZE_CUTOFF)

#
# Number of cells in each population
#

table(table_full$group)

write.table(x=table_full, file=paste(ROOT_FOLDER,"../engine/data-processed/HIPS_ZSD_",CONFIG_PLOT,".csv", sep=""), row.names=F, sep=";")

##################################################
#--------------------- PLOTS ---------------------
##################################################

#
# Size and Intensity
#

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_volume, fill=condition)) +
  facet_wrap(~structure) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_volume", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_volume", fac_var="condition")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_intensity_mean, fill=condition)) +
  facet_wrap(~structure, scales="free") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_intensity_mean", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_intensity_mean", fac_var="condition")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_intensity_sum, fill=condition)) +
  facet_wrap(~structure, scales="free") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_intensity_sum", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_intensity_sum", fac_var="condition")

#
# IO Intensity
#

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_io_intensity_ratio, fill=condition)) +
  facet_wrap(~structure) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_io_intensity_ratio", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_io_intensity_ratio", fac_var="condition")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_io_intensity_ratio_slice, fill=condition)) +
  facet_wrap(~structure) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_io_intensity_ratio_slice", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_io_intensity_ratio_slice", fac_var="condition")

#
# Bright Spots
#

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_bright_spots_number, fill=condition)) +
  facet_wrap(~structure) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_bright_spots_number", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_bright_spots_number", fac_var="condition")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_bright_spots_intensity_mean, fill=condition)) +
  facet_wrap(~structure, scales="free") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_bright_spots_intensity_mean", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_bright_spots_intensity_mean", fac_var="condition")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_bright_spots_intensity_mean_norm, fill=condition)) +
  facet_wrap(~structure) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_bright_spots_intensity_mean_norm", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_bright_spots_intensity_mean_norm", fac_var="condition")

ggplot(table_full) +
  geom_point(aes(dna_volume, dna_bright_spots_number, col=condition)) +
  facet_wrap(~structure, scales="free") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_color_brewer(palette=COLOR_PALLETE)

#
# Texture
#

feature_name = c("dna_texture_haralick_contrast",
                 "dna_texture_haralick_entropy",
                 "dna_texture_haralick_variance")

feature_name = c(
                 "dna_texture_haralick_variance")

for (fname in feature_name) {
  
  result_table <- get_single_texture_feature_table(table_full=table_full, feature_name=fname, pixel_size=PIXEL_SIZE)
  
  if (CONFIG_PLOT=="set1") {
    result_table$condition <- factor(result_table$condition, levels=c("Control","DMSO","TSA"))
  } else {
    result_table$condition <- factor(result_table$condition, levels=c("Control","TSA"))
  }
  
  fig <- ggplot(result_table, aes(distance, feature, group=condition, col=condition)) +
    geom_line() +
    geom_errorbar(aes(x=distance, ymin=feature-1.96*sd/sqrt(count), ymax=feature+1.96*sd/sqrt(count)), width=0.1) +
    ylab(fname) +
    geom_vline(aes(xintercept=0.75), linetype=2) +
    facet_wrap(~structure) +
    theme_bw() +
    scale_color_brewer(palette=COLOR_PALLETE)
  
  print(fig)
}

#
# MOCK FOR HIPS and CARDIO
#

table_full_cardio <- table_full
table_full_cardio$cell_type <- "cardio"
table_full <- rbind(table_full, table_full_cardio)

ggplot(table_full) +
  geom_boxplot(aes(group, dna_volume, fill=group)) +
  facet_wrap(~cell_type) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_volume, fill=condition)) +
  facet_wrap(cell_type~structure) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_volume, fill=condition)) +
  facet_wrap(cell_type~structure) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)
