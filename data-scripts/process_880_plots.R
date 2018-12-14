rm(list=ls())
library(tidyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(gridExtra)

PIXEL_SIZE <- 0.05

CONFIG_PLOT <- "set2"

SIZE_CUTOFF = list(hiPS=0*(PIXEL_SIZE**3),
                   cardiomyocyte=0*(PIXEL_SIZE**3))

COLOR_PALLETE <- "Dark2" #Dark2

ROOT_FOLDER <- "/home/matheus.viana/projects/qcb/qcb-analysis/data-raw/"

setwd(ROOT_FOLDER)

#
# Functions
#

get_single_texture_feature_table <- function(table_full, feature_name, pixel_size, distances) {
  
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
    
    if (length(distances) != length(npixels)) {
      
      aggTable <- rbind(aggTable, data.frame(subTable_long,
                                             structure=strsplit(g," ")[[1]][1],
                                             condition=strsplit(g," ")[[1]][2],
                                             cell_type=strsplit(g," ")[[1]][3]))
      
    } else {
      
      aggTable <- rbind(aggTable, data.frame(aggregate(feature~npixels, data=subTable_long, FUN=mean),
                                             sd=aggregate(.~npixels, data=subTable_long, FUN=sd)[,2],
                                             count=aggregate(.~npixels, data=subTable_long, FUN=length)[,2],
                                             structure=strsplit(g," ")[[1]][1],
                                             condition=strsplit(g," ")[[1]][2],
                                             cell_type=strsplit(g," ")[[1]][3]))
    }
    
  }
  
  aggTable$npixels <- as.numeric(as.character(aggTable$npixels))
  
  aggTable <- within(aggTable, distance<-pixel_size*npixels)
  
  return (aggTable)
}


is_different_anova <- function(Table, num_var,fac_var) {
  aov_table <- aov(Table[,which(names(Table)==num_var)]~Table[,which(names(Table)==fac_var)])
  print(summary(aov_table))
  result <- TukeyHSD(aov_table, p.adjust.method ="none")$`Table[, which(names(Table) == fac_var)]`
  print(result)
}

#
# Loading data
#

conditions <- c("NUCLEUS_HIPS_880", "NUCLEUS_HIPS_H2B_880", "NUCLEUS_CARDIO_880")

table_full <- NULL
for (cond in conditions) {
  table_meta <- data.frame(read.table(paste(cond,"_meta.csv", sep=""), header=T, sep=",", stringsAsFactors=FALSE))
  table_feat <- data.frame(read.table(paste(cond,"_feature.csv", sep=""), header=T, sep=","))
  table_feat <- table_feat[,-which(names(table_feat)=="cell_id")]
  table_full <- rbind(table_full, cbind(table_meta, table_feat))
}

#
# Pre process the table and create variable
#

table_full$condition[which(table_full$condition=="TSA-50nm")] <- "TSA"
table_full <- subset(table_full, condition!="DMSO")
table_full$condition <- factor(table_full$condition)

table_full <- within(table_full, group <- paste(structure, condition, cell_type))
table_full <- within(table_full, dna_volume <- (PIXEL_SIZE**3)*dna_volume)
table_full <- within(table_full, dna_io_intensity_ratio <- dna_io_intensity_outer_mean/dna_io_intensity_mid_mean)
table_full <- within(table_full, dna_io_intensity_ratio_slice <- dna_io_intensity_outer_slice_mean/dna_io_intensity_mid_slice_mean)
table_full <- within(table_full, dna_bright_spots_intensity_mean_norm <- dna_bright_spots_intensity_mean/dna_intensity_mean)

#
# Cel type-dependent size cutoff
#

table_full_cutoff <- NULL
for (ct in names(SIZE_CUTOFF)) {
  table_full_cutoff <- rbind(table_full_cutoff, subset(table_full, cell_type==ct & dna_volume>SIZE_CUTOFF[[ct]]))
}
table_full <- table_full_cutoff
rm(table_full_cutoff)

#
# Number of cells in each population
#

table(table_full$group)

write.table(x=table_full, file=paste(ROOT_FOLDER,"../engine/data-processed/880_RAW.csv", sep=""), row.names=F, sep=";")

ggplot(table_full) +
  geom_bar(aes(group, fill=condition), position="dodge") +
  facet_grid(vars(structure),vars(cell_type), scales="free", space="free") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE, drop=TRUE) +
  xlab("")

##################################################
#--------------------- PLOTS ---------------------
##################################################

#
# Size and Intensity
#

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_volume, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type)) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE) +
  xlab("")

is_different_anova(Table=table_full, num_var="dna_volume", fac_var="group")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_intensity_mean, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type), scale="free") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE) +
  xlab("")

is_different_anova(Table=table_full, num_var="dna_intensity_mean", fac_var="group")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_intensity_sum, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type), scales="free_y") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE) +
  xlab("")

is_different_anova(Table=table_full, num_var="dna_intensity_sum", fac_var="group")

#
# IO Intensity
#

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_io_intensity_ratio, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type)) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE) +
  xlab("")

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_io_intensity_ratio", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_io_intensity_ratio", fac_var="condition")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_io_intensity_ratio_slice, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type)) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLORars(structure),vars(cell_type), scales="free", space="free") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE, drop=TRUE) +
  xlab("")

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_io_intensity_ratio_slice", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_io_intensity_ratio_slice", fac_var="condition")

#
# Bright Spots
#

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_bright_spots_number, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type), scales="free_y") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE) +
  xlab("")

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_bright_spots_number", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_bright_spots_number", fac_var="condition")

ggplot(subset(table_full, cell_type=="hiPS")) +
  geom_boxplot(aes(condition, dna_bright_spots_intensity_mean, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type)) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)
  xlab("")

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_bright_spots_intensity_mean", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_bright_spots_intensity_mean", fac_var="condition")

ggplot(subset(table_full, cell_type=="hiPS")) +
  geom_boxplot(aes(condition, dna_bright_spots_intensity_mean_norm, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type)) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)
  xlab("")

is_different(Table=subset(table_full, structure="H2B"), num_var="dna_bright_spots_intensity_mean_norm", fac_var="condition")
is_different(Table=subset(table_full, structure="CBX1"), num_var="dna_bright_spots_intensity_mean_norm", fac_var="condition")

ggplot(subset(table_full, cell_type=="hiPS")) +
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
    facet_grid(vars(structure),vars(cell_type)) +
    theme_bw() +
    scale_color_brewer(palette=COLOR_PALLETE)
  
  print(fig)
}

#
# Texture distance = 1, 2
#

for (fname in feature_name) {
  
  result_table <- get_single_texture_feature_table(table_full=table_full, feature_name=fname, pixel_size=PIXEL_SIZE, distances=c(1,2))
  
  result_table$structure <- factor(result_table$structure, levels=c("H2B","CBX1"))
  
  if (CONFIG_PLOT=="set1") {
    result_table$condition <- factor(result_table$condition, levels=c("Control","DMSO","TSA"))
  } else {
    result_table$condition <- factor(result_table$condition, levels=c("Control","TSA"))
  }
  
  result_table <- subset(result_table, npixels==1 | npixels==4)
  
  fig <- ggplot(result_table, aes(interaction(distance,condition), feature, fill=condition)) +
    geom_boxplot() +
    ylab(fname) +
    facet_grid(vars(structure),vars(cell_type)) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    scale_fill_brewer(palette=COLOR_PALLETE) +
    xlab("")
  
  result_table <- within(result_table, group <- paste(structure, condition, cell_type))
  
  is_different_anova(Table=result_table, num_var="feature", fac_var="group")
  
  print(fig)
}

