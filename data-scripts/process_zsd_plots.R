rm(list=ls())
library(tidyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(gridExtra)

PIXEL_SIZE <- 0.108

CONFIG_PLOT <- "set2"

SIZE_CUTOFF = list(hiPS=230000*(PIXEL_SIZE**3),
                   cardiomyocyte=0*(PIXEL_SIZE**3))

COLOR_PALLETE <- "Dark2" #Dark2

ROOT_FOLDER <- "/home/matheus.viana/projects/qcb/qcb-analysis/data-raw/"

setwd(ROOT_FOLDER)

#
# Functions
#

get_single_texture_feature_table <- function(table_full, feature_name, pixel_size, first_distance_only=FALSE) {
  
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

    if (first_distance_only) {

      subTable_long <- subset(subTable_long, npixels==npixels[1])
      
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

is_different_ttest <- function(Table, num_var,fac_var) {
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

is_different_anova <- function(Table, num_var,fac_var) {
  aov_table <- aov(Table[,which(names(Table)==num_var)]~Table[,which(names(Table)==fac_var)])
  print(summary(aov_table))
  result <- TukeyHSD(aov_table, p.adjust.method ="none")$`Table[, which(names(Table) == fac_var)]`
  print(result)
}

#
# Loading data
#

conditions <- c("NUCLEUS_HIPS_H2B-TSA", "NUCLEUS_HIPS_H2B+TSA", "NUCLEUS_HIPS_CBX1-TSA", "NUCLEUS_HIPS_CBX1+TSA", "NUCLEUS_CARDIO")

table_full <- NULL
for (cond in conditions) {
  table_meta <- data.frame(read.table(paste(cond,"_meta.csv", sep=""), header=T, sep=","))
  table_feat <- data.frame(read.table(paste(cond,"_feature.csv", sep=""), header=T, sep=","))
  table_feat <- table_feat[,-which(names(table_feat)=="cell_id")]
  table_full <- rbind(table_full, cbind(table_meta, table_feat))
}

#
# Setting specific configuration set
#

config <- list(
  set1=
    list(
      hiPS=list(
        list("Control","20181005_R02_002"),
        list("DMSO","20181005_R02_001"),
        list("TSA","20181005_R01_001"),
        list("TSA","20181005_R01_002"),
        list("TSA","20181005_R01_003"),
        list("Control","20181016_R02"),
        list("DMSO","20181016_R03"),
        list("TSA","20181016_R01")
      ),
      cardiomyocyte=list(
        list("DMSO","20181114_R010"),
        list("TSA","20181114_R08"),
        list("Control","20181114_R09"),
        list("TSA","20181114_R05"),
        list("Control","20181114_R06"),
        list("DMSO","20181114_R07")
      )
    ),
  set2=
    list(
      hiPS=list(
        list("Control","20181005_R02_002"),
        list("Control","20181012_R02"),
        list("Control","20181026_R02"),
        list("TSA","20181005_R01_001"),
        list("TSA","20181005_R01_002"),
        list("TSA","20181005_R01_003"),
        list("TSA","20181012_R01"),
        list("TSA","20181026_R01"),
        list("Control","20181016_R02"),
        list("Control","20181019_R04"),
        list("TSA","20181016_R01"),
        list("TSA","20181019_R03")
      ),
      cardiomyocyte=list(
        list("TSA","20181107_R04"),
        list("TSA","20181114_R08"),
        list("Control","20181114_R09"),
        list("Control","20181107_R06"),
        list("TSA","20181107_R01"),
        list("Control","20181107_R03"),
        list("TSA","20181114_R05"),
        list("Control","20181114_R06")
      ))
)

table_full$valid <- FALSE
table_full$condition <- as.character(table_full$condition)
for (ct in unique(table_full$cell_type)) {
  for (exp_id in unique(table_full$czi)) {
    exp_id_char = as.character(exp_id)
    for (n in seq(1,length(config[[CONFIG_PLOT]][[ct]]),1)) {
      if (grepl(config[[CONFIG_PLOT]][[ct]][[n]][2], exp_id_char) == TRUE) {
        valids <- which(table_full$czi==exp_id)
        table_full$valid[valids] <- TRUE
        table_full$condition[valids] <- unlist(config[[CONFIG_PLOT]][[ct]][[n]][1])
      }
    }
  }
}
table_full <- subset(table_full, valid==TRUE & dna_intensity_mean>400)
table_full$condition <- as.factor(table_full$condition)

#
# Pre process the table and create variable
#

table_full <- within(table_full, dna_intensity_mean<-(dna_intensity_mean-400))
table_full <- within(table_full, dna_intensity_sum<-(dna_intensity_sum-(400*dna_volume)))
table_full <- within(table_full, dna_bright_spots_intensity_mean<-dna_bright_spots_intensity_mean-400)

table_full <- within(table_full, group <- paste(structure, condition, cell_type))
table_full <- within(table_full, dna_volume <- (PIXEL_SIZE**3)*dna_volume)
table_full <- within(table_full, dna_io_intensity_ratio <- dna_io_intensity_outer_mean/dna_io_intensity_mid_mean)
table_full <- within(table_full, dna_io_intensity_ratio_slice <- (dna_io_intensity_outer_slice_mean-400)/(dna_io_intensity_mid_slice_mean-400))
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

write.table(x=table_full, file=paste(ROOT_FOLDER,"../engine/data-processed/ZSD_",CONFIG_PLOT,".csv", sep=""), row.names=F, sep=";")

ggplot(table_full) +
  geom_bar(aes(condition, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type), scales="free_y") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE) +
  xlab("")

##################################################
#--------------------- PLOTS ---------------------
##################################################

#
# Size and Intensity

#

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_volume, fill=condition)) +
  #stat_summary(aes(condition, dna_volume, fill=condition), fun.y = mean, geom="point") +
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

is_different_anova(Table=table_full, num_var="dna_io_intensity_ratio", fac_var="group")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_io_intensity_ratio_slice, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type)) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE) +
  xlab("") + ylim(c(0.6,1.1))

is_different_anova(Table=table_full, num_var="dna_io_intensity_ratio_slice", fac_var="group")

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

is_different_anova(Table=table_full, num_var="dna_bright_spots_number", fac_var="group")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_bright_spots_intensity_mean, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type), scales="free_y") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE)
  xlab("")

is_different_anova(Table=table_full, num_var="dna_bright_spots_intensity_mean", fac_var="group")

ggplot(table_full) +
  geom_boxplot(aes(condition, dna_bright_spots_intensity_mean_norm, fill=condition)) +
  facet_grid(vars(structure),vars(cell_type), scales="free_y") +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1)) +
  scale_fill_brewer(palette=COLOR_PALLETE) +
  xlab("") + ylim(1,5)

is_different_anova(Table=subset(table_full,dna_bright_spots_intensity_mean_norm<5), num_var="dna_bright_spots_intensity_mean_norm", fac_var="group")

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
# Texture distance = 1
#

for (fname in feature_name) {
  
  result_table <- get_single_texture_feature_table(table_full=table_full, feature_name=fname, pixel_size=PIXEL_SIZE, first_distance_only=TRUE)
  
  if (CONFIG_PLOT=="set1") {
    result_table$condition <- factor(result_table$condition, levels=c("Control","DMSO","TSA"))
  } else {
    result_table$condition <- factor(result_table$condition, levels=c("Control","TSA"))
  }
  
  fig <- ggplot(result_table, aes(condition, feature, fill=condition)) +
    geom_boxplot() +
    ylab(fname) +
    facet_grid(vars(structure),vars(cell_type)) +
    theme_bw() +
    scale_fill_brewer(palette=COLOR_PALLETE) +
    xlab("distance 1px (0.108um)")
  
  result_table <- within(result_table, group <- paste(structure, condition, cell_type))
  
  is_different_anova(Table=result_table, num_var="feature", fac_var="group")
  
  print(fig)
}

