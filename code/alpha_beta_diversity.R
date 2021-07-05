#Record everything and make them transparency 
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(broom)
library(svglite)
library(extrafont)
library(vegan)
library(labdsv)
library(scales)

setwd("/Users/biofuture/Documents/Desktop-contents-2021-04-28/Hpolori-Eradication/All_samples_2021_June/")
otable <- read.table(file = "S14_Q10_16s_subtype.output_table.txt", sep="\t")
tb <- t(otable)
colnames(tb) <- otable[,1]
tb[-1,] 
cater <- read.table(file = "meta_data_firstline.txt",sep="\t", header=TRUE)
inner_join(otable, cater, by = "Name") 

vd <- vegdist(otable,method="jaccard")
vd.pco <- pco(vd, k=4)
ccav <- cca(otable)
 
ccav
library(pairwiseAdonis)
adonis2(otable ~ Sample_time, data = cater, permutations = 999)
