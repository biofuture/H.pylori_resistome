##Sharing examples for SCFA 
##Developed by Dr Xiaotao Jiang 2021-07-21

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
library(funr)
library(pheatmap)
library(pairwiseAdonis)

#install packages with the following commands 
#install.packages("funr")
#colours for your groups 
cols <- c("#4DAF4A","#E41A1C", "#984EA3", "#377EB8")

otable <- read.table(file = "tables/scfa_uc_test.txt", sep="\t", header=TRUE, row.names = 1)
datasetname <- "scfa"
draw_pcoa(datasetname, otable, cols)
draw_boxplot_all_columns(otable, "scfa", "uc-example", cols)

draw_boxplot_all_columns <- function(alphad, nametofile, dirnames, cols){
  library(tibble)
  #library(funr)
  #script.dir <- funr::get_script_path()
  # 
  #ifelse(!dir.exists(file.path(script.dir,dirnames)), dir.create(file.path(script.dir, dirnames)), FALSE)
  ifelse(!dir.exists(file.path(dirnames)), dir.create(file.path(dirnames)), FALSE)
  alphad$total_resistome <- rowSums(alphad[, 2:ncol(alphad)])
  for (i in 2:(dim(alphad)[2])) {
    #alphad[,i] <- log2(alphad[,i] + 0.00000001) 
    nametodraw <- colnames(alphad)[i]
    nametodraw
    pv <- tidy(with(alphad, pairwise.wilcox.test(alphad[,i],Sample_time, p.adjust.method = "none")))
    pv
    maxvlue <- max(alphad[,i])
    minvlue <- max(alphad[,i]) - 0.2
    comp <-  combn(levels(as.factor(alphad[,c("Sample_time")])),2,list)
    
    ggboxplot(alphad, x = "Sample_time", y = nametodraw, notch = FALSE,main=nametodraw, fill = "Sample_time") + geom_jitter()  + 
      theme(text=element_text(size=12, family="Arial"),axis.text.x = element_text(angle = 90), legend.position = "none") +
      scale_colour_manual(values = cols) + 
      scale_fill_manual(values = cols) +
      #geom_signif(
      # data=pv[c(1, 2,3, 4,5, 6),], 
      #aes(xmin = group1, xmax = group2, annotations =  sprintf("p = %.2g", p.value), y_position =  c(maxvlue, maxvlue + 0.2, maxvlue+0.4, maxvlue+0.6, maxvlue + 0.8, maxvlue+1), xmin = c(1, 1, 1, 2, 3, 4), xmax = c(2, 3, 4, 3, 4,5)),
      #manual= TRUE
      #) +
      stat_compare_means(comparisons = comp,aes(label = paste0("p = ", ..p.format..)), method = "wilcox.test",  p.adjust.method = "bonf", hide.ns =TRUE) +
      xlab("") +  ylab("SCFA abundance") 
    namesout <- paste0(dirnames, "/", nametofile, "_", nametodraw, ".boxplot.pdf")
    ggsave(file=namesout, width = 3, height = 4.5, device='pdf')
  }
}

draw_pcoa <- function(datasetname, otable, cols) {
  ###--------------------PcoA--------------------------------------##
  cater <- ""
  cater$Sample_time  <- otable$Sample_time
  otable <- subset(otable, select = -Sample_time)
  #vd <- vegdist(otable,method="jaccard")
  #vd.pco <- pco(vd, k=4)
  #ccav <- cca(otable)
  #library(pairwiseAdonis)
  #adonis2(otable ~ Sample_time, data = cater, permutations = 999)
  #pairwise.adonis(vd, cater$Sample_time)
  
  ##beta dispersal 
  #vare.cap <- betadisper(vd, group = cater$Sample_time)
  #anova(vare.cap, permutations = 999)
  #plot(mod, hull = FALSE, ellipse = TRUE)
  ##PCOA with prcomp 
  #vd <- vegdist(otable,method="bray")
  vd <- vegdist(otable,method="bray")
  #vd.pco <- pco(vd, k=4)
  
  ##refer to scripts https://github.com/Tom-Jenkins/utility_scripts/blob/master/visualise_dapc_pca.
  df.pca <- prcomp(vd, scale. =T)
  pca.results <- data.frame(df.pca$x)[, 1:2]
  #variance explained
  #eigs <- df.pca$sdev^2
  #eigs[1] / sum(eigs)
  sumpca <- summary(df.pca)
  PC1label <- paste("PC1 ", sumpca$importance[,1][2] * 100, "%", sep = "")  
  PC2label <- paste("PC2 ", sumpca$importance[,2][2] * 100, "%", sep = "")
  
  pca.results$target <- cater$Sample_time
  centroid <- aggregate(cbind(PC1, PC2) ~ target, data = pca.results,  FUN = mean)
  
  library(RColorBrewer)
  #show_col(brewer.pal(4, "Set1"))
  #cols1 = brewer.pal(12, "Set3")
  #cols2 = brewer.pal(8,"Set2")
  #?brewer.all
  #cols = c(cols, cols2)
  #cols = brewer.pal(4, "Set1")
  #cols <- c("#4DAF4A","#E41A1C", "#984EA3", "#377EB8")
  
  pca.results1 <- dplyr::left_join(pca.results, centroid, by = "target", suffix = c("", ".cen"))
  
  ggplot(data=pca.results1, aes(PC1, PC2, color=target)) + 
    #geom_hline(yintercept = 0) +
    #geom_vline(xintercept = 0) +
    geom_point(size = 3)+ stat_ellipse(geom = "polygon", type="euclid", alpha=0, level = 0.3) + 
    geom_segment(aes(xend=PC1.cen, yend=PC2.cen, color=target), show.legend = T) + 
    scale_fill_manual(values = cols)+
    scale_colour_manual(values = cols)+
    ggtheme +  theme(legend.position = "top") + 
    geom_label(data = centroid, 
               aes(label = target, fill = target), size = 5, show.legend = F, color = "white") +
    xlab(PC1label) + ylab(PC2label)  
  pcaname <- paste("PcoA", datasetname, ".pdf")
  ggsave(file=pcaname, height=8, width = 12)

  #-----------------------------Above PcoA-------------------------##
  
  adonisout <- pairwise.adonis(otable, cater$Sample_time)
  adonisname <- paste("Adonis-test-perm999", datasetname, ".csv")
  write.csv(file=adonisname, adonisout)
  
}

# Custom ggplot2 theme
ggtheme = theme(legend.title = element_blank(),
                axis.text.y = element_text(colour="black", size=14),
                axis.text.x = element_text(colour="black", size=14),
                axis.title = element_text(colour="black", size=14),
                legend.position = "right",
                legend.text = element_text(size=15),
                legend.key = element_rect(fill = NA),
                legend.key.size = unit(0.7, "cm"),
                legend.box.spacing = unit(0, "cm"),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                # title centered
                plot.title = element_text(hjust=0.5, size=25) 
)
