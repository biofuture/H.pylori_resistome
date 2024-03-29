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
#install.packages("funr")
library(funr)
library(pheatmap)

setwd("/Users/biofuture/Documents/Desktop-contents-2021-04-28/Hpolori-Eradication/All_samples_2021_June/")

otable <- read.table(file = "./Q10_16s_subtype_B/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1)
datasetname <- "Q10_16s_subtype_B"
draw_pcoa(datasetname, otable)

Q10_species <- read.table(file = "./Q10_16s_species/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1)
datasetname <- "Q10_16s_species"
draw_pcoa(datasetname, Q10_species)
Q10_genus <- read.table(file = "./Q10_16s_genus/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1)
datasetname <- "Q10_16s_genus"
draw_pcoa(datasetname, Q10_genus)
Q10_phylum <- read.table(file = "./Q10_16s_phylum/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1)
datasetname <- "Q10_16s_phylum"
draw_pcoa(datasetname, Q10_phylum)

S14_subtype <- read.table(file = "./S14_16s_subtype_A/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1)
S14_resistome_subtype <- "S14_16s_subtype"
draw_pcoa(S14_resistome_subtype, S14_subtype)
S14_species <- read.table(file = "./S14_16s_species/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1)
S14_speciesname <- "S14_species"
draw_pcoa(S14_speciesname, S14_species)
S14_genus <- read.table(file = "./S14_16s_genus/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1)
S14_genusname <- "S14_genus"
draw_pcoa(S14_genusname, S14_genus)
S14_phylum <- read.table(file = "./S14_16s_phylum/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1)
S14_phylumname <- "S14_phylum"
draw_pcoa(S14_phylumname, S14_phylum)

##Secondlines
#BQ
BQ_16s_subtype <- read.table(file = "./BQ_16s_subtype/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1, quote = "")
draw_pcoa("BQ_16s_subtype", BQ_16s_subtype)
BQ_16s_species <- read.table(file = "./BQ_16s_species/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1, quote = "")
draw_pcoa("BQ_16s_species", BQ_16s_species)
BQ_16s_genus <- read.table(file = "./BQ_16s_genus/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1, quote = "")
draw_pcoa("BQ_16s_genus", BQ_16s_genus)
BQ_16s_phylum <- read.table(file = "./BQ_16s_phylum/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1, quote = "")
draw_pcoa("BQ_16s_phylum", BQ_16s_phylum)

#EAML
EAML_16s_subtype <- read.table(file = "./EAML_16s_subtype/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1, quote = "")
draw_pcoa("EAML_16s_subtype", EAML_16s_subtype)
EAML_16s_species <- read.table(file = "./EAML_16s_species/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1, quote = "")
draw_pcoa("EAML_16s_species", EAML_16s_species)
EAML_16s_genus <- read.table(file = "./EAML_16s_genus/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1, quote = "")
draw_pcoa("EAML_16s_genus", EAML_16s_genus)
EAML_16s_phylum <- read.table(file = "./EAML_16s_phylum/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header=TRUE, row.names = 1, quote = "")
draw_pcoa("EAML_16s_phylum", EAML_16s_phylum)


draw_pcoa <- function(datasetname, otable) {
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
cols <- c("#4DAF4A","#E41A1C", "#984EA3", "#377EB8")

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

##------------------CAP analysis-------------------------------------------------##

vare.cap <- capscale(otable ~ Sample_time, cater, dist = "bray")
anova <- anova(vare.cap, permutations = 999)
pairwise.adonis(otable, cater$Sample_time)

dms <- vare.cap$CCA$wa
dms <- as.data.frame(dms)
dms$Sample_time <- cater$Sample_time
#vare.cap$CCA$u
um <- vare.cap$CCA$u   ##centroid of samples 
um  <- as.data.frame(um)
um$Sample_time<- cater$Sample_time

CAPlabelx <- paste("CAP1 ",  ((vare.cap$CCA$eig[1] / sum(vare.cap$CCA$eig))* 10000) %/% 100 , "%", sep = "") 
CAPlabely <- paste("CAP1 ",  ((vare.cap$CCA$eig[2] / sum(vare.cap$CCA$eig))* 10000) %/% 100 , "%", sep = "") 
library(RColorBrewer)
#show_col(brewer.pal(4, "Set1"))
#cols = brewer.pal(4, "Set1")
p <- ggplot(data=dms, aes(x=CAP1,y=CAP2, color=Sample_time))
q <- p+geom_point(size=3) 
q + scale_colour_manual(values = cols) + stat_ellipse(geom = "polygon", type="euclid", alpha=0, level = 0.3)+
  scale_fill_manual(values = cols) + 
  geom_segment(aes(xend=um$CAP1, yend=um$CAP2, color=Sample_time), show.legend = F) + 
  ggtheme +  theme(legend.position = "top") +
  labs(fill="Time Points") + xlab(CAPlabelx) + ylab(CAPlabely) 

capname <- paste("CAP", datasetname, ".pdf")
ggsave(file=capname, height=8, width = 12)

#-----------------------------Above PcoA-------------------------##
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

##------------Alpha diversity resistome and microbiome---------------------##

cols <- c("#4DAF4A","#E41A1C", "#984EA3", "#377EB8")
#resistome subtype 
Q10_16s_subtype_alpha <- read.table(file = "Q10_16s_subtype_B/out/Krakenalpha/Sample_time_filter_table.txt", sep="\t", header = TRUE, row.names = 1)
S14_16s_subtype_alpha <- read.table(file = "S14_16s_subtype_A/out/Krakenalpha/Sample_time_filter_table.txt", sep="\t", header = TRUE, row.names = 1)
BQ_16s_subtype_alpha <- read.table(file = "BQ_16s_subtype/out/Krakenalpha/Sample_time_filter_table.txt", sep="\t", header = TRUE, row.names = 1)
EAML_16s_subtype_alpha <- read.table(file = "EAML_16s_subtype/out/Krakenalpha/Sample_time_filter_table.txt", sep="\t", header = TRUE, row.names = 1)
alphadiv_resistome(Q10_16s_subtype_alpha, "Q10_16s_subtype_alpha")
alphadiv_resistome(S14_16s_subtype_alpha, "S14_16s_subtype_alpha")
alphadiv_resistome(BQ_16s_subtype_alpha, "BQ_16s_subtype_alpha")
alphadiv_resistome(EAML_16s_subtype_alpha, "EAML_16s_subtype_alpha")

Q10_species_alpha <- read.table(file = "Q10_16s_species/out/Krakenalpha/Sample_time_filter_table.txt", sep="\t", header = TRUE, row.names = 1)
S14_species_alpha <- read.table(file = "S14_16s_species/out/Krakenalpha/Sample_time_filter_table.txt", sep="\t", header = TRUE, row.names = 1)
BQ_16s_species_alpha <- read.table(file = "BQ_16s_species/out/Krakenalpha/Sample_time_filter_table.txt", sep="\t", header = TRUE, row.names = 1)
EAML_16s_species_alpha <- read.table(file = "EAML_16s_species/out/Krakenalpha/Sample_time_filter_table.txt", sep="\t", header = TRUE, row.names = 1)

alphadiv(Q10_species_alpha, "Q10_16s_species_alpha")
alphadiv(S14_species_alpha, "S14_16s_species_alpha")
alphadiv(BQ_16s_species_alpha, "BQ_16s_species_alpha")
alphadiv(EAML_16s_species_alpha, "EAML_16s_species_alpha")

remove(alphad)
remove(namealpha)

alphadiv <- function(alphad, namealpha) {
  
library(tidyverse)
pv <- tidy(with(alphad, pairwise.wilcox.test(observed,Sample_time, p.adjust.method = "none")))
pv
nametodraw <- colnames(alphad)[2]
nametodraw
maxvlue <- max(alphad[,2])
minvlue <- max(alphad[,2]) - 50
ggboxplot(alphad, x = "Sample_time", y = nametodraw, notch = FALSE,main=nametodraw, fill = "Sample_time") + geom_jitter()  + 
  theme(text=element_text(size=12, family="Arial"),axis.text.x = element_text(angle = 90), legend.position = "none") +
  scale_colour_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  geom_signif(
    data=pv[c(1, 2,3, 4,5, 6),], 
    aes(xmin = group1, xmax = group2, annotations =  sprintf("p = %.2g", p.value), y_position =  c(maxvlue, maxvlue + 100, maxvlue+150, maxvlue+200, maxvlue + 250, maxvlue+300), xmin = c(1, 1, 1, 2, 3, 4), xmax = c(2, 3, 4, 3, 4,5)),
    manual= TRUE
  ) + xlab("") +  ylab("Oberved # of Species") 
namesout <- paste(namealpha, "_observed.boxplot.svg")
ggsave(file=namesout, width = 3, height = 4.5, device='svg')
pv <- tidy(with(alphad, pairwise.wilcox.test(shannon,Sample_time, p.adjust.method = "bonf")))
pv
nametodraw <- colnames(alphad)[3]
nametodraw
maxvlue <- max(alphad[,3])
minvlue <- max(alphad[,3]) - 50
ggboxplot(alphad, x = "Sample_time", y = nametodraw, notch = FALSE,main=nametodraw, fill = "Sample_time") + geom_jitter()  + 
  theme(text=element_text(size=12, family="Arial"),axis.text.x = element_text(angle = 90), legend.position = "none") +
  scale_colour_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  geom_signif(
    data=pv[c(1, 2, 3, 4, 5, 6),], 
    aes(xmin = group1, xmax = group2, annotations =  sprintf("p = %.2g", p.value),  y_position =  c(maxvlue, maxvlue + 0.2, maxvlue+0.4,maxvlue+0.6, maxvlue+0.8,  maxvlue+1), xmin = c(1, 1, 1, 2, 3, 4), xmax = c(2, 3, 4, 3,4,5)),
    manual= TRUE
  ) + xlab("") +  ylab("Shannon diversity") 
namesout <- paste(namealpha, "_shannon.boxplot.svg")
ggsave(file=namesout, width = 3, height = 4.5, device='svg')
}
?geom_signif
alphadiv_resistome <- function(alphad, namealpha) {

  library(tidyverse)
  pv <- tidy(with(alphad, pairwise.wilcox.test(observed,Sample_time, p.adjust.method = "none")))
  pv
  nametodraw <- colnames(alphad)[2]
  nametodraw
  maxvlue <- max(alphad[,2])
  minvlue <- max(alphad[,2]) - 50
  ggboxplot(alphad, x = "Sample_time", y = nametodraw, notch = FALSE,main=nametodraw, fill = "Sample_time") + geom_jitter()  + 
    theme(text=element_text(size=12, family="Arial"),axis.text.x = element_text(angle = 90), legend.position = "none") +
    scale_colour_manual(values = cols) + 
    scale_fill_manual(values = cols) +
    geom_signif(
      data=pv[c(1, 2,3, 4,5, 6),], 
      aes(xmin = group1, xmax = group2, annotations =  sprintf("p = %.2g", p.value), y_position =  c(maxvlue, maxvlue + 20, maxvlue+40, maxvlue+60, maxvlue + 80, maxvlue+100), xmin = c(1, 1, 1, 2, 3, 4), xmax = c(2, 3, 4, 3, 4,5)),
      manual= TRUE, textsize = 3.0
    ) + xlab("") +  ylab("Oberved # of ARGs subtypes") 
  namesout <- paste(namealpha, "_observed.boxplot.svg")
  ggsave(file=namesout, width = 3, height = 4.5, device='svg')
  pv <- tidy(with(alphad, pairwise.wilcox.test(shannon,Sample_time, p.adjust.method = "none")))
  pv
  nametodraw <- colnames(alphad)[3]
  nametodraw
  maxvlue <- max(alphad[,3])
  minvlue <- max(alphad[,3]) - 50
  ggboxplot(alphad, x = "Sample_time", y = nametodraw, notch = FALSE,main=nametodraw, fill = "Sample_time") + geom_jitter()  + 
    theme(text=element_text(size=12, family="Arial"),axis.text.x = element_text(angle = 90), legend.position = "none") +
    scale_colour_manual(values = cols) + 
    scale_fill_manual(values = cols) +
    geom_signif(
      data=pv[c(1, 2, 3, 4, 5, 6),], 
      aes(xmin = group1, xmax = group2, annotations =  sprintf("p = %.2g", p.value),  y_position =  c(maxvlue, maxvlue + 0.2, maxvlue+0.4,maxvlue+0.6, maxvlue+0.8,  maxvlue+1), xmin = c(1, 1, 1, 2, 3, 4), xmax = c(2, 3, 4, 3,4,5)),
      manual= TRUE,textsize = 3.0
    ) + xlab("") +  ylab("Shannon diversity") 
  namesout <- paste(namealpha, "_shannon.boxplot.svg")
  ggsave(file=namesout, width = 3, height = 4.5, device='svg')
}

##Mechanisms of resistance in each time points 
Q10_16s_Mechanism <- read.table(file = "Q10_16s_Mechanism/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header = TRUE, row.names = 1)
S14_16s_Mechanism <- read.table(file = "S14_16s_Mechanism/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header = TRUE, row.names = 1)
BQ_16s_Mechanism <- read.table(file = "BQ_16s_Mechanism/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header = TRUE, row.names = 1)
EAML_16s_Mechanism <- read.table(file = "EAML_16s_Mechanism/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header = TRUE, row.names = 1)

draw_boxplot_all_columns(Q10_16s_Mechanism, "Q10_16s_Mechanism", "Q10_16s_Mechanism_boxplot")
draw_boxplot_all_columns(S14_16s_Mechanism, "S14_16s_Mechanism", "S14_16s_Mechanism_boxplot")
draw_boxplot_all_columns(BQ_16s_Mechanism, "BQ_16s_Mechanism", "BQ_16s_Mechanism_boxplot")
draw_boxplot_all_columns(EAML_16s_Mechanism, "EAML_16s_Mechanism", "EAML_16s_Mechanism_boxplot")

remove(alphad)

draw_boxplot_all_columns <- function(alphad, nametofile, dirnames){
  
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
   xlab("") +  ylab("# of ARGs per 16S rDNA") 
  namesout <- paste0(dirnames, "/", nametofile, "_", nametodraw, ".boxplot.svg")
  ggsave(file=namesout, width = 3, height = 4.5, device='svg')
  }
}

?stat_compare_means
##Type level resistome boxplot plotting for each treatments 
Q10_16s_type_boxplot <- read.table(file = "Q10_16s_type/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header = TRUE, row.names = 1)
S14_16s_type_boxplot <- read.table(file = "S14_16s_type/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header = TRUE, row.names = 1)
BQ_16s_type_boxplot <- read.table(file = "BQ_16s_type/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header = TRUE, row.names = 1)
EAML_16s_type_boxplot <- read.table(file = "EAML_16s_type/out/Kraken/TaxProfileWithMeta.txt", sep="\t", header = TRUE, row.names = 1)

draw_boxplot_all_columns(Q10_16s_type_boxplot, "Q10_16s_type", "Q10_16s_type_boxplot")
draw_boxplot_all_columns(S14_16s_type_boxplot, "S14_16s_type", "S14_16s_type_boxplot")
draw_boxplot_all_columns(BQ_16s_type_boxplot, "BQ_16s_type", "BQ_16s_type_boxplot")
draw_boxplot_all_columns(EAML_16s_type_boxplot, "EAML_16s_type", "EAML_16s_type_boxplot")

##-----------check the types of resistance and their respond to treatment--------------------##

##S14 ammoxlin used beta lactam type  clarimycin is macrolide antibiotics 
##Q10 tetracycline is tetracycline type 
#EAML used ammoxlin and levofloxacin which are beta-lactam and quinolone class 
#BQ used tetracycline and bismutch 

#When the background resistance type/class is high before treatment, it is more likely to amplify the resistance burden while use antibiotics 
#to treat diseases, this is because the resistance bugs could take over the niche of other suspective bacteria that were killed during the 
#treatment, hence increase the level of resistance burden after treatment. If background resistance is low for that type, it is not easily 
#accumulated the resistance burden for that type, hence will be better treated for the diseases, however with bigger influence on microbial community

#heatmap of all differntial abundant subtypes 

##differntial abundant antimicrobial resistance types/subtypes with time 
Q10_16s_subtype_heatmap <- read.table(file="./Q10_16s_subtype_B/out/Krakendiff-abundance/Sample_time_filter_table.txt", sep="\t",row.names = 1, header = TRUE, quote = "")
Q10_16s_subtype_sig <- read.table(file="./Q10_16s_subtype_B/out/Krakendiff-abundance/Sample_time_kw_stat_summary.txt", header=TRUE, sep="\t")
S14_16s_subtype_heatmap <- read.table(file="./S14_16s_subtype_A/out/Krakendiff-abundance/Sample_time_filter_table.txt", sep="\t",row.names = 1, header = TRUE, quote = "")
S14_16s_subtype_sig <- read.table(file="./S14_16s_subtype_A/out/Krakendiff-abundance/Sample_time_kw_stat_summary.txt", header=TRUE, sep="\t")
BQ_16s_subtype_heatmap <- read.table(file="./BQ_16s_subtype//out/Krakendiff-abundance/Sample_time_filter_table.txt", sep="\t",row.names = 1, header = TRUE, quote = "")
BQ_16s_subtype_sig <- read.table(file="./BQ_16s_subtype//out/Krakendiff-abundance/Sample_time_kw_stat_summary.txt", header=TRUE, sep="\t")
EAML_16s_subtype_heatmap <- read.table(file="./EAML_16s_subtype/out/Krakendiff-abundance/Sample_time_filter_table.txt", sep="\t",row.names = 1, header = TRUE, quote = "")
EAML_16s_subtype_sig <- read.table(file="./EAML_16s_subtype//out/Krakendiff-abundance/Sample_time_kw_stat_summary.txt", header=TRUE, sep="\t")

draw_heatmap(Q10_16s_subtype_heatmap, Q10_16s_subtype_sig, "Q10_16s_subtype_sig")
draw_heatmap(S14_16s_subtype_heatmap, S14_16s_subtype_sig, "S14_16s_subtype_sig")
draw_heatmap(BQ_16s_subtype_heatmap, BQ_16s_subtype_sig, "BQ_16s_subtype_sig")
draw_heatmap(EAML_16s_subtype_heatmap, EAML_16s_subtype_sig, "EAML_16s_subtype_sig")

draw_heatmap <- function(q10subtype, sigq10, nametofile){
  
sig <- sigq10 %>% filter(kw_pval < 0.05 & FDR < 0.05)
q10subtype <- q10subtype %>% arrange(Sample_time)
cater <- q10subtype %>% select(Sample_time)
q10subtype <- q10subtype %>% select(!Sample_time)
sigsubtype <- t(q10subtype)[c(sig$Row.names),]

trans <- sigsubtype
trans <- trans[rowSums(trans==0, na.rm=TRUE)<ncol(q10subtype), ] 
trans <- trans[ order(row.names(trans)), ]

namesout <- paste0(nametofile, ".heatmap.pdf")
pdf(namesout, width = 18, height = 24)
pheatmap(log10(trans + 0.000001), cluster_rows = FALSE,cluster_cols = FALSE, border_color=NA, annotation_col = cater,
         scale = "row", show_colnames = FALSE, fontsize_row = 9,  
         #cutree_cols = 4, 
         color = colorRampPalette(c("steelblue","snow","tomato3"))(1000))
dev.off()
}



