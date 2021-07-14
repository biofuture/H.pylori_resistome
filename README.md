# H.pylori_resistome
This is the R code for H.pylori eradication resistome

This code is specially for longitudinal microbiome data graphics plotting

##alpha diversity in the codes

After geting the big table, subsample to run analysis for each Trail

Boxplot showing temporal dynamics of alpha diversity

```
library(ggplot2)
perl /srv/scratch/mrcbio/scripts/sub_sample_with_metadata.pl merge_16s_subtype.all.update.txt meta_data_firstline.txt Firstline S14_16s_subtype Group:S14

#go to conda env with appropriate R dependent software
conda activate shotgunpipe
python /srv/scratch/mrcbio/bin/visulization/visualisation.py --data_matrix S14_16s_subtype.output_table.txt --metadata ../meta_data_firstline.txt --treatment_groups Sample_time --out_dir S14_16s_subt


```



