# Bulk-RNA-seq-analysis-on-Gout-data
A short bulk RNA-seq analysis performed in R based on a class-provided gout dataset. Code is adapted and expanded to facilitate a basic RNA-seq workflow starting from DEseq2 data, including analysis such as volcano plots, PCA, immune cell deconvolution and heatmap creation.

Hiya, this is just a relatively straightforward project that aims to perform some basic, but useful analysis on a pre-existing gout transcript dataset. Although this is not my first project, it is my first repo made in github so forgive me for any transgressions I may commit.


# Introduction
This analysis will be predominantly performed in R version 4.4.1 and will use these initial data files: gout_sample_sheet, gout_em, gout_DE, and gout_annotations (see in the dataset folder "gout").
This will hopefully walk through the steps required for the analysis of bulk RNA-sequencing data, including pre-processing and data visualisation (note that differential gene analysis by DESeq2 has already been conducted on the supplied data).

## Packages used
```
install.packages("reshape")  
install.packages("extrafont")  
install.packages("cowplot")  
install.packages("amap")  
install.packages("pubr")  
install.packages("ggsignif")  
install.packages("pheatmap")  
install.packages("BiocManager")  
BiocManager::install("clusterProfiler")  
BiocManager::install("org.Hs.eg.db")  
BiocManager::install("STRINGdb")  

library(ggplot2)                  
library(ggrepel)                  
library(extrafont)                                
library(cowplot)                 
library(reshape2)                 
library(ggsignif)                
library(amap)                    
library(pheatmap)                
library(clusterProfiler)          
library(org.Hs.eg.db)             
library(STRINGdb)                
```
================================================================================
## Font import
This is an optional step that imports additional fonts onto R, it's not really required but does make the text on the graphs look better (personally).
```
font_import()      
loadfonts(font_import(paths = "C:/Windows/Fonts", prompt = FALSE))  
windowsFonts()                     #Just to check if it has worked (it should)
```
## Theme creation
Obviously you can use any desired theme, but I decided to create my own that is modified off "cowplot".  
```
theme_supreme <- function()  {
  font <- "Arial Nova"             
  theme_cowplot(12) %+replace%      
  theme(
    panel.grid.major = element_blank(),                           #strip the gridlines
    panel.grid.minor = element_blank(),                           
    panel.background = element_rect(),
    panel.border = element_rect(color="black", size=1),           # Have a nice border around the plots
    strip.background = element_blank(),
    strip.text = element_text(color = "black"), 
    axis.text = element_text(family=font),
    axis.title.x = element_text(family=font, 
                                face="bold"),                     #Want text bold
    axis.title.y = element_text(family=font, 
                                angle=90, 
                                margin = margin(r=5),            #the y-axis usually looks a bit cramped, so increase margin
                                face="bold")
  )
}
```

# Data pre-processing
#Required:
GOUT_em, GOUT_de, GOUT_annotations, GOUT_ss
```
##Merge steps
names(GOUT_annotations) = c("Gene", "Chromosome", "Start", "Stop")        #Change labels to something more useful
master_temp = merge(GOUT_em, GOUT_annotations, by.x = 0, by.y = 0)
master = merge(master_temp, GOUT_de, by.x = 1, by.y = 0)                  

#Cleaning up
row.names(master) = master[,"Gene"]                                      #Replace ENSEMBLE with gene name
master = master[,-c(1,30)]                                               #Remove redundant columns
master = na.omit(master)                                                 #Remove N/A values
sorted_order = order(master[,"p.adj"], decreasing=FALSE)
master = master[sorted_order,]                                           #Sort it out
```

Now to modify the table to include important information (like minus log10p or the directionality of the change)
```
#Create new columns (mlop10p, sig, direction)
master$mean_expression = rowMeans(master[,2:28])
master$mlog10p = -log10(master$p.adj)                                     #This is for mlog10p
master$sig = as.factor(master$p.adj<0.05 & abs(master$log2fold) > 1.0)    #This is for sig

master_sig_up = subset(master, p.adj<0.05 & log2fold>1)                   #Sig genes that are upregulated
master_sig_down = subset(master, p.adj<0.05 & log2fold< -1)               #Sig genes that are downregulated
master_non_sig = subset(master, sig==FALSE)                               #Genes that are not differentially expressed

master_non_sig$direction = "No change"
master_sig_down$direction = "Down"
master_sig_up$direction = "Up"
master=rbind(master_non_sig, master_sig_down, master_sig_up)              #Adding directionality (up, down, non)
sorted_order = order(master[,"p.adj"], decreasing=FALSE)
master = master[sorted_order,]

#Labelling the "top" genes (down, up, and top 6)
master_sig_up_top3 = master_sig_up[1:3,]  
master_sig_down_top3 = master_sig_down[1:3,]
master_top6=rbind(master_sig_up_top3, master_sig_down_top3)               #Creating a top 6 list
```

Putting it all together (making more useful tables (em_scaled, sig, etc))
```
GOUT_em_symbols = master[,GOUT_ss$sample] 
GOUT_em_scaled = na.omit(data.frame(t(scale(t(GOUT_em_symbols)))))

master_sig = subset(master, p.adj<0.05 & abs(master$log2) > 1.0)
GOUT_sig_genes = c(rownames(master_sig,))

GOUT_em_symbols_sig = GOUT_em_symbols[GOUT_sig_genes,]
GOUT_em_scaled_sig = GOUT_em_scaled[GOUT_sig_genes,]
```

# Data visualisation
## Expression density plot
Required
GOUT_em, ggplot, reshape2, and theme_supreme()

Making the required melted table
```
GOUT_em.m = melt(GOUT_em)
GOUT_em.m = na.omit(GOUT_em.m)
GOUT_em.m = merge(GOUT_em.m, GOUT_ss, by.x = "variable", by.y = "sample")
```
#Generating the expression density plot/facet
```
exp_dens = ggplot(GOUT_em.m, aes(x=log10(value+0.01), fill=sample_group)) +     #Add the "value + 0.01" to include smaller values
  geom_density() +
  scale_fill_manual("Sample Group", breaks = c("HC", "GOUT"),
                    values =c("#1B9E77", "#F4A261"), 
                    labels=c("HC", "GOUT")) +
  facet_wrap(~variable, ncol=7) +                                               #Dealing with 28 samples so rows of 7
  theme_supreme()+ theme(strip.text = element_text(size=5), 
                         axis.text = element_text(size=8)) +
  labs(x="Expression (log10)", y="Density") +
  xlim(c(-2,5))

exp_dens
```
And it should look like this. Note the bimodal distribution of peaks, which depicts a fairly good quality of transcript reads. 
<img width="789" height="584" alt="image" src="https://github.com/user-attachments/assets/0755a14f-6d45-41bc-945f-41825e98d4fb" />
