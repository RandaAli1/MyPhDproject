# Analysing and displaying the enrichment of RTEs in the euchromatin domain of Roadmap cell types. 


# This code calculates the emprical P-value for the enrichment of RTEs in the euchromatin domains of the 127 cell types analysed 
# by the Roadmap project (Roadmap Epigenomics Consortium, 2015). Next, the number of cell types enriched with one or more RTEs are counted by tissue group.
# Finally, a bar chart representing counts of cell types enriched with RTEs per tissue group is plotted.  

Mydir = c("C:/Users/rali1/OneDrive - University of Plymouth/Thesis_2019/MyAnalysis/")

setwd(paste0(Mydir))

# Load required libraries:

library(readxl)
library(dplyr)
library(plyr)
library(ggplot2)
## Load files
Data = read.delim("Enrichment_in_Euchromatin_regions/Results/TotalRTEsInEuchromatinRegions.txt", header = T)
Eidnames = read.delim("Enrichment_in_Euchromatin_regions/Results/RoadmapCellNames_0502020.txt", header = T, sep = "\t")

Random_L1 = read.delim("Enrichment_in_Euchromatin_regions/Results/RandomL1Hs_1k_H3K4Euchromatin.txt", header = T, sep = "\t")
Random_Alu = read.delim("Enrichment_in_Euchromatin_regions/Results/RandomAluY_1k_H3K4Euchromatin.txt", header = T, sep = "\t")
Random_SVA = read.delim("Enrichment_in_Euchromatin_regions/Results/RandomSVA_1k_H3K4Euchromatin.txt", header = T, sep = "\t")

## Create loop for calculating the empirical p-value for the enrichment of L1s per cell type:

Sheet1 = data.frame("Total" = Data$L1, "EID" = Data$EID)

Results = NA

for (row in 1:nrow(Sheet1)) {
  rep = 1000
  n = as.data.frame(Random_L1[row,])
  n$EID = NULL
  v = data.frame(as.numeric(t(n)))
  bar <- subset(v, v$as.numeric.t.n.. >= Sheet1$Total[row])
  obs = nrow(bar)
  empPval = (obs+1)/(rep+1)
  df = data.frame("EID"= Sheet1$EID[row], "P-val"=empPval)
  Results <<- rbind(df, Results)
}

# Extract cell types with a significant enrichment of L1s based on P-value <= 0.05:

sign <- subset(Results, Results$P.val <= 0.05)
L1_sign = data.frame(Type = "L1", sign)
# Match EID with cell names:

L1names = inner_join(sign, Eidnames, by = "EID")

# Count the number of cell types enriched in L1s per anatomical region:

df_L1 = data.frame(aggregate(cbind(count = Eid_id) ~ ANATOMY, 
                             data = L1names, 
                             FUN = function(x){NROW(x)}))

### Repeat above steps for Alus and SVAs ##

## Create loop for calculating the empirical p-value for the enrichment of Alus per cell type:

Sheet1 = data.frame("Total" = Data$Alu, "EID" = Data$EID)
Results = NA

for (row in 1:nrow(Sheet1)) {
  rep = 1000
  n = as.data.frame(Random_Alu[row,])
  n$EID = NULL
  v = data.frame(as.numeric(t(n)))
  bar <- subset(v, v$as.numeric.t.n.. >= Sheet1$Total[row])
  obs = nrow(bar)
  empPval = (obs+1)/(rep+1)
  df = data.frame("EID"= Sheet1$EID[row], "P-val"=empPval)
  Results <<- rbind(df, Results)
}
1/1001
# Extract cell types with a significant enrichment of Alus based on P-value <= 0.05:

sign <- subset(Results, Results$P.val <= 0.05)
Alu_sign = data.frame(Type = "Alu", sign)
# Match EID with cell names:

Alunames = inner_join(sign, Eidnames, by = "EID")

# Count the number of cell types enriched in Alus per anatomical region:

df_Alu = data.frame(aggregate(cbind(count = Eid_id) ~ ANATOMY, 
                             data = Alunames, 
                             FUN = function(x){NROW(x)}))

## Create loop for calculating the empirical p-value for the enrichment of SVAs per cell type:

Sheet1 = data.frame("Total" = Data$SVA, "EID" = Data$EID)
Results = NA

for (row in 1:nrow(Sheet1)) {
  rep = 1000
  n = as.data.frame(Random_SVA[row,])
  n$EID = NULL
  v = data.frame(as.numeric(t(n)))
  bar <- subset(v, v$as.numeric.t.n.. >= Sheet1$Total[row])
  obs = nrow(bar)
  empPval = (obs+1)/(rep+1)
  df = data.frame("EID"= Sheet1$EID[row], "P-val"=empPval)
  Results <<- rbind(df, Results)
}
# Extract cell types with a significant enrichment of SVAs based on P-value <= 0.05:

sign <- subset(Results, Results$P.val <= 0.05)
SVA_sign = data.frame(Type = "SVA", sign)
# Match EID with cell names:

SVAnames = inner_join(sign, Eidnames, by = "EID")

# Count the number of cell types enriched in SVAs per anatomical region:

df_SVA = data.frame(aggregate(cbind(count = Eid_id) ~ ANATOMY, 
                             data = SVAnames, 
                             FUN = function(x){NROW(x)}))


# Group all P-values into one data frame:

Pval_SVA_Alu = merge(SVA_sign, Alu_sign, by = "EID" , all.x = TRUE)
Pval_all = merge(Pval_SVA_Alu,L1_sign, by ="EID" , all.x = TRUE )
colnames(Pval_all) = c("EID", "SVA", "Pval_SVA", "Alu", "Pval_Alu","L1","Pval_L1")

## Write data frame of P-values to file:

#write.table(Pval_all, file = "Enrichment_in_Euchromatin_regions/Results/RTEEnrichmentInEuchromatin_PvalByEID.txt", sep="\t", row.names=F,col.names=T, quote = F)


## Group all counts into one data frame:

merged_SVA_Alu = merge(df_SVA, df_Alu, by = "ANATOMY" , all.x = TRUE)
merged_all = merge(merged_SVA_Alu, df_L1, by = "ANATOMY" , all.x = TRUE)
colnames(merged_all) = c("Anatomy", "SVA", "Alu", "L1")
merged_all[is.na(merged_all)] <- 0
all_counts = data.frame("Anatomy" = merged_all$Anatomy, "L1" = merged_all$L1, "Alu" = merged_all$Alu, "SVA" = merged_all$SVA)

## Write data frame to file:

#write.table(all_counts, file = "Enrichment_in_Euchromatin_regions/Results/RTEEnrichmentInEuchromatin_CountsByAnatomyGroup.txt", sep="\t", row.names=F,col.names=T, quote = F)


##### Plot counts of RTEs enriched in euchromatin regions of various tissue types:

# Create dataframe to define X and Y datapoints for ggplot:

df_L1$Type =  c("L1")
df_Alu$Type =  c("Alu")
df_SVA$Type = c("SVA")

L1_Alu = full_join(df_L1,df_Alu)
All_RTE = full_join(L1_Alu, df_SVA)

# Create the plot: 

p = ggplot(data=All_RTE, aes(x=ANATOMY, y=count, fill=Type)) +
  # To add a black line around the bars and change the width
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.5) + 
  xlab('Tissue anatomy') +
  ylab('Count of cell types enriched with one or more RTEs') +
  scale_y_continuous(limits=c(0, 20), breaks = seq(0, 20, by = 2))+
  scale_fill_manual(values=c('#00FF55','#FFFF00','#2317D5'))

# Modify the appearance of plot:

p2 = p + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  #panel.grid.major = element_line(colour = "white"),
  panel.grid.minor = element_line(colour = "white"),
  # Change axis line
  axis.line = element_line(colour = "black"), 
  # Change background to white
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = c('#f2f2f2')),
  # to add legend box background color
  legend.background = element_rect(fill="white", size=0.05, linetype="solid"),
  legend.position = c(0.8, 0.9),
  legend.direction = ("horizontal"),
  # To remove legand title:
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=8, face="bold"),
  axis.title.x = element_text(size=11, face="bold"),
  axis.title.y = element_text(size=11, face="bold"),
  axis.text = element_text(size = 9, face = "bold"),
  axis.text.x = element_text(angle = 45, hjust=0.9,vjust=0.95))

p2

# To save the plot:
# click on Export > Save as image > change Width to 625 and Height to 550 > Change image name > Save 

hg = merge(L1names,Alunames, by = "EID", all.x = TRUE, all.y = T)
hv = merge(hg,SVAnames,by = "EID", all.x = TRUE, all.y = T)
write.table(hv, file = "Enrichment_in_Euchromatin_regions/Results/results2807.txt", sep="\t", row.names=F,col.names=T, quote = F)
