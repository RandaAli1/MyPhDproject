# Regression analysis for the relationship between the density of Alus per chromosome and gene density

# This code conducts a simple linear regression analysis to determine whether there is 
# a significant relationship between the density of reference and non-reference Alus
# per chromosome and chromosomal gene densities. 

# Change Mydir to path to MyAnalysis:

Mydir = c("Path_to_MyAnalysis/")

setwd(paste0(Mydir))

#install.packages("ggpubr")
#install.packages("ggplot2")

library(ggpubr)
library(ggplot2)

# Load Input file:

Sheet1 <- read.csv("Regression_analysis/RegressionAnalysis_inputFile.txt", sep = "\t", fileEncoding="latin1")


chrname <- Sheet1$Chromosome               # Chromosome names 
chrlength <- Sheet1$ChrLength_Mb           # Chromosome size in megabases (Mb)
GeneDensity = Sheet1$Genes_Density.Mbp     # Gene density per chromosome

nonref_Aludensity = Sheet1$Nonref_Alu.Mb   # Density of non-reference Alus per chromosome
Ref_Aludensity = Sheet1$Ref_Alu.Mb         # Density of reference Alus per chromosome

# Regression analysis for non-reference Alus density vs. gene density

# Create regression model for gene density vs. density of Alu/chr:

linRegLine <- lm(nonref_Aludensity ~ GeneDensity)
summary(linRegLine)

# Calculate 95% prediction intervals:

temp_var <- predict(linRegLine, interval="prediction", level = 0.95)
new_df <- cbind(Sheet1, temp_var)

# Draw regression plot: 

p = ggplot(new_df, aes(x=GeneDensity,y=nonref_Aludensity)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  geom_text(aes(label=as.character(chrname)),hjust=0.5, vjust=-0.6, size = 3.5) +
  geom_smooth(method=lm, se = T) +
  #Add axis label
  xlab('Gene Density (genes/million bases)') +
  ylab('Non-reference Alus Density (Alus/million bases)') +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed") 

p

# Modify the appearance of plot:

p1 = p + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black"), 
  # Change background to white
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  axis.title.x = element_text(size=11, face="bold"),
  axis.title.y = element_text(size=11, face="bold"),
  axis.text = element_text(size = 9, face = "bold")
)

p1

## Regression analysis for reference Alus density vs. gene density

# Use same steps as non-reference Alus:

nrlinRegLine <- lm(Ref_Aludensity ~ GeneDensity)
summary(nrlinRegLine)

temp_var <- predict(nrlinRegLine, interval="prediction", level = 0.95)
new_df <- cbind(Sheet1, temp_var)

## 
nrp = ggplot(new_df, aes(x=GeneDensity,y=Ref_Aludensity)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  geom_text(aes(label=as.character(chrname)),hjust=0.5, vjust=-0.6, size = 3.5) +
  geom_smooth(method=lm, se = T) +
  #Add axis label
  xlab('Gene Density (genes/million bases)') +
  ylab('Reference Alus Density (Alus/million bases)') +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed") 

nrp

nrp1 = nrp + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black"), 
  # Change background to white
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  axis.title.x = element_text(size=11, face="bold"),
  axis.title.y = element_text(size=11, face="bold"),
  axis.text = element_text(size = 9, face = "bold")
)

nrp1

## Use ggpubr library to plot regression plot of reference next to non-reference L1s:


plot_L1 = ggarrange(nrp1, p1, labels = c("A","B"), ncol = 2, nrow = 1, 
                    font.label = list(size = 10), hjust = -0.8)
plot_L1




