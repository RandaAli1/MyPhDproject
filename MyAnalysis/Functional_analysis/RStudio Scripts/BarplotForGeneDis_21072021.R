# Displaying the distribution of reference and non-reference retrotransposable elements (RTEs) in genic and intergenic regions of the human genome.

# This code displays the frequency distribution of RTEs in genic regions defined as intronic and exonic regions, as well as intergenic regions defined as regions between genes.
# Gene regions are based on the NCBI RefSeq genes, obtained from the UCSC table browser.
# The code is executable from RStudio, requiring the ggplot2 and ggpubr libraries.  

# Change Mydir to path to MyAnalysis:

Mydir = c("C:/Users/rali1/OneDrive - University of Plymouth/Thesis_2019/MyAnalysis/")

setwd(paste0(Mydir))

#install.packages("ggplot2")
#install.packages("ggpubr")

library(ggplot2)
library(plyr) 

# Load input file: 

Sheet1 = read.table("Functional_analysis/Results_Copy/GeneRegionsFrequencies_21072021.txt", header = T)

data = data.frame(Group = Sheet1$Region, Category = Sheet1$Category, Freq = round(Sheet1$Frequency, digits = 2), Type = Sheet1$RTE_Type)

neworder <- c("L1","Alu","SVA")

data2 <- arrange(transform(data,
                           Type=factor(Type,levels=neworder)),Type)

## Plot the Distribution of Reference and non-reference RTEs in genic vs. intergenic regions:

p = ggplot(data=data2, aes(x=Group, y=Freq, fill=Category)) +
  facet_wrap(~Type)+
  # To add a black line around the bars and change the width
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.5) +
  xlab('Genomic group') +
  ylab('Frequency of RTE (%)') +
  scale_y_continuous(limits=c(0, 70), breaks = seq(0, 70, by = 5))
p

# To order the Regions:
p1 = p + scale_fill_manual(values=c('#006FFF','#2059A4')) + scale_x_discrete(limits=c("Intergenic", "Intronic", "Exonic"))
p1

p2 = p1 + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black"), 
  # Change background to white
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = c('#F0F0F0')),
  # to add legend box background color
  legend.background = element_rect(fill="white", size=0.05, linetype="solid"),
  legend.position = c(0.9,0.9),
  # To remove legand title:
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=8, face="bold"),
  axis.title.x = element_text(size=11, face="bold"),
  axis.title.y = element_text(size=11, face="bold"),
  axis.text = element_text(size = 10, face = "bold"),
  strip.text = element_text(size=15))

p2

# To save the plot:
# click on Export > Save as image  
