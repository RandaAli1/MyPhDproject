# Displaying the distribution of reference and non-reference retrotransposable elements (RTEs) in enhancer and enhancer-free regions of the human genome.

# This code displays the frequency distribution of RTEs in enhancer vs. enhancer-free regions of the human genome.
# Enhancer regions are based on GeneHancer database, obtained from the UCSC table browser. 
# The code is executable from RStudio, requiring the ggplot2 and ggpubr libraries.

# Change Mydir to path to MyAnalysis:

Mydir = c("Path_to_MyAnalysis/")

setwd(paste0(Mydir))

library(ggplot2)

# Load input file: 

Sheet1 = read.table("Functional_analysis/Results/EnhancerRegionsFrequencies.txt", header = T)

data = data.frame(Group = Sheet1$Region, Category = Sheet1$Category, Freq = round(Sheet1$Frequency, digits = 2), Type = Sheet1$RTE_Type)

L1 = subset(data,Type== "L1")
Alu = subset(data,Type == "Alu")
SVA = subset(data,Type == "SVA")

## Plot the Distribution of reference and non-reference L1s in enhancer vs. enhancer-free regions:

p = ggplot(data=L1, aes(x=Group, y=Freq, fill=Category)) +
  # To add a black line around the bars and change the width
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.5) +
  #To add frequencies inside the bar
  geom_text(aes(label=paste0(L1$Freq,"%")), vjust=-0.9,color="black", position = position_dodge(0.5), size=2.5)+
  xlab('Enhancer region') +
  ylab('Frequency of L1s (%)') +
  scale_y_continuous(limits=c(0, 100), breaks = seq(0, 100, by = 20))

# To change the colur:
p1 = p + scale_fill_manual(values=c('#FFE600','#11D118'))

# Modify the appearance of plot:

p2 = p1 + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black"), 
  # Change background to white
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  # to add legend box background color
  legend.background = element_rect(fill="white", size=0.05, linetype="solid"),
  legend.position = "none",
  # To remove legand title:
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=8, face="bold"),
  axis.title.x = element_text(size=10, face="bold"),
  axis.title.y = element_text(size=10, face="bold"),
  axis.text = element_text(size = 9, face = "bold"))

L1_enhancer <- p2

# Plot the Distribution of reference and non-reference Alus in enhancer vs. enhancer-free regions: 

p = ggplot(data=Alu, aes(x=Group, y=Freq, fill=Category)) +
  # To add a black line around the bars and change the width
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.5) +
  #To add frequencies inside the bar
  geom_text(aes(label=paste0(Alu$Freq,"%")), vjust=-0.9,color="black", position = position_dodge(0.5), size=2.5)+
  xlab('Enhancer region') +
  ylab('Frequency of Alus (%)')+
  scale_y_continuous(limits=c(0, 100), breaks = seq(0, 100, by = 20))

# To change the colur:
p1 = p + scale_fill_manual(values=c('#FFE600','#11D118'))

# Modify the appearance of plot:

p2 = p1 + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black"), 
  # Change background to white
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  # to add legend box background color
  legend.background = element_rect(fill="gray", size=0.05, linetype="solid"),
  legend.position = c(0.8, 0.8),
  # To remove legand title:
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=8, face="bold"),
  axis.title.x = element_text(size=10, face="bold"),
  axis.title.y = element_text(size=10, face="bold"),
  axis.text = element_text(size = 9, face = "bold"))

Alu_enhancer <- p2

# Plot the Distribution of reference and non-reference SVAs in enhancer vs. enhancer-free regions: 

p = ggplot(data=SVA, aes(x=Group, y=Freq, fill=Category)) +
  # To add a black line around the bars and change the width
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.5) +
  #To add frequencies inside the bar
  geom_text(aes(label=paste0(SVA$Freq,"%")), vjust=-0.9,color="black", position = position_dodge(0.5), size=2.5)+
  xlab('Enhancer region') +
  ylab('Frequency of SVAs (%)')+
  scale_y_continuous(limits=c(0, 100), breaks = seq(0, 100, by = 20))

# To change the colur:
p1 = p + scale_fill_manual(values=c('#FFE600','#11D118'))
p1
# Modify the appearance of plot:

p2 = p1 + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black"), 
  # Change background to white
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  # to add legend box background color
  legend.background = element_rect(fill="white", size=0.05, linetype="solid"),
  legend.position = "none",
  # To remove legand title:
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=8, face="bold"),
  axis.title.x = element_text(size=10, face="bold"),
  axis.title.y = element_text(size=10, face="bold"),
  axis.text = element_text(size = 9, face = "bold"))

SVA_enhancer <- p2

# To plot all the figures in one:

#install.packages("ggpubr")
library("ggpubr")

plot_all = ggarrange(L1_enhancer, Alu_enhancer, SVA_enhancer, labels = c("A","B","C"), font.label = list(size = 10), ncol = 2, nrow = 2)

plot_all

# To save the plot:
# click on Export > Save as image > change Width to 750 and Height to 450 > Change image name > Save 
