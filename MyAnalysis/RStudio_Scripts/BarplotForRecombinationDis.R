# Displaying the distribution of reference and non-reference retrotransposable elements (RTEs) in various recomination regions of the human genome.

# This code displays the frequency distribution of RTEs in various recombination regions defined as cold, 
# intermediate, and hot based on the the standardized recombination rate of the 2010 sex-averaged 
# Decode recombination map (Kong et al., 2010).
# The recombination map was obtained from the UCSC table browser.
# The code is executable from RStudio, requiring the ggplot2 and ggpubr libraries.  

# Change Mydir to path to MyAnalysis:

# Change Mydir to path to MyAnalysis:

Mydir = c("Path_to_MyAnalysis/")

setwd(paste0(Mydir))


# Load required libraries 
library(ggplot2)
library("ggpubr")


# Load input file: 

Sheet1 = read.table("Recombination_rate_analysis/Results/RecombinationFrequencies.txt", header = T)

data = data.frame(Group = Sheet1$Region, Ref = Sheet1$Freq_Reference, NonRef = Sheet1$Freq_Non.reference, Type = Sheet1$RTE_type)

## Plot the Distribution of Reference RTEs in the various recombination regions:

p = ggplot(data=data, aes(x=Group, y=Ref, fill=Type)) +
  # To add a black line around the bars and change the width
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.5) + 
  xlab('Recombination Region') +
  ylab('Frequency of Reference RTEs (%)') +
  scale_y_continuous(limits=c(0, 60), breaks = seq(0, 60, by = 5))

# To order the Regions:
p1 = p + scale_x_discrete(limits=c("Cold", "Intermediate", "Hot"), labels=c("Cold", "Intermediate", "Hot"))+
  scale_fill_manual(values=c('#00FF55','#FFFF00','#2317D5'))

# Modify the appearance of plot:

p2 = p1 + theme(
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
  legend.position = "none",
  # To remove legand title:
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=8, face="bold"),
  axis.title.x = element_text(size=14, face="bold"),
  axis.title.y = element_text(size=14, face="bold"),
  axis.text = element_text(size = 11, face = "bold"))

RTE_ref <- p2

# Plot the Distribution of Non-reference RTEs in the various recombination regions: 

p = ggplot(data=data, aes(x=Group, y=NonRef, fill=Type)) +
  # To add a black line around the bars and change the width
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.5) + 
  xlab('Recombination Region') +
  ylab('Frequency of Non-reference RTEs (%)') +
  scale_y_continuous(limits=c(0, 60), breaks = seq(0, 60, by = 5))
p
# To order the Regions:
p1 = p + scale_fill_manual(values=c('#00FF55','#FFFF00','#4237E5')) + scale_x_discrete(limits=c("Cold", "Intermediate", "Hot"))
p1

p2 = p1 + theme(
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
  legend.position = c(0.9, 0.9),
  # To remove legand title:
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=8, face="bold"),
  axis.title.x = element_text(size=14, face="bold"),
  axis.title.y = element_text(size=14, face="bold"),
  axis.text = element_text(size = 11, face = "bold"))

RTE_nonref <- p2

# To plot all the figures in one:

plot_all = ggarrange(RTE_ref, RTE_nonref, labels = c("A","B"), ncol = 2, nrow = 1)
plot_all

# To save the plot:
# click on Export > Save as image > change Width to 700 and Height to 500 > Change image name > Save 

