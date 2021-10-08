
# Plot and statistical analysis for the GC distribution of RTEs

# This code plots the GC-content distribution of reference and non-reference L1s
# in 2% bin. The GC-content distribution for reference and non-reference L1s against 
# each other and against the genome were also compared using the Kolmogorov-Smirnov test.


# Change Mydir to path to MyAnalysis:

Mydir = c("C:/Users/rali1/OneDrive - University of Plymouth/Thesis_2019/MyAnalysis/")

setwd(paste0(Mydir))


# Load required library:
library(ggplot2)


# Section 1: Plots a linear graph overlaying the GC-content distributions: 

# Load the text files containing the GC contnent of the genome
# plus the reference and non-reference RTE windows: 

xnonref <- read.table("GC_content_analysis/Results/RTEdb_L1Hs_GCresults20kb_26072019.txt", sep = "\t", header = T)
xref <- read.table("GC_content_analysis/Results/RefL1PA2PA5_GCresults20kb_29052019.txt",sep = "\t", header = T)
xgenome <- read.table("GC_content_analysis/Results/GCcontent_wholeGenome_20kb_09022018.txt", sep = "\t", header = T)
xgenome1 <- na.omit(xgenome)

## Group the GC contents into 2% bins and calculate the frequency of
# each bin for each input file type

# Genome

his_Genome <- as.data.frame(table(cut(xgenome1$x, breaks=seq(20,75, by=2))))
his_Genome$Pct <- round((his_Genome$Freq / sum(his_Genome$Freq))*100,3)

# Non-Reference

his_nonRef <- as.data.frame(table(cut(xnonref$x, breaks=seq(20,75, by=2))))
his_nonRef$Pct <- round((his_nonRef$Freq / sum(his_nonRef$Freq))*100,3)

# Reference
his_Ref <- as.data.frame(table(cut(xref$x, breaks=seq(20,75, by=2))))
his_Ref$Pct <- round((his_Ref$Freq / sum(his_Ref$Freq))*100,3)

# Plot the GC-content distribushions:

p = ggplot(his_nonRef, aes(x=Var1,y=Pct, group = 1)) + 
  geom_line(aes(color="Non-Reference L1s"), size=1.5) +
  geom_line(data=his_Ref,aes(color="Reference L1s"), size=1.5)+
  geom_line(data=his_Genome,aes(color="Genome"), size=1.5)+
  xlab('GC content in 2% bins') +
  scale_color_manual(values=c('#bfbfbf','#ffcc00','#39ac39'))+
  scale_y_continuous(name="Frequency (%)", limits=c(0, 30), breaks = seq(0, 30, by = 5))
p 

# Rename X-axis using the 2% bin sizes:

p1 = p + scale_x_discrete(labels=c("20,22","22,24",	"24,26",	"26,28",	"28,30",	"30,32",
                                   "32,34",	"34,36",	"36,38",	"38,40",	"40,42",	"42,44",	"44,46",	"46,48",	"48,50",	
                                   "50,52",	"52,54",	"54,56",	"56,58",	"58,60",	"60,62",	"62,64",	"64,66",	"66,68",	
                                   "68,70",	"70,72",	"72,74"))
# Modify the appearance of plot:

p2 = p1 + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black"), 
  # Change background to white
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  # to add legend box background color
  legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
  legend.position = c(0.15, 0.9),
  # To remove legand title:
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=11, face="bold"),
  axis.title.x = element_text(size=12, face="bold"),
  axis.title.y = element_text(size=12, face="bold"),
  axis.text = element_text(size = 10, face = "bold"))

p2

# Section 2: Statistical analysis:

# Find the average GC-content: 

summary(xref$x)
summary(xnonref$x)
summary(xgenome$x)

# K-S test for ref vs. genome: 

ks.test(xgenome$x, xref$x)

# K-S test for nonref vs. genome:

ks.test(xgenome$x, xnonref$x)

# K-S test for ref vs. non-ref: 

ks.test(xref$x, xnonref$x)

