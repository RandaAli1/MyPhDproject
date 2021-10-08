# Regression analysis for the relationship between the frequency of L1s per chromosome and chromosome size

# This code conducts a simple linear regression analysis to determine whether there is 
# a significant relationship between the frequency of reference and non-reference L1s per chromosome and chromosome size. 


# Change Mydir to path to MyAnalysis:

Mydir = c("Path_to_MyAnalysis/")

setwd(paste0(Mydir))

library(ggplot2)
library(ggpubr)  

# Load Input file:

Sheet1 <- read.csv("Regression_analysis/RegressionAnalysis_inputFile.txt", sep = "\t", fileEncoding="latin1")

chrname <- Sheet1$Chromosome                # Chromosome namee
chrlength <- Sheet1$ChrLength_Mb            # Chromosome size in megabases (Mb)
Ref <- Sheet1$RefL1.chr_Freq                # Frequency of reference L1s per chromosome 
NonRef <- Sheet1$NonRefL1.chr_Freq          # Frequency of non-reference L1s per chromosome

# Regression analysis for frequency of reference L1s per chromosome vs. chromosome length:

# Create regression model for total L1s vs. chromosome length:

linRegLine <- lm(Ref ~ chrlength)
summary(linRegLine)

# Calculate 95% preditction intervals:

temp_var <- predict(linRegLine, interval="prediction", level = 0.95)
new_df <- cbind(Sheet1, temp_var)

p = ggplot(new_df, aes(x=chrlength,y=Ref)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  geom_text(aes(label=as.character(chrname)),hjust=0.5, vjust=-0.6) +
  geom_smooth(method=lm, se = T) +
  #Add axis label
  xlab('Chromosome length (Mb)') +
  ylab('Reference L1s (%)') +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed")+
  scale_y_continuous(limits=c(0, 10), breaks = seq(0, 10, by = 2))
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
  axis.title.x = element_text(size=14, face="bold"),
  axis.title.y = element_text(size=14, face="bold"),
  axis.text = element_text(size = 11, face = "bold")
)
p1

ref = p1

# Regression analysis for frequency of non-reference L1s per chromosome vs. chromosome length:

linRegLine <- lm(NonRef ~ chrlength)
summary(linRegLine)

temp_var <- predict(linRegLine, interval="prediction", level = 0.95)
new_df <- cbind(Sheet1, temp_var)

p = ggplot(new_df, aes(x=chrlength,y=NonRef)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  geom_text(aes(label=as.character(chrname)),hjust=0.5, vjust=-0.6) +
  geom_smooth(method=lm, se = T) +
  # Add axis label: 
  xlab('Chromosome length (Mb)') +
  ylab('Non-reference L1Hs (%)') +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed")+
  scale_y_continuous(limits=c(0, 10), breaks = seq(0, 10, by = 2))
p

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
  axis.title.x = element_text(size=14, face="bold"),
  axis.title.y = element_text(size=14, face="bold"),
  axis.text = element_text(size = 11, face = "bold")
)
p1

nonref = p1

## Use ggpubr library to plot regression plot of reference next to non-reference L1s:

plot_L1 = ggarrange(ref, nonref, labels = c("A","B"), ncol = 2, nrow = 1)
plot_L1

