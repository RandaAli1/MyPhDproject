#Intersecting each set of random LD-blocks with Polymorphic RTEs and plot density plot

# This code overlaps the position of RTEs with 1,000 sets of LD-blocks, and creates a density plot 
# of the resulting frequency of overlap between the two variants.

# Change Mydir to path to MyAnalysis:

Mydir = c('Path_to_RTE_enrichment_in_GWAS_loci/')
setwd(paste0(Mydir))

# Load file of common RTEs:
# Note: Load files of different RTE classes seperatly to obtain a density plot for each RTE type

RTE = read.delim("DataFile/Demo_RTECommonInEUR.bed", header = F)

bed.Intersect <-function(functionstring="bedtoolsintersect",LDblock.noHLA,common_L1,opt.string="") {
  a.file <- tempfile()
  b.file <- tempfile()
  out   <- tempfile()
  options(scipen =99)
  write.table(LDblock.noHLA,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(common_L1,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  command <- paste(functionstring,"-a",a.file,"-b",b.file,"-wa -wb | cut -f2 | uniq -c | wc -l",opt.string,">",out,sep="\t")
  cat(command,"\n")
  try(system(command))
  res <- read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res) }

R = 5
results_freq <- NA
datalist = list()

for (i in 2:R) {
  # Load the random LD-blocks file for each set (results of GenerateLDblocks.R code):
  similar = read.delim(paste0("1kLDBlocks/LDblock.noHLA_set_",i,".bed"), header = F)
  nrow = nrow(similar)
  # Identify frequency of overlap using BEDtools function:
  # Set RTE as L1, Alu or SVA variable:
  total_overlap = bed.Intersect("bedtools intersect",RTE,similar)
  Freq = round((total_overlap/nrow)*100,digits = 3) 
  results_freqL1_i <- cbind(nrow,total_overlap,Freq)
  results_freqL1_i$i = i
  datalist[[i]] = results_freqL1_i }
results_freqL1 = do.call(rbind,datalist)
colnames(results_freqL1) = c("nrow","Total","Freq")

# Write results into table: 

write.table(results_freqL1, "RandomLDblocks_RTEelements.txt", 
            sep="\t", row.names=F,col.names=T, quote = F)

# Load Freq column into a variable:
f_1 <- as.data.frame(as.numeric(results_freqL1$Freq))
colnames(f_1) = c("Freq")
# Create a file to write the Density plot in: 
png(filename="SNP_RTE_Density.png")
plot(density(f_1$Freq))
# Add vertical line through the plot to indicate the fraction (%) TAS of LD-blocks with RTE elements and close the plot
# abline(v=3.90, col="red")
dev.off()
