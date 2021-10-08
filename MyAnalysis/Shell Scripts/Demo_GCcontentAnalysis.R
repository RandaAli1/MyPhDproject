# Analysing the Guanine-cytosine (GC) content of genomic features.

# This code analysis the GC content of a sample hg19 genome and coordinates.
# It is executable as a script from the command line. 
# The code counts the number of nucleotides within each feature. 
# Next, the percentage of GC nucleotides is calculated. 

# Change Mydir to path to MyAnalysis:

Mydir = c("Path_to_MyAnalysis/")

setwd(paste0(Mydir))

# Load input file:

locs <- read.delim("GC_content_analysis/Raw_data/SampleFile.bed", header=F, sep="\t")

# Create a variable containing the number of rows of the input file:

len = nrow(locs)

# Function for Calculating the total number of each nucleotide
# using bedtools nuc: 

bedTools.gcContent <-function(x)
{
  b.file <- tempfile()
  out <- tempfile()
  options(scipen =99)
  write.table(locs,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  command <- paste("bedtools nuc","-fi","GC_content_analysis/Raw_data/Demo_hg19.fa","-bed",b.file,">",out,sep="\t");
  cat(command,"\n")
  try(system(command))
  res1 <- read.table(out)
  unlink(b.file);unlink(out);
  return(res1)
}

# Function for Calculating the percentage of G and C nucleotides: 

GCcontent <- function(data, SampleSize){
  d <- data[sample(dim(data)[[1]],SampleSize),]
  results <- bedTools.gcContent(d)
  Gc.content <- (round(c((results$V7+results$V8)/(results$V6+results$V7+results$V8+results$V9)),3))*100
  return(Gc.content)
}

# Run the GC content function:

results <- GCcontent(locs, len)
summaryGC <- summary(results, na.rm=TRUE)

# Write GC content results in a txt file: 

write.table(results, file = "GC_content_analysis/Results/Demo_GCcontent.txt", sep="\t")
