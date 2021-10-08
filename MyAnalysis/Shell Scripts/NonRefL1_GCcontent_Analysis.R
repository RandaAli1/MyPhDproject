# Analysing the Guanine-cytosine (GC) content of genomic features.

# This code analysis the GC content of non-reference L1s in 20 kilobase (kb) windows.
# It is executable as a script from the command line. 
# The code counts the number of nucleotides within each feature. 
# Next, the percentage of GC nucleotides is calculated. 

# Change Mydir to path to MyAnalysis:

Mydir = c("Path_to_MyAnalysis/")

setwd(paste0(Mydir))

# Load input file:

locs <- read.delim("GC_content_analysis/Raw_data/RTEdb_L1Hs_all_200bpmerged_20kbWindow.bed", header=F, sep="\t")

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
  command <- paste("bedtools nuc","-fi","GC_content_analysis/Raw_data/hg19_NoChr_08112017.fa","-bed",b.file,">",out,sep="\t");
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
  Gc.content <- (round(c((results$V15+results$V16)/(results$V14+results$V15+results$V16+results$V17)),3))*100
  return(Gc.content)
}

# Run the GC content function:

results <- GCcontent(locs, len)
summaryGC <- summary(results, na.rm=TRUE)

# Write GC content results in a txt file: 

write.table(results, file = "GC_content_analysis/Results/RTEdb_L1Hs_GCresults20kb_26072019.txt", sep="\t")
