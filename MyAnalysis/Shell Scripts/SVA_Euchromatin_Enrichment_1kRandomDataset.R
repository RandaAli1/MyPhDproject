# Analysing the distribution of a control dataset based on the size of non-reference SVA elements
# in euchromatin regions

# This code creates a control dataset of random insertions. The control dataset is used to subsample 1,000 sets of 
# 1,888 random insertions corrosponding to the size of non-reference curated SVA elements. Next, the frequency 
# of random insertions interrupting euchromatin regions in each of the 127 H3K4me3 maps obtained from the roadmap website is calculated to create a background
# distribution against which the enrichment of SVA elements in euchromatin regions can be determined via the empirical P-value. 

# Change Mydir to path to MyAnalysis:

Mydir = c("Path_to_MyAnalysis/")

setwd(paste0(Mydir))

# Create a function for generating the control dataset:
bedtoolsRandomFun = function(functionstring="bedtools random",SampleSize,opt.string="") {
  a.file =tempfile()
  out = tempfile()
  
  command = paste(functionstring,"-g Enrichment_in_Euchromatin_regions/Raw_data/chr_lengthbp.genome -n" ,SampleSize,opt.string,">",out,sep="\t")
  cat(command,"\n")
  try(system(command))
  
  res <- read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}

# Run function to generate the control dataset:

v = bedtoolsRandomFun("bedtools random","100000")

# Create a function for subdampling a set of random insertions:

EnhancersFrequency <- function(data, SampleSize){
  d <- data[sample(dim(data)[[1]],SampleSize),] 
  return(d)
}

# Create function for overlapping the set of random insertions with Roadmap H3K4me3 peaks:

bedTools.2in <- function(functionstring="bedtools intersect",d,opt.string="")
{
  a.file <- tempfile()
  out   <- tempfile()
  write.table(d,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  command <- paste(functionstring,"-a",a.file,"-b Enrichment_in_Euchromatin_regions/Raw_data/RoadmapAll_H3K4me3.broadPeak_sorted_27062019.bed -wa -wb | sort -k1,1V -k2,2n | cut -f13 | sort -k1,1V -k2,2n | uniq -c",opt.string,">",out,sep="\t")
  cat(command,"\n")
  try(system(command))
  
  res <- read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}

# Run one set to creare a dataframe of "EID" names

v1 <- EnhancersFrequency(v,1888)
results_1 <- bedTools.2in("bedtools intersect",v1)

Results2 = data.frame(results_1$V2)
colnames(Results2) = c("EID")

# Set R to the number of iterations for the loop function:

R = 1000

# Run loop function for investigating the frequency of a random insertions 
# in the euchromatin regions of each cell type:

for (i in 1:R) {
  print(i)
  v1 <- EnhancersFrequency(v,1888) 
  results_1 <- bedTools.2in("bedtools intersect",v1)
  results_n = data.frame(results_1$V1)
  colnames(results_n) = c(i)
  Results2 <<- cbind(Results2, results_n)
}

# Extract results in a tab delimited text file:

write.table(Results2, file = "Enrichment_in_Euchromatin_regions/Results/RandomSVA_1k_H3K4Euchromatin.txt", sep="\t", row.names=F,col.names=T, quote = F)

