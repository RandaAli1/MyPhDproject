# Analysing the distribution of a control dataset based on the size of non-reference AluY elements
# in enhancer regions.

# This code creates a control dataset of random insertions. The control dataset is used to subsample 1,000 sets of 
# 27,699 random insertions corrosponding to the size of non-reference curated AluY elements. Next, the frequency 
# of random insertions interrupting enhancer regions is calculated for each set, to create a background
# distribution against which the frequency of AluY elements in enhancer regions can be
# compared via z statistics. 

args = commandArgs()
setwd(args[6])

# Create a function for generating the control dataset:

bedtoolsRandomFun = function(functionstring="bedtools random",SampleSize,opt.string="") {
a.file =tempfile()
out = tempfile()

command = paste(functionstring,"-g Enrichment_In_Functional_regions/Raw_data/chr_lengthbp.genome -n" ,SampleSize,opt.string,">",out,sep="\t")
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

# Create function for overlapping the set of random insertions with GeneHancer regions:

bedTools.2in<- function(functionstring="bedtools intersect",d,opt.string="")
{
  a.file <- tempfile()
  out   <- tempfile()
  write.table(d,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  command <- paste(functionstring,"-a",a.file,"-b Functional_analysis/GeneHancer/AllChromosomes_GeneHancerRegulators_26062019.bed -wao | sort -k1,1V -k2,2n",opt.string,">",out,sep="\t")
  cat(command,"\n")
  try(system(command))
  
  res <- read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}

Results1 = NA

# Set R to the number of iterations for the loop function:

R = 1000

# Run loop function for investigating the distribushion of a random set of insertions 
# in gene regions, and calculating the frequency of insertions in enhancer regions

for (i in 1:R) {
	 
	  v1 <- EnhancersFrequency(v,27699) 
	 results <- bedTools.2in("bedtools intersect",v1)
	 len <- nrow(results)
     overlappingEnhancers = subset(results, results$V8 != "-1")
     overlapping_enhancers = length(unique(overlappingEnhancers[[2]]))
     frequency_enhancers <- c((overlapping_enhancers/len)*100)
     Results1 <<- rbind(Results1, frequency_enhancers)
}

write.table(Results1, file = "Enrichment_In_Functional_regions/Results/AluY_1kRandomInsertions_Enhancerregions.txt", sep="\t", row.names=F,col.names=T, quote = F)
