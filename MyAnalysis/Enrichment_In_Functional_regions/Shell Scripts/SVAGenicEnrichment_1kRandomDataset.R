# Analysing the distribution of a control dataset based on the size of non-reference SVA_E/F elements
# in genic vs. intergenic regions.


# This code creates a control dataset of random insertions. The control dataset is used to subsample 1,000 sets of 
# 1,888 random insertions corrosponding to the size of non-reference curated SVA_E/F elements. Next, the frequency 
# of random insertions interrupting an intergenic, intronic, or exonic region is calculated for each set to create a background
# distribution against which the frequency of SVA_E/F elements in each of the genomic regions can be
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

genic_subsampling <- function(data, SampleSize) {
  d <- data[sample(dim(data)[[1]],SampleSize),]
  return(d)
}

# Create function for overlapping the set of random insertions with RefSeq genes:

bedTools.2in <-function(functionstring="bedtools intersect",d,opt.string="")
{
  a.file <- tempfile()
  out   <- tempfile()
  options(scipen =99) 
  write.table(d,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  command <- paste(functionstring,"-a",a.file,"-b Functional_analysis/RefSeq/Refseq_geneIntervalsWithGeneNames_02102017.bed -loj | sort -k1,1V -k2,2n -k13 -k12 | bedtools merge -c 6,7,8,9,10,11,12,13,14 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct",opt.string,">",out,sep="\t")
  cat(command,"\n")
  try(system(command))
  
  res <- read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}


results_exon <- NA
results_intron <- NA
results_Intergenic <- NA

# Set R to the number of iterations for the loop function:

R = 1000

# Run loop function for investigating the distribushion of a random set of insertions 
# in gene regions, and calculating the frequency of insertions in intergenic, intronic, and exonic regions

for (i in 1:R) {
  
  v1 = genic_subsampling(v,1888)  
  results = bedTools.2in("bedtools intersect",v1)
  len <- nrow(results)
  intron <- nrow(results[grep("^Intron", results$V10), ])
  exon <- nrow(results[grep("Exon", results$V10, invert = F), ])
  Intergenic = nrow(results[grep("\\.", results$V10, invert = F), ])
  
  freq_exon <- (exon/len)*100
  freq_intron <- (intron/len)*100
  freq_Intergenic = (Intergenic/len)*100
  results_exon <<- rbind(results_exon, freq_exon)
  results_intron <<- rbind(results_intron, freq_intron)
  results_Intergenic <<- rbind(results_Intergenic,freq_Intergenic)
  
}

all_results <- cbind(results_Intergenic,results_intron,results_exon)
colnames(all_results) <- c("intergenic","intron","exon")

write.table(all_results, file = "Enrichment_In_Functional_regions/Results/SVAEF_1kRandomInsertionsIngeneregions.txt", sep="\t", row.names=F,col.names=T, quote = F)
