# Generating Linkage disequilibrium (LD) blocks for 1,000 sets of SNPs

# This code takes an output file from SNPsnap and generates LD-blocks around each
# SNP using its first 5' and last 3' proxy SNPs.
# Note: Only takes SNPs identified in European populations, and located on chromosomes 1-22


# Load required library:
library("dplyr")

# Change Mydir to path to RTE_enrichment_in_GWAS_loci:

Mydir = c('Path_To_RTE_enrichment_in_GWAS_loci/')
setwd(paste0(Mydir))

# Load SNPsnap matches:

matchedSNPs.annotated <- read.delim("DataFile/Demo_matched_snps.txt", header = F) 

# Load .bim file for all SNPs identified in European super population:

onekPGph3 <- read.delim("DataFile/allchromosome_european_1kPGph3.bim", header = F, sep = "\t")
colnames(onekPGph3) <- c("CHR", "BP", "SNP","A1","A2")
onekPGph3[,3] <- as.character(onekPGph3[,3])

# Load HLA region file:

HLA_reg <- read.delim("DataFile/HLA_region.bed", header = F, sep = "\t")

# The function to create the LD-blocks 1k times:

LD_blocks_1k <- function (x) {
  matching.set <- x[i] # extract the set 
  # Seperate the chr:pos format of SNPsnap matches into tab delimited and sorted file:
  file_randomSNP <-function(x)
  {
    a.file <- tempfile()
    out   <- tempfile()
    options(scipen =99) 
    write.table(matching.set,file=a.file,quote=F,col.names=F,row.names=F)
    command <- paste("tr ':' '\t' <",a.file,"|sort -k1,1V -k2,2n >",out)
    cat(command,"\n")
    try(system(command))
    res <- read.table(out,header=F)
    unlink(a.file);unlink(out)
    return(res) }
  randomSNPs_loci <- file_randomSNP(matching.set)
  colnames(randomSNPs_loci) <- c("CHR", "BP") 
  
  # Matched SNPs rs# with SNP loci:
  
  SNP_rsID <- merge(randomSNPs_loci, onekPGph3, by=c("CHR", "BP"))
  SNP_rsIDordered <- arrange(SNP_rsID, CHR, BP)
  #Extract matched SNPs rs#
  rsId_only <- as.data.frame(SNP_rsIDordered$SNP)
  colnames(rsId_only) <- NULL 
  rsId_only[,1] <- as.character(rsId_only[,1])
  rsId_only = unique(rsId_only)
  
  # Run PLINK command to obtain tagging SNPs:
  
  taggingSNPs <-function(x)
  {
    a.file <- tempfile()
    options(scipen =99) 
    write.table(rsId_only,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
    command <- paste("plink --bfile allchromosome_european_1kPGph3 --show-tags",a.file,"--list-all --tag-r2 0.8 --noweb --out taggedSNPs")
    cat(command,"\n")
    try(system(command))
    unlink(a.file)
    return("done") }
  taggingSNPs(x)
  
  # Processing the output file to create the LD block file:
  
  taggedSNPs_results <- read.table("taggedSNPs.tags.list", header = T) 
  taggedSNPs_results[,1] <- as.character(taggedSNPs_results[,1])
  nrow(unique(rsId_only))
  
  #Seperate results into groups as follows: not tagged, no left tags, no right tags, all matched with #left/right for the median:
  #No tagging SNPs:
  tagged_nozero <- taggedSNPs_results[!(taggedSNPs_results$NTAG == "0"),] 
  tagged_zero <- taggedSNPs_results[(taggedSNPs_results$NTAG == "0"),]
  #No left:
  tagged_noleft <- tagged_nozero[(tagged_nozero$BP==tagged_nozero$LEFT),]
  #No right:
  tagged_noright <- tagged_nozero[(tagged_nozero$BP==tagged_nozero$RIGHT),]
  #Defined LD block:
  tagged_withright <- tagged_nozero[!(tagged_nozero$BP==tagged_nozero$RIGHT),]
  tagged_defined <- tagged_withright[!(tagged_withright$BP==tagged_withright$LEFT),]
  nrow(tagged_defined)
  #Define median size:
  median_bp <- (median(tagged_defined$KBSPAN))*1000
  half_median <- round(median_bp/2, digits = 0)
  
  #Add half the median to the left,right,and left/right for tagged_zero:
  #left
  tagged_mediantoleft_file <- data.frame (CHR = tagged_noleft$CHR, Start = tagged_noleft$LEFT, End = tagged_noleft$RIGHT, SNP_ID = tagged_noleft$SNP, SNP_loci =tagged_noleft$BP)
  bed.slop <-function(functionstring="bedtools slop",x,opt.string="")
  {
    a.file <- tempfile()
    out   <- tempfile()
    options(scipen =99) # not to use scientific notation when writing out
    write.table(x,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
    command <- paste(functionstring,"-i",a.file,"-g chr_lengthbp.genome -l",half_median,"-r 0",opt.string,">",out,sep="\t")
    cat(command,"\n")
    try(system(command))
    res3 <- read.table(out,header=F)
    unlink(a.file);unlink(out) 
    return(res3) }
  tagged_mediantoleft <- bed.slop("bedtools slop",tagged_mediantoleft_file)
  tagged_mediantoleft[,1] <- as.numeric(tagged_mediantoleft[,1])
  tagged_mediantoleft[,2] <- as.numeric(tagged_mediantoleft[,2])
  tagged_mediantoleft[,3] <- as.numeric(tagged_mediantoleft[,3])
  #right
  tagged_mediantorightfile <- data.frame (CHR = tagged_noright$CHR, Start = tagged_noright$LEFT, End = tagged_noright$RIGHT, SNP_ID = tagged_noright$SNP, SNP_loci =tagged_noright$BP)
  bed.slop <-function(functionstring="bedtools slop",x,opt.string="")
  {
    a.file <- tempfile()
    out   <- tempfile()
    options(scipen =99) # not to use scientific notation when writing out
    write.table(x,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
    command <- paste(functionstring,"-i",a.file,"-g chr_lengthbp.genome -r",half_median,"-l 0",opt.string,">",out,sep="\t")
    cat(command,"\n")
    try(system(command))
    res3 <- read.table(out,header=F)
    unlink(a.file);unlink(out) 
    return(res3) }
  tagged_mediantoright <- bed.slop("bedtools slop",tagged_mediantorightfile)
  tagged_mediantoright[,1] <- as.numeric(tagged_mediantoright[,1])
  tagged_mediantoright[,2] <- as.numeric(tagged_mediantoright[,2])
  tagged_mediantoright[,3] <- as.numeric(tagged_mediantoright[,3])
  #left/right
  tagged_mediantozerofile <- data.frame (CHR = tagged_zero$CHR, Start = tagged_zero$LEFT, End = tagged_zero$RIGHT, SNP_ID = tagged_zero$SNP, SNP_loci = tagged_zero$BP)
  bed.slop <-function(functionstring="bedtools slop",x,opt.string="")
  {
    a.file <- tempfile()
    out   <- tempfile()
    options(scipen =99) # not to use scientific notation when writing out
    write.table(x,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
    command <- paste(functionstring,"-i",a.file,"-g chr_lengthbp.genome -b",half_median,opt.string,">",out,sep="\t")
    cat(command,"\n")
    try(system(command))
    res3 <- read.table(out,header=F)
    unlink(a.file);unlink(out) 
    return(res3) }
  tagged_mediantozero <- bed.slop("bedtools slop",tagged_mediantozerofile)
  tagged_mediantozero[,1] <- as.numeric(tagged_mediantozero[,1])
  tagged_mediantozero[,2] <- as.numeric(tagged_mediantozero[,2])
  tagged_mediantozero[,3] <- as.numeric(tagged_mediantozero[,3])
  #Combined defined, left, right, left/right and untagged into one file (sorted):
  #tagged_defined in the same format:
  tagged_defined_formated <- data.frame (CHR = tagged_defined$CHR, Start =tagged_defined$LEFT, End = tagged_defined$RIGHT, SNP_ID = tagged_defined$SNP, SNP_loci = tagged_defined$BP)
  tagged_defined_formated[,1] <- as.numeric(tagged_defined_formated[,1])
  tagged_defined_formated[,2] <- as.numeric(tagged_defined_formated[,2])
  tagged_defined_formated[,3] <- as.numeric(tagged_defined_formated[,3])
  #Combine files:
  all_SNPLDblock <- bind_rows(tagged_defined_formated,tagged_mediantoleft,tagged_mediantoright,tagged_mediantozero)
  all_SNPLDblock_ord <- arrange(all_SNPLDblock, CHR, Start)
  colnames(all_SNPLDblock_ord) <- NULL 
  #Merge LD-blocks:
  bed.merge <-function(functionstring="bedtools merge",all_SNPLDblock_ord,opt.string="")
  {
    a.file <- tempfile()
    out   <- tempfile()
    options(scipen =99) # not to use scientific notation when writing out
    write.table(all_SNPLDblock_ord,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
    command <- paste(functionstring,"-i",a.file,"-c 4,5 -o distinct,distinct",opt.string,">",out,sep="\t")
    cat(command,"\n")
    try(system(command))
    res3 <- read.table(out,header=F)
    unlink(a.file);unlink(out) 
    return(res3) }
  merged_SNPldblock <- bed.merge("bedtools merge",all_SNPLDblock_ord)
  #Remove HLA region:
  bed.subtract <-function(functionstring="bedtools subtract",merged_SNPldblock,HLA_reg,opt.string="")
  {
    a.file <- tempfile()
    b.file <- tempfile()
    out   <- tempfile()
    options(scipen =99) 
    write.table(merged_SNPldblock,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
    write.table(HLA_reg,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
    command <- paste(functionstring,"-a",a.file,"-b",b.file,"-A | sort -k1,1V -k2,2n",opt.string,">",out,sep="\t")
    cat(command,"\n")
    try(system(command))
    res <- read.table(out,header=F)
    unlink(a.file);unlink(b.file);unlink(out)
    return(res) }
  LDblock.noHLA <- bed.subtract("bedtools subtract",merged_SNPldblock,HLA_reg)
  return(LDblock.noHLA)
}
R <- 1000 
for (i in 1:R) {
  results1 <- LD_blocks_1k(matchedSNPs.annotated)
  myfile = paste0("1kLDBlocks/LDblock.noHLA_set_",i,".bed")
  write.table(results1, file = myfile, sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE)
}
