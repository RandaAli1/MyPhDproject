# Z-statistics calculation for the overrepresentation of RTEs in gene regions. 

# This code analyses the overrepresentation of RTEs in gene regions by 
# comparing the frequency distribution of reference and non-reference RTEs in 
# gene regions against the the population mean of 1,000 sets of random 
# iterations using z-statstics.

# Change Mydir to path to MyAnalysis:

Mydir = c("Path_to_MyAnalysis/")

setwd(paste0(Mydir))

# Load frequency distribution of reference and non-reference RTEs in each genomic region:
data = read.table("Enrichment_In_Functional_regions/Results/GeneRegionsFrequencies.txt", header = T)

# Load L1 and corrosponding random results:
results_L1 = read.delim("Enrichment_In_Functional_regions/Results/L1_1kRandomInsertionsIngeneregions.txt", header = T)
L1 = subset(data,RTE_type=="L1")
Ref_L1 = L1$Freq_Reference
NonRef_L1 = L1$Freq_Non.reference

# Load Alu and corrosponding random results:
results_Alu = read.delim("Enrichment_In_Functional_regions/Results/AluY_1kRandomInsertionsIngeneregions.txt", header = T)
Alu = subset(data,RTE_type=="Alu")
Ref_Alu = Alu$Freq_Reference
NonRef_Alu = Alu$Freq_Non.reference

# Load SVA and corrosponding random results:

results_SVA = read.delim("Enrichment_In_Functional_regions/Results/SVAEF_1kRandomInsertionsIngeneregions.txt", header = T)
SVA = subset(data,RTE_type== "SVA")
Ref_SVA <- SVA$Freq_Reference
NonRef_SVA <- SVA$Freq_Non.reference


#### L1 Z-score and P-value calculation: 

# Remove NA entries from random dataset:
results_L1 = results_L1[complete.cases(results_L1 * 0), , drop=FALSE]

# Set loop for calculating summary statistics and standard deviation of random dataset:

results = NA
standdv = NA

for (i in 1:3) {
  v1 = results_L1[,i]
  sumry = summary(as.numeric(as.character(v1)),na.rm = TRUE)
  deviation = sd(as.numeric(as.character(v1)),na.rm = TRUE)
  results <<- rbind(results, sumry)
  standdv <<- rbind(standdv, deviation)
}

# Remove NA from results:
results1 = results[complete.cases(results * 0), , drop=FALSE]
standdv1 = standdv[complete.cases(standdv * 0), , drop=FALSE]

# Create loop for calculating Z-score and P-value for Reference L1s:

zscore = NA
pvalue1sided = NA

for (i in 1:3) {
  z = ((Ref_L1[i]-results1[i,4])/standdv1[i,])
  pvalue = 2*pnorm(-abs(z))
  zscore <<- rbind(zscore, z)
  pvalue1sided <<- rbind(pvalue1sided, pvalue)
}

# Remove NA and round to 3 significant figures (sf):
zscore = round(zscore[complete.cases(zscore * 0), , drop=FALSE], digits = 3)
pvalue1sided = format.pval(pvalue1sided[complete.cases(pvalue1sided * 0), , drop=FALSE], digits = 3, eps = 0)

# Bind z-score and P-value in dataframe:

results_Zscore <- NA
results_Pval <- NA

results_Zscore <<- rbind(results_Zscore, zscore)
results_Pval <<- rbind(results_Pval, pvalue1sided)

refL1_results <- cbind("Reference","L1",zscore,pvalue1sided)
colnames(refL1_results) <- c("RTE_Category","RTE_type","Z_Score","P-value")
refL1_results = noquote(refL1_results)
rownames(refL1_results) <- c("Intergenic", "Intronic","Exonic")


# Calculate Z-score and P-value for Non-reference L1s:

zscore = NA
pvalue1sided = NA

for (i in 1:3) {
  z = ((NonRef_L1[i]-results1[i,4])/standdv1[i,])
  pvalue = 2*pnorm(-abs(z))
  zscore <<- rbind(zscore, z)
  pvalue1sided <<- rbind(pvalue1sided, pvalue)
}

# Remove NA and round to 3 sf:
zscore = round(zscore[complete.cases(zscore * 0), , drop=FALSE], digits = 3)
pvalue1sided = format.pval(pvalue1sided[complete.cases(pvalue1sided * 0), , drop=FALSE], digits = 3, eps = 0)

# Bind results
NonrefL1_results <- cbind("Non-Reference","L1",zscore,pvalue1sided)
NonrefL1_results = noquote(NonrefL1_results)
rownames(NonrefL1_results) <- c("Intergenic", "Intronic","Exonic")
colnames(NonrefL1_results) = NULL

####### Below code (code lines 107-251) is repeating the Z-score analysis for Alu and SVA RTEs:
 
## Alu Z-score and P-value calculation: 

# Remove NA entries from random dataset:
results_Alu = results_Alu[complete.cases(results_Alu * 0), , drop=FALSE]

# Set loop for calculating summary statistics and standard deviation of random dataset:

results = NA
standdv = NA

for (i in 1:3) {
  v1 = results_Alu[,i]
  sumry = summary(as.numeric(as.character(v1)),na.rm = TRUE)
  deviation = sd(as.numeric(as.character(v1)),na.rm = TRUE)
  results <<- rbind(results, sumry)
  standdv <<- rbind(standdv, deviation)
}

# Remove NA from results:
results1 = results[complete.cases(results * 0), , drop=FALSE]
standdv1 = standdv[complete.cases(standdv * 0), , drop=FALSE]

# Create loop for calculating Z-score and P-value for Reference Alus:

zscore = NA
pvalue1sided = NA

for (i in 1:3) {
  z = ((Ref_Alu[i]-results1[i,4])/standdv1[i,])
  pvalue = 2*pnorm(-abs(z))
  zscore <<- rbind(zscore, z)
  pvalue1sided <<- rbind(pvalue1sided, pvalue)
}

# Remove NA and round to 3 sf:
zscore = round(zscore[complete.cases(zscore * 0), , drop=FALSE], digits = 3)
pvalue1sided = format.pval(pvalue1sided[complete.cases(pvalue1sided * 0), , drop=FALSE], digits = 3, eps = 0)

# Bind z-score and P-value in dataframe:

results_Zscore <- NA
results_Pval <- NA

results_Zscore <<- rbind(results_Zscore, zscore)
results_Pval <<- rbind(results_Pval, pvalue1sided)

refAlu_results <- cbind("Reference","Alu",zscore,pvalue1sided)
refAlu_results = noquote(refAlu_results)
rownames(refAlu_results) <- c("Intergenic", "Intronic","Exonic")
colnames(refAlu_results) = NULL

# Calculate Z-score and P-value for Non-reference Alus:

zscore = NA
pvalue1sided = NA

for (i in 1:3) {
  z = ((NonRef_Alu[i]-results1[i,4])/standdv1[i,])
  pvalue = 2*pnorm(-abs(z))
  zscore <<- rbind(zscore, z)
  pvalue1sided <<- rbind(pvalue1sided, pvalue)
}

# Remove NA and round to 3 sf:
zscore = round(zscore[complete.cases(zscore * 0), , drop=FALSE], digits = 3)
pvalue1sided = format.pval(pvalue1sided[complete.cases(pvalue1sided * 0), , drop=FALSE], digits = 3, eps = 0)

# Bind results
NonrefAlu_results <- cbind("Non-Reference","Alu",zscore,pvalue1sided)
NonrefAlu_results = noquote(NonrefAlu_results)
rownames(NonrefAlu_results) <- c("Intergenic", "Intronic","Exonic")
colnames(NonrefAlu_results) = NULL

#SVA Z-score and P-value calculation: 

# Remove NA entries from random dataset:
results_SVA = results_SVA[complete.cases(results_SVA * 0), , drop=FALSE]

# Set loop for calculating summary statistics and standard deviation of random dataset:

results = NA
standdv = NA

for (i in 1:3) {
  v1 = results_SVA[,i]
  sumry = summary(as.numeric(as.character(v1)),na.rm = TRUE)
  deviation = sd(as.numeric(as.character(v1)),na.rm = TRUE)
  results <<- rbind(results, sumry)
  standdv <<- rbind(standdv, deviation)
}

# Remove NA from results:
results1 = results[complete.cases(results * 0), , drop=FALSE]
standdv1 = standdv[complete.cases(standdv * 0), , drop=FALSE]

# Create loop for calculating Z-score and P-value for Reference SVAs:

zscore = NA
pvalue1sided = NA

for (i in 1:3) {
  z = ((Ref_SVA[i]-results1[i,4])/standdv1[i,])
  pvalue = 2*pnorm(-abs(z))
  zscore <<- rbind(zscore, z)
  pvalue1sided <<- rbind(pvalue1sided, pvalue)
}

# Remove NA and round to 3 sf:
zscore = round(zscore[complete.cases(zscore * 0), , drop=FALSE], digits = 3)
pvalue1sided = format.pval(pvalue1sided[complete.cases(pvalue1sided * 0), , drop=FALSE], digits = 3, eps = 0)

# Bind z-score and P-value in dataframe:

results_Zscore <- NA
results_Pval <- NA

results_Zscore <<- rbind(results_Zscore, zscore)
results_Pval <<- rbind(results_Pval, pvalue1sided)

refSVA_results <- cbind("Reference","SVA",zscore,pvalue1sided)
refSVA_results = noquote(refSVA_results)
rownames(refSVA_results) <- c("Intergenic", "Intronic","Exonic")
colnames(refSVA_results) = NULL

# Calculate Z-score and P-value for Non-reference SVAs:

zscore = NA
pvalue1sided = NA

for (i in 1:3) {
  z = ((NonRef_SVA[i]-results1[i,4])/standdv1[i,])
  pvalue = 2*pnorm(-abs(z))
  zscore <<- rbind(zscore, z)
  pvalue1sided <<- rbind(pvalue1sided, pvalue)
}

# Remove NA and round to 3 sf:
zscore = round(zscore[complete.cases(zscore * 0), , drop=FALSE], digits = 3)
pvalue1sided = format.pval(pvalue1sided[complete.cases(pvalue1sided * 0), , drop=FALSE], digits = 3, eps = 0)

# Bind results
NonrefSVA_results <- cbind("Non-Reference","SVA",zscore,pvalue1sided)
NonrefSVA_results = noquote(NonrefSVA_results)
rownames(NonrefSVA_results) <- c("Intergenic", "Intronic","Exonic")
colnames(NonrefSVA_results) = NULL

# Bind all Z-score and P-value results for reference and non-reference L1s, Alus, and SVAs:

all_results = noquote(rbind(refL1_results,NonrefL1_results,refAlu_results,NonrefAlu_results,
                            refSVA_results,NonrefSVA_results))

all_results = noquote(cbind(Region = rownames(all_results), all_results))

# Write results to file:

write.table(all_results, file = "Enrichment_In_Functional_regions/Results/Zstatistics_RTEsInGeneRegionsVsRandom.txt", sep="\t", row.names=F,col.names=T, quote = F)

