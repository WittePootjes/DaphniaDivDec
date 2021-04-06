#Packages
library(dada2); packageVersion("dada2")
library(vegan)
library(igraph)
library(SpiecEasi)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw())


setwd('C:/Users/Anna/Desktop/Daphnia_files/Complete analysis')
tt=read.table('meta_data.txt', sep ='\t', header=TRUE) # First I upload the meta-data table from the source

# The meta-data contains a lot of useless information which I am going to remove with a gsub command:
meta_data <- subset(tt, select = -c(1:3,5:6,8:10))
meta_data$experiment_title <- gsub('Illumina MiSeq sequencing; 16S rRNA sequencing of Daphnia magna microbiota - ', '', meta_data$experiment_title)
meta_data$experiment_title <- gsub('Illumina MiSeq sequencing; 16S rRNA sequencing of negative control - ', '', meta_data$experiment_title)
row.names(meta_data) <- meta_data$run_accession
meta_data$run_accession <- NULL # This leaves me with the sample names as row names and 1 column that defines the treatment

path <- "C:/Users/Anna/Desktop/Daphnia_files/Complete analysis/Sequences" # defining path as a directory to all my FASTQ files

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq.gz and SAMPLENAME_2.fastq.gz
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE)) #creating a variable "fnFs" that stores my forward reads fastq files
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE)) #creating a variable "FnRs" that stores my reverse reads fastq files

#Quality check of my sequences:
plotQualityProfile(fnFs)   #visualizing the phred quality profiles of the forward reads
plotQualityProfile(fnRs)   #visualizing the phred quality profiles of the forward reads
plotQualityProfile(fnRs[1:2]) 
plotQualityProfile(fnFs[1:2])
# forward reads were overall of good quality but reverse reads were much worse especially towards the end. Sequences of low quality will have to be trimmed.

#saving a funciton that will extract my sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) 
length(sample.names) # just checking if they are all there, in total we have 93 samples to analyze

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # the variable filFts will store my filtered/trimmed forward reads
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # the variable filtRs will store my filtered/trimmed reverse reads
names(filtFs) <- sample.names #forward reads names
names(filtRs) <- sample.names #reverse reads names
print(names(filtFs))
print(names(filtRs))
length(filtRs) #93 forward reads
length(filtFs) #93 reverse reads

# Because I have a full coverage of a sequence from my forward and reverse primers (16S v4 region is about 254bp) 
#and my sequences are around 250bp long,it is fine to be aggressive with trimming (I only need 20bp overlap for the alignment to work!) 
#truncating the sequences based on quality scores (truncLen= 200 for forward,150 for reverse)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen= 200, 150,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

#Checking if trimming was successful with quality profiles:
plotQualityProfile(filtRs[1:2]) 
plotQualityProfile(filtFs[1:2]) 

#Generating an error model of my data, required by DADA2 for inffering ASVs
errF <- learnErrors(filtFs, multithread=TRUE) # error model for filtered forward reads
plotErrors(errF, nominalQ=TRUE) # visualsing the error model
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errR, nominalQ=TRUE)
# For the plotted graphs, the black curve should follow the red line, if it doesn't it means that the quality of our sequences is not good enough. 
# Here the error graphs were fine, so I carried on to the next step.

#dereplication# this is an optional step in the work-flow, I chose to do it to speed up the computation
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

#Sample interference. I added pooling for increased sensitivity for singletons
dada_forward <- dada(derep_forward, err=errF, DETECT_SINGLETONS=TRUE, pool=TRUE, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=errR, DETECT_SINGLETONS=TRUE, pool=TRUE, multithread=TRUE)

#Merging our forward and reverse sequencing with a mergePairs function:
mergers <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, minOverlap=170, verbose=TRUE)
# I added minOverlap=120 as I trimmed the reverse primes at 150
length(mergers) # 93 elements = 93 samples
head(mergers[[1]]) #there were no mismatches 

# Contructing an amplicon sequence variant table (ASV) with makeSequenceTable
seqtab <- makeSequenceTable(mergers)
class(seqtab)
dim(seqtab)# In our 93 samples we have 2627 ASVs 
table(nchar(getSequences(seqtab))) # Some sequences are much shorter than amplified V4 area and are result of non-specific priming
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 248:252] # cutting out these short sequences and keeping only those between 250-256.
table(nchar(getSequences(seqtab2))) 

#Chimeras removal from the seqtab2 table
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim) # there are 1234 ASVs in 93 samples
sum(seqtab.nochim)/sum(seqtab2) # Here chimeras account for only about 5% of the merged reads, it shouldn't affect abundance
seqtab.nochim
#Assigning taxonomies:
taxa <- assignTaxonomy(seqtab.nochim,"C:/Users/Anna/Desktop/Daphnia_files/All files/filtered/silva_nr99_v138_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/Users/Anna/Desktop/Daphnia_files/All files/filtered/silva_species_assignment_v138.fa")
taxa.print <- taxa
rownames(taxa.print) <- NULL # Removing sequence rownames for display only

# giving our seq headers more manageable names (ASV1, ASV2, etc.) instead of sequences
asv_seq <- colnames(seqtab.nochim) #this variable will store the sequences
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

#A for loop that will replace the sequence with and ASV character and number
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste("ASV", i, sep="_")
}
colnames(seqtab.nochim) <- asv_headers
dim(seqtab.nochim)

# count table:
getwd()
asv_tab <- t(seqtab.nochim) #transposing rows and columns to match taxa table
write.table(asv_tab, "ASV_table.csv", sep=",", quote=F, col.names=NA) #saving as a CSV file

#taxa table
taxa
row.names(taxa) <- gsub(pattern=">", replacement="", x=asv_headers) #substituting row.names for asv.headers
write.table(taxa, "ASV_taxonomy.csv", sep = ",", quote=F, col.names=NA, row.names = TRUE) #saving as a CSV file

#Meta data table
meta_data
write.table(meta_data, "meta_data.csv", sep = ",", col.names=NA, row.names = TRUE) #saving as a CSV file

#Creating phyloseq elements:
ps <- phyloseq(otu_table(t(asv_tab), taxa_are_rows=FALSE), 
               sample_data(meta_data), 
               tax_table(taxa))

#Now that the samples are saved as phyloseq objects in ps I can investigate them a bit before the analysis:
length(sample_names(ps)) #there are 93 samples
rank_names(ps) #ranks start from Kingdom to Species
sample_variables(ps) # there is only one variable which is treatment

#Checking sequencing depth
sample_sums(ps)
write.table(sample_sums(ps), 'Sequencing depth raport.csv', sep = ',')
max(sample_sums(ps)) #Maximum number of sequences is 62411
min(sample_sums(ps)) #Minimum number of sequences is 6 

# The sequencing depth is not equal and samples will need to be normalised
#Normalisation by the sum
ASV
ASV = data.frame(otu_table(ps))
sums <-apply(ASV, 2, sum)
norm.asv <- ASV
for(i in 1:ncol(ASV))
{norm.asv[,i] <- ASV[,i]/sums[i]}

#Saving normalised otu.table for diversity analysis
write.table(norm.asv, 'Normalised ASV table.csv', sep = ',', col.names =NA, row.names =TRUE)






