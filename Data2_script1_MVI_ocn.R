#####################################
#Assignment 1
#Oivia Nieves, Hayden Dickerson, Maria Ingersoll

library(dada2); #packageVersion("dada2"); citation("dada2")
library(ShortRead); #packageVersion("ShortRead")
library(ggplot2); #packageVersion("ggplot2")
library(phyloseq); #packageVersion("phyloseq")

#Set path to unzipped, renamed fastq files

path <- "C:/Users/Olivia Nieves/Documents/BU undergrad/BI586/Apo-Sym-Oculina_2021_assign1/Oculina_control_16S_sym-apo" # CHANGE ME to the directory containing the fastq files after unzipping.
fns <- list.files(path)
#Let's make sure that all of our files are there
fns

################################
##### Trimming/Filtering #######
################################

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures reads are in same order
fnFs <- sort(list.files(path, pattern = "_R1_16S_final.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_16S_final.fastq", full.names = TRUE))
#fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files- these are old 454 data but most data are paired end

# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq; OTHERWISE MODIFY
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)
sample.names

#########Visualize Raw data

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs[c(1,2,3,4,5,6,7,8,9)])
plotQualityProfile(fnFs[c(10,11,12,13,14,15,16,17,18)])
plotQualityProfile(fnFs[c(19,20,21,22,23,24,25,26)])

#Recommend trimming where quality profile crashes - in this case, forward reads mostly fine up to 300
#For common ITS amplicon strategies with paired end reads, it is undesirable to truncate reads to a fixed length due to the large amount of length variation at that locus. That is OK, just leave out truncLen. Make sure you removed the forward and reverse primers from both the forward and reverse reads though! 

#Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
filtRs
filtFs

# Filter
# Truncating at 225, trimLeft=25, and changing matchID to TRUE per Nicola's code (this enforces matching between id-line seq identifiers of F and R reads)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen= 200, #end of single end reads = approx. 300 bp
                     maxN=0, #DADA does not allow Ns
                     maxEE=1, #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     trimLeft=25, #N nucleotides to remove from the start of each read: ITS2 primer = F 20bp
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchID=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out)
tail(out)

#A word on Expected Errors vs a blanket quality threshold
#Take a simple example: a read of length two with quality scores Q3 and Q40, corresponding to error probabilities P=0.5 and P=0.0001. The base with Q3 is much more likely to have an error than the base with Q40 (0.5/0.0001 = 5,000 times more likely), so we can ignore the Q40 base to a good approximation. Consider a large sample of reads with (Q3, Q40), then approximately half of them will have an error (because of the P=0.5 from the Q2 base). We express this by saying that the expected number of errors in a read with quality scores (Q3, Q40) is 0.5.
#As this example shows, low Q scores (high error probabilities) dominate expected errors, but this information is lost by averaging if low Qs appear in a read with mostly high Q scores. This explains why expected errors is a much better indicator of read accuracy than average Q.

################################
##### Learn Error Rates #######
################################
#DADA2 learns its error model from the data itself by alternating estimation of the error rates and the composition of the sample until they converge on a jointly consistent solution (this is similar to the E-M algorithm)
#As in many optimization problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Maximum cycles was set to 30, but Convergence was found after 4 rounds

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 

################################
##### Dereplicate reads #######
################################
#Dereplication combines all identical sequencing reads into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

################################
##### Infer Sequence Variants #######
################################

#Must change some of the DADA options b/c original program optomized for ribosomal data, not ITS - from github, "We currently recommend BAND_SIZE=32 for ITS data." leave as default for 16S/18S
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]
dadaRs[[1]]

#~############################~#
##### Merge paired reads #######
#~############################~#

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

summary((mergers[[1]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
head(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 244-264 window
#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing


################################
##### Remove chimeras #######
################################
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# Identified 125 bimeras out of 508 input sequences.
# [1] 13 383

sum(seqtab.nochim)/sum(seqtab)
# 0.5749017

#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)


write.csv(seqtab,file="Apo_Sym_Oc_seqtab.csv")
write.csv(seqtab.nochim,file="Apo_Sym_Oc_nochim.csv")
################################
##### Track Read Stats #######
################################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="ReadFilterStats_AllData_final.csv",row.names=TRUE,quote=FALSE)

#THIS IS ITS2 PART FROM ORIGONAL SCRIPT
# ################################
# ##### Assign Taxonomy #######
# ################################
# 
# #It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to classify sequence variants taxonomically. 
# #DADA2 provides a native implementation of the RDP's naive Bayesian classifier. The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, and outputs the taxonomic assignments with at least minBoot bootstrap confidence.
# #Here, I have supplied a modified version of the GeoSymbio ITS2 database listing more taxonomic info as phyloseq requires (Franklin et al. 2012)
# #For example: GeoSymbio data (taken from "all clades" at https://sites.google.com/site/geosymbio/downloads):
# #>A1.1
# #modified version for phyloseq looks like this instead:
# #>Symbiodinium; Clade A; A1.1
# 
# taxa <- assignTaxonomy(seqtab.nochim, "GeoSymbio_ITS2_LocalDatabase_verForPhyloseq.fasta", minBoot=5,multithread=TRUE,tryRC=TRUE,outputBootstraps=FALSE)
# #minboot should be higher
# #Obtain a csv file for the taxonomy so that it's easier to map the sequences for the heatmap.
# write.csv(taxa, file="taxa.csv",row.name=TRUE,quote=FALSE)
# unname(head(taxa, 30))
# unname(taxa)
# 
# #Now, save outputs so can come back to the analysis stage at a later point if desired
# saveRDS(seqtab.nochim, file="final_seqtab_nochim.rds")
# saveRDS(taxa, file="final_taxa_blastCorrected.rds")
# 
# #If you need to read in previously saved datafiles
# seqtab.nochim <- readRDS("final_seqtab_nochim.rds")
# taxa <- readRDS("final_taxa_blastCorrected.rds")
# head(taxa)
#~############################~#
##### Assign Taxonomy ##########
#~############################~#

# #Using package DECIPHER as an alternate to 'assignTaxonomy'
BiocManager::install("DECIPHER")
library(DECIPHER); packageVersion("DECIPHER")
# #citation("DECIPHER")
# 
# #http://DECIPHER.codes/Downloads.html. Download the SILVA SSU r132 (modified) file to follow along.
# Olivia: I downloaded r138
dna <- DNAStringSet(getSequences(seqtab.nochim)) 
# Create a DNAStringSet from the ASVs
load("C:/Users/Olivia Nieves/Documents/BU undergrad/BI586/Apo-Sym-Oculina_2021_assign1/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE, threshold=50) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
   taxa[startsWith(taxa, "unclassified_")] <- NA
   taxa
 }))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

#also doing other taxonomy method:
#Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Olivia Nieves/Documents/BU undergrad/BI586/Apo-Sym-Oculina_2021_assign1/silva_nr99_v138_train_set.fa",tryRC=TRUE)
unname(head(taxa))
taxa.plus <- addSpecies(taxa, "C:/Users/Olivia Nieves/Documents/BU undergrad/BI586/Apo-Sym-Oculina_2021_assign1/silva_species_assignment_v138.fa",tryRC=TRUE,verbose=TRUE)
# 18 out of 383 were assigned to the species level.
# Of which 10 had genera consistent with the input table.

saveRDS(taxa.plus, file="Oculina16s_taxaplus.rds")
saveRDS(taxa, file="Oculina16s_taxa.rds")
write.csv(taxa.plus, file="Oculina16s_taxaplus.csv")
write.csv(taxa, file="Oculina16s_taxa.csv")

saveRDS(seqtab.nochim, file="Oculina16s_seqtab.nochim.rds")
write.csv(seqtab.nochim, file="Oculina16s_seqtab.nochim.csv")
write.csv(seqtab.nochim, file="Oculina16s_seqtab.nochim_renamed.csv")


################################
##### handoff 2 phyloseq #######
################################

#import dataframe holding sample information
#have your samples in the same order as the seqtab file in the rows, variables as columns (put the metadata file in samdf)
samdf<-read.delim("Oculina_16S_metadata.txt")
head(samdf)
head(seqtab.nochim)
head(taxa)
rownames(samdf) <- samdf$SampleID
head(samdf)

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

### Remove mitochondria and chloroplasts ###
# First, get an idea of how many taxa you must remove
ps.mito <- subset_taxa(ps, (Family=="Mitochondria"))
ps.mito
    # phyloseq-class experiment-level object
    # otu_table()   OTU Table:         [ 51 taxa and 13 samples ]
    # sample_data() Sample Data:       [ 13 samples by 6 sample variables ]
    # tax_table()   Taxonomy Table:    [ 51 taxa by 6 taxonomic ranks ]

ps.chlor <- subset_taxa(ps, (Order=="Chloroplast"))
ps.chlor
    # phyloseq-class experiment-level object
    # otu_table()   OTU Table:         [ 23 taxa and 13 samples ]
    # sample_data() Sample Data:       [ 13 samples by 6 sample variables ]
    # tax_table()   Taxonomy Table:    [ 23 taxa by 6 taxonomic ranks ]

# Now sequentially remove the 51 Mitochondria taxa and the 23 chloroplast taxa
ps.nomito <- subset_taxa(ps, (Family!="Mitochondria") | is.na(Family))
ps.nomito
ps.clean <- subset_taxa(ps.nomito, (Order!="Chloroplast") | is.na(Order))
ps.clean
# End up with:
    # phyloseq-class experiment-level object
    # otu_table()   OTU Table:         [ 309 taxa and 13 samples ]
    # sample_data() Sample Data:       [ 13 samples by 6 sample variables ]
    # tax_table()   Taxonomy Table:    [ 309 taxa by 6 taxonomic ranks ]

#replace sequences with shorter names (correspondence table output below)
ids<-taxa_names(ps.clean)
ids <- paste0("sq",seq(1, length(colnames(seqtab.nochim))))
colnames(seqtab.nochim) <- ids

#Bar-plots
top90 <- names(sort(taxa_sums(ps.clean), decreasing=TRUE))[1:90]
ps.top90 <- transform_sample_counts(ps.clean, function(OTU) OTU/sum(OTU))
ps.top90 <- prune_taxa(top90, ps.top90)

plot_bar(ps.top90, x="Sample", fill="Class") 

#visusalize via counts rather than abundances:
plot_bar(ps.clean, x = "Sample", fill= "Class") #+ facet_wrap("tank")
#
#Obtain a csv file for the phyloseq data. Will give you the abundances for each sample and class. Useful for constructing the heatmap. Also, enables you to use ggplot, and construct more aesthetically pleasing plot.
psz <- psmelt(ps.top90)
write.csv(psz, file="Phyloseqoutputfinal.csv")
p <- ggplot(psz, aes(x=Sample, y=Abundance, fill=Class))
p + geom_bar(stat="identity", colour="black")
