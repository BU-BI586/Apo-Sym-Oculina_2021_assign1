#Nicola's Moorea 16S analysis
#Almost entirely based on DADA2 Pipeline 1.8 Walkthrough:
#https://benjjneb.github.io/dada2/tutorial.html
#with edits by Carly D. Kenkel and modifications for my data by Nicola Kriefall
#12/30/18

#~########################~#
##### PRE-PROCESSING #######
#~########################~#

#fastq files should have R1 & R2 designations for PE reads
#Also - some pre-trimming. Retain only PE reads that match amplicon primer. Remove reads containing Illumina sequencing adapters

#in Terminal home directory:
#following instructions of installing BBtools from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
#1. download BBMap package, sftp to installation directory
#2. untar: 
#tar -xvzf BBMap_(version).tar.gz
#3. test package:
#cd bbmap
#~/bin/bbmap/stats.sh in=~/bin/bbmap/resources/phix174_ill.ref.fa.gz

# my adaptors for 16S, which I saved as "adaptors.fasta"
# >forward
# AATGATACGGCGACCAC
# >forwardrc
# GTGGTCGCCGTATCATT
# >reverse
# CAAGCAGAAGACGGCATAC
# >reverserc
# GTATGCCGTCTTCTGCTTG

#primers for 16S: 
# >forward
# GTGYCAGCMGCCGCGGTA
# >reverse
# GGACTACHVGGGTWTCTAAT

##Still in terminal - making a sample list based on the first phrase before the underscore in the .fastq name
#ls *R1_001.fastq | cut -d '_' -f 1 > samples.list

##cuts off the extra words in the .fastq files
#for file in $(cat samples.list); do  mv ${file}_*R1*.fastq ${file}_R1.fastq; mv ${file}_*R2*.fastq ${file}_R2.fastq; done 

##gets rid of reads that still have the adaptor sequence, shouldn't be there, I didn't have any
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adaptors.fasta out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq; done &>bbduk_NoIll.log

##getting rid of first 4 bases (degenerate primers created them)
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq ftl=4 out1=${file}_R1_NoIll_No4N.fastq out2=${file}_R2_NoIll_No4N.fastq; done &>bbduk_No4N.log

##only keeping reads that start with the 16S primer
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq restrictleft=20 k=10 literal=GTGYCAGCMGCCGCGGTA,GGACTACHVGGGTWTCTAAT copyundefined=t outm1=${file}_R1_NoIll_No4N_16S.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_No4N_16S.fastq outu2=${file}_R2_check.fastq; done &>bbduk_16S.log
##higher k = more reads removed, but can't surpass k=20 or 21

##using cutadapt to remove primer
# for file in $(cat samples.list)
# do
# cutadapt -g GTGYCAGCMGCCGCGGTA -a ATTAGAWACCCVHGTAGTCC -G GGACTACHVGGGTWTCTAAT -A TACCGCGGCKGCTGRCAC -n 2 --discard-untrimmed -o ${file}_R1.fastq -p ${file}_R2.fastq ${file}_R1_NoIll_No4N_16S.fastq ${file}_R2_NoIll_No4N_16S.fastq
# done &> clip.log
##-g regular 5' forward primer 
##-G regular 5' reverse primer
##-o forward out
##-p reverse out
##-max-n 0 means 0 Ns allowed
##this overwrote my original renamed files 

# did sftp of *_R1.fastq & *_R2.fastq files to the folder to be used in dada2

#~########################~#
##### DADA2 BEGINS #########
#~########################~#

#installing/loading packages:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")
library(dada2); packageVersion("dada2")
#I have version 1.10.0 - tutorial says 1.8 but I think that's OK, can't find a version 1.10 walkthrough
library(ShortRead)
#packageVersion("ShortRead")
library(Biostrings)
#packageVersion("Biostrings")
path <- "~/Desktop/mr16s/files" # CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)
sample.names

#### check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[3]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[3]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[3]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[3]]))
#no primers - amazing

#### Visualizing raw data ####

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs.filtN[c(1,2,3,4)])
plotQualityProfile(fnFs.filtN[c(90,91,92,93)])
#looks mostly good up to 200

#Then look at quality profile of R2 reads
plotQualityProfile(fnRs.filtN[c(1,2,3,4)])
plotQualityProfile(fnRs.filtN[c(90,91,92,93)])
#180

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#changing a bit from default settings - maxEE=1 (1 max expected error, more conservative), truncating length at 200 bp for both forward & reverse [leaves ~50bp overlap], added "trimleft" to cut off primers [18 for forward, 20 for reverse]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(200,180), #leaves ~50bp overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     #trimLeft=c(18,20), #N nucleotides to remove from the start of each read
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)

#~############################~#
##### Learn Error Rates ########
#~############################~#

#setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 

#~############################~#
##### Dereplicate reads ########
#~############################~#
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#~###############################~#
##### Infer Sequence Variants #####
#~###############################~#

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

#~##################################~#
##### Construct sequence table #######
#~##################################~#
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 244-264 window

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(244,264)] #again, being fairly conservative wrt length

#~############################~#
##### Remove chimeras ##########
#~############################~#
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Identified 105 bimeras out of 2485 input sequences.

sum(seqtab.nochim)/sum(seqtab)
#0.9959795
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 

#~############################~#
##### Track Read Stats #########
#~############################~#

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="mr16s_readstats.csv",row.names=TRUE,quote=FALSE)

#~############################~#
##### Assign Taxonomy ##########
#~############################~#

# #Using package DECIPHER as an alternatie to 'assignTaxonomy'
# 
# #BiocManager::install("DECIPHER")
# library(DECIPHER); packageVersion("DECIPHER")
# #citation("DECIPHER")
# 
# #http://DECIPHER.codes/Downloads.html. Download the SILVA SSU r132 (modified) file to follow along.
# 
# dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# load("~/Downloads/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
# ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE, threshold=50) # use all processors
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid <- t(sapply(ids, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

#also doing other taxonomy method:
#Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/silva_nr_v132_train_set.fa.gz",tryRC=TRUE)
unname(head(taxa))
taxa.plus <- addSpecies(taxa, "~/Downloads/silva_species_assignment_v132.fa.gz",tryRC=TRUE,verbose=TRUE)
# 237 out of 2380 were assigned to the species level.
# Of which 214 had genera consistent with the input table.

saveRDS(taxa.plus, file="mr16s_taxaplus.rds")
saveRDS(taxa, file="mr16s_taxa.rds")
write.csv(taxa.plus, file="mr16s_taxaplus.csv")
write.csv(taxa, file="mr16s_taxa.csv")

saveRDS(seqtab.nochim, file="mr16s_seqtab.nochim.rds")
write.csv(seqtab.nochim, file="mr16s_seqtab.nochim.csv")
write.csv(seqtab.nochim, file="mr16s_seqtab.nochim_renamed.csv")

#### Read in previously saved datafiles ####
setwd("~/moorea_holobiont/mr_16S")
seqtab.nochim <- readRDS("mr16s_seqtab.nochim.rds")
taxa <- readRDS("mr16s_taxa.rds")
taxa.plus <- readRDS("mr16s_taxaplus.rds")

#~############################~#
##### handoff 2 phyloseq #######
#~############################~#

#BiocManager::install("phyloseq")
library('phyloseq')
library('ggplot2')
library('Rmisc')
library(cowplot)

#import dataframe holding sample information
samdf<-read.csv("mr16s_sampledata.csv")
head(samdf)
rownames(samdf) <- samdf$id

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.plus))

ps

#first look at data
ps_glom <- tax_glom(ps, "Family")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, x="site_zone", fill="Family")+
  theme(legend.position="none")

#phyloseq object with shorter names - doing this one instead of one above
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
#making output fasta file for lulu step & maybe other things, before giving new ids to sequences
#path='~/moorea_holobiont/mr_16s/mr16s.fasta'
#uniquesToFasta(seqtab.nochim, path, ids = ids, mode = "w", width = 20000)

colnames(seqtab.nochim)<-ids
taxa2 <- cbind(taxa.plus, rownames(taxa.plus)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

#phyloseq object with new taxa ids
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa2))

ps #2380 taxa

#### remove mitochondria, chloroplasts, non-bacteria #### 
ps.mito <- subset_taxa(ps, (Family=="Mitochondria"))
ps.mito #52 taxa to remove
ps.chlor <- subset_taxa(ps, (Order=="Chloroplast"))
ps.chlor #75 taxa to remove
ps.notbact <- subset_taxa(ps, (Kingdom!="Bacteria") | is.na(Kingdom))
ps.notbact #85 taxa to remove

ps.nomito <- subset_taxa(ps, (Family!="Mitochondria") | is.na(Family))
ps.nomito #2328 taxa
ps.nochlor <- subset_taxa(ps.nomito, (Order!="Chloroplast") | is.na(Order))
ps.nochlor #2253 taxa
ps.clean <- subset_taxa(ps.nochlor, (Kingdom=="Bacteria"))
ps.clean #2168 taxa

#just archaea
ps.arch <- subset_taxa(ps.nomito, (Kingdom=="Archaea"))
ps.arch #17 taxa

seqtab.clean <- data.frame(otu_table(ps.clean))
write.csv(seqtab.clean,file="seqtab.cleaned.csv")
seqtab.clean <- read.csv("seqtab.cleaned.csv",row.names=1)

#### rarefy #####
library(vegan)

rarecurve(seqtab.clean,step=100,label=FALSE) #after removing contaminats

total <- rowSums(seqtab.clean)
subset(total, total <12000)
#9 samples
#B5 & F9 identified by MCMC.OTU below as being too low

row.names.remove <- c("B5","F9") #being less strict
row.names.remove <- c("A8","B5","B8","C2","D4","E8","F9","H4","H5")
seqtab.less <- seqtab.clean[!(row.names(seqtab.clean) %in% row.names.remove),]
samdf.rare <- samdf[!(row.names(samdf) %in% row.names.remove), ]
#84 samples left being strict
#91 samples left being less strict

seqtab.rare <- rrarefy(seqtab.less,sample=12000)
rarecurve(seqtab.rare,step=100,label=FALSE)

#phyloseq object but rarefied
ps.rare <- phyloseq(otu_table(seqtab.rare, taxa_are_rows=FALSE), 
                    sample_data(samdf.rare), 
                    tax_table(taxa2))
ps.rare #2253 taxa first time through
#2168 another time

#removing missing taxa - lost after rarefying
ps.rare <- prune_taxa(taxa_sums(ps.rare) > 0, ps.rare)
ps.rare #2101 taxa
#2018 taxa round 2

seqtab.rare <- data.frame(otu_table(ps.rare))

#saving
write.csv(seqtab.rare, file="mr16s_seqtab.rare_12k_rd2.csv")
write.csv(samdf.rare, file="mr16s_samdf.rare_12k.csv")
#re-reading
samdf.rare <- read.csv("mr16s_samdf.rare_12k.csv",row.names=1)
seqtab.rare <- read.csv("mr16s_seqtab.rare_12k_rd2.csv",row.names=1)

#### alpha diversity ####
#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).
plot_richness(ps.rare, x="site", measures=c("Shannon", "Simpson"), color="zone") + theme_bw()

df <- data.frame(estimate_richness(ps.rare, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))
df <- data.frame(estimate_richness(ps.clean, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))
df

df$id <- rownames(df)
df.div <- merge(df,samdf,by="id") #add sample data
df.div <- merge(df,samdf.rare,by="id") #add sample data

write.csv(df.div,file="mr16s_diversity_rare12k.csv") #saving
df.div <- read.csv("mr16s_diversity_rare12k.csv",row.names=1,header=TRUE) #reading back in 

quartz()
gg.site.sha <- ggplot(df.div,aes(x=site,y=Shannon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5)+
  ylab("Shannon Diversity")+
  xlab("Site")+
  theme_cowplot()

gg.site.sim <- ggplot(df.div,aes(x=site,y=InvSimpson))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5)+
  ylab("Inv. Simpson diversity")+
  xlab("Site")+
  theme_cowplot()

gg.site.obs <- ggplot(df.div,aes(x=site,y=Observed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5)+
  ylab("OTU richness")+
  xlab("Site")+
  theme_cowplot()

ggarrange(gg.site.obs,gg.site.sha,gg.site.sim,nrow=1)

gg.sh <- ggplot(df.div, aes(x=zone, y=Shannon,color=zone,shape=zone))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Reef zone")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("BR","FR"))+
  scale_x_discrete(labels=c("BR","FR"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  #guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  theme(legend.position="none")+
  geom_jitter(alpha=0.5)+
  facet_wrap(~site)
gg.sh

gg.si <- ggplot(df.div, aes(x=zone, y=InvSimpson,color=zone,shape=zone))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Reef zone")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  #guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  theme(legend.position="none")+
  geom_jitter(alpha=0.5)+
  facet_wrap(~site)+
  scale_x_discrete(labels=c("BR","FR"))
gg.si

gg.obs <- ggplot(df.div, aes(x=zone, y=Observed,color=zone,shape=zone))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Reef zone")+
  ylab("OTU richness")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  #guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  theme(legend.position="none")+
  geom_jitter(alpha=0.5)+
  facet_wrap(~site)+
  scale_x_discrete(labels=c("BR","FR"))
gg.obs

ggarrange(gg.obs,gg.sh,gg.si,nrow=1)

#stats
library(car)

shapiro.test(df.div$Observed) #nope
df.div$obs.log <- log(df.div$Observed)
shapiro.test(df.div$obs.log) #yessss
leveneTest(df.div$obs.log~site*zone,data=df.div) #fine

a.div <- aov(obs.log~site,data=df.div)
summary(a.div)
TukeyHSD(a.div) #Tahiti different from the other two 

shapiro.test(df.div$Shannon) #fine
leveneTest(df.div$Shannon~site*zone,data=df.div) #fine

a.div <- aov(Shannon~site*zone,data=df.div)
summary(a.div)
TukeyHSD(a.div) #nothing

shapiro.test(df.div$InvSimpson) #not normal
df.div$si.log <- log(df.div$InvSimpson) 
shapiro.test(df.div$si.log) #NORMAL
leveneTest(df.div$si.log~site*zone,data=df.div) #fine 

a.div <- aov(si.log~site*zone,data=df.div)
summary(a.div)
TukeyHSD(a.div) #nothing

#### diversity + size stats ####
size <- read.csv("~/moorea_holobiont/mr_ITS2/mr_size.csv",header=TRUE,row.names=1)

row.names.remove <- c(109)
size2 <- size[!(row.names(size) %in% row.names.remove),]

row.names(df.div) <- df.div$sample
row.names(size2) <- size2$coral_id

df.div.size <- merge(df.div,size2,by=0)

lm.sh.size <- lm(size~Shannon*zone.y*site.y,data=df.div.size)
lm.sh.size <- lm(size~Shannon,data=df.div.size)
summary(lm.sh.size) #not sig, with or without zone & site

lm.si.size <- lm(size~si.log*zone.y*site.y,data=df.div.size)
lm.si.size <- lm(size~si.log,data=df.div.size)
summary(lm.si.size) #not sig, with or without zone & site

lm.ob.size <- lm(size~obs.log*zone.y*site.y,data=df.div.size)
lm.ob.size <- lm(size~obs.log,data=df.div.size)
summary(lm.ob.size)

#### trim underrepresented ASVs ####
library(MCMC.OTU)

seq.formcmc <- read.csv("seqtab.cleaned_formcmc.csv")
#seq.formcmc <- read.csv("mr16s_seqtab.rare_12k_formcmc.csv")
seq.formcmc <- read.csv("mr16s_seqtab.rare_12k_rd2_formcmc.csv")

#regular
seq.trim <- purgeOutliers(seq.formcmc,count.columns=3:2170,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
#2 bad samples - B5 & F9
#207 ASVs passing cutoff for reads
#195 ASVs show up in 2% of samples 
#remove bad samples from sample data frame
row.names.remove <- c("B5","F9")
samdf.trim <- samdf[!(row.names(samdf) %in% row.names.remove), ]
#rarefied
seq.trim <- purgeOutliers(seq.formcmc,count.columns=3:2020,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
#no bad samples
#223 ASVs passing cutoffs - rd1 
#209 show up in 2% of samples - rd1
#221 ASVs passing cutoffs - rd2
#207 show up in 2% of samples - rd2

#rename rows
rownames(seq.trim) <- seq.trim$sample
#remove sample info
#seq.trim <- seq.trim[,3:211] #rarefied
seq.trim <- seq.trim[,3:209] #rarefied rd2
seq.trim <- seq.trim[,3:197] #regular

write.csv(seq.trim,file="seq.trim.csv")
seq.trim <- read.csv("seq.trim.csv",row.names=1)

write.csv(seq.trim,file="seq.rare12k.trim_rd2.csv")
seq.trim <- read.csv("seq.rare12k.trim_rd2.csv",row.names=1)

#remake phyloseq objects
ps.trim <- phyloseq(otu_table(seq.trim, taxa_are_rows=FALSE), 
                    sample_data(samdf.trim), 
                    tax_table(taxa2))
ps.trim #195 taxa

ps.rare.trim <- phyloseq(otu_table(seq.trim, taxa_are_rows=FALSE), 
                         sample_data(samdf), 
                         tax_table(taxa2))
ps.rare.trim #207 taxa rd2

#### pcoa plots ####
library(stats)
library(MCMC.OTU)

#transform to relative abundance
ps.trim.rel <- transform_sample_counts(ps.trim, function(x) x / sum(x))
ord <- ordinate(ps.trim.rel, "PCoA", "bray")
plot_ordination(ps.trim.rel, ord,color="site", shape="site")+
  stat_ellipse()

#saving
seq.trim.rel <- data.frame(otu_table(ps.trim.rel))

ord <- ordinate(ps.rare.trim, "PCoA", "bray")
quartz()
gg.pcoa.site.rare <- plot_ordination(ps.rare.trim, ord,color="site", shape="site")+
  stat_ellipse()+
  theme_cowplot()+
  scale_color_manual(name="Site",values=c("darkslategray3","darkslategray4","#000004"))+
  scale_shape_manual(name="Site",values=c(8,4,9))+
  xlab("Axis 1 (32.3%)")+
  ylab("Axis 2 (20.1%)")+
  annotate(geom="text", x=0.35, y=0.65, label="p < 0.01**",size=4)+
  ggtitle("Rarefied")
gg.pcoa.site.rare

ord.rel <- ordinate(ps.trim.rel, "PCoA", "bray")
gg.pcoa.site <- plot_ordination(ps.trim.rel, ord.rel,color="site", shape="site")+
  stat_ellipse()+
  theme_cowplot()+
  scale_color_manual(name="Site",values=c("darkslategray3","darkslategray4","#000004"))+
  scale_shape_manual(name="Site",values=c(8,4,9))+
  xlab("Axis 1 (32.4%)")+
  ylab("Axis 2 (21%)")+
  annotate(geom="text", x=0.35, y=0.65, label="p < 0.01**",size=4)+
  ggtitle("Relative abundance")
gg.pcoa.site

quartz()
ggarrange(gg.pcoa.site.rare,gg.pcoa.site,labels="AUTO",common.legend=TRUE,legend="right")

#other color options
#scale_color_manual(name="Site",values=c("darkturquoise","darkslategray4","#000004"))+

#now by site
ps.mnw <- subset_samples(ps.trim.rel,site=="MNW")
#ps.mnw <- subset_samples(ps.rare.trim,site=="MNW")
ord.mnw <- ordinate(ps.mnw, "PCoA", "bray")
gg.mnw <- plot_ordination(ps.mnw, ord.mnw,color="zone", shape="zone")+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  ggtitle("Moorea NW")+
  #xlab("Axis 1 (33.2%)")+#rarefied
  #ylab("Axis 2 (22.5%)")+#rarefied
  #xlab("Axis 1 (26.2%)")+#non-rarefied
  #ylab("Axis 2 (21.3%)")+#non-rarefied
  theme(axis.text=element_text(size=10))
gg.mnw

ps.mse <- subset_samples(ps.trim.rel,site=="MSE")
#ps.mse <- subset_samples(ps.trim,site=="MSE")
ord.mse <- ordinate(ps.mse, "PCoA", "bray")
gg.mse <- plot_ordination(ps.mse, ord.mse,color="zone", shape="zone")+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  ggtitle("Moorea SE")+
  annotate(geom="text", x=-0.35, y=0.55, label="p < 0.01**",size=4)+#change for non-rarefied
  #xlab("Axis 1 (33.1%)")+ #rarefied
  #ylab("Axis 2 (22.2%)")+ #rarefied
  xlab("Axis 1 (27.9%)")+#non-rarefied
  ylab("Axis 2 (19%)")+#non-rarefied
  theme(axis.text=element_text(size=10))
gg.mse

ps.tnw <- subset_samples(ps.trim.rel,site=="TNW")
#ps.tnw <- subset_samples(ps.trim,site=="TNW")
ord.tnw <- ordinate(ps.tnw, "PCoA", "bray")
gg.tnw <- plot_ordination(ps.tnw, ord.tnw,color="zone", shape="zone")+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  ggtitle("Tahiti NW")+
  #annotate(geom="text", x=-0.4, y=0.6, label="p < 0.01**",size=4)+ #rarefied
  annotate(geom="text", x=-0.25, y=0.6, label="p < 0.01**",size=4)+ #not rarefied
  #xlab("Axis 1 (42.6%)")+ #rarefied
  #ylab("Axis 2 (27.9%)")+ #rarefied
  xlab("Axis 1 (34.3%)")+#non-rarefied
  ylab("Axis 2 (26.5%)")+#non-rarefied
  theme(axis.text=element_text(size=10))
gg.tnw

#### pcoa plot ####
library(ggpubr)

quartz()
ggarrange(gg.mnw,gg.mse,gg.tnw,nrow=1,common.legend=TRUE,legend="right",labels="AUTO")

### carly's method
df.seq <- as.data.frame(seq.trim)
all.log=logLin(data=df.seq)

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
all.dist=vegdist(all.log,method="bray")
all.pcoa=pcoa(all.dist)

# plotting:
scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
scorez$id <- rownames(scorez)
pcoa.all <- merge(scorez,samdf.rare,by='id')
pcoa.all <- merge(scorez,samdf.trim,by='id')

quartz()
ggplot(pcoa.all,aes(x=Axis.1,y=Axis.2,color=zone,shape=zone))+
  geom_point()+
  stat_ellipse()+
  facet_wrap(~site)+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Back reef","Fore reef"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  xlab("Axis 1 (54.14%)")+
  ylab("Axis 2 (33.58%)")
#looks the same as before rarefying

#### Stats ####
library(vegan)
#help on adonis here:
#https://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html#more

#not rarefied
#rel abundance or not
ps.trim.rel <- transform_sample_counts(ps.trim, function(x) x / sum(x))
seq.trim.rel <- data.frame(otu_table(ps.trim.rel))

dist.seqtab <- vegdist(seq.trim.rel)
dist.seqtab <- vegdist(seq.trim)
bet.all <- betadisper(dist.seqtab,samdf.trim$zone)
bet.all <- betadisper(dist.seqtab,samdf.trim$site)
anova(bet.all)
#not sig - rel or not
plot(bet.all)
adonis(seq.trim ~ site, data=samdf.trim, permutations=999)
adonis(seq.trim ~ site*zone, data=samdf.trim, permutations=999)
adonis(seq.trim ~ zone, strata=samdf.trim$site,data=samdf.trim, permutations=999)
adonis(seq.trim.rel ~ zone, strata=samdf.trim$site, data=samdf.trim, permutations=999)
adonis(seq.trim.rel ~ site*zone, data=samdf.trim, permutations=999)
adonis(seq.trim.rel ~ site, data=samdf.trim, permutations=999)

#site pairs
samdf.mnw.tnw <- subset(samdf.trim,site==c("MNW","TNW"))
sams.mnw.tnw <- c(rownames(samdf.mnw.tnw))
seq.mnw.tnw <- seq.trim.rel[(row.names(seq.trim.rel) %in% sams.mnw.tnw),]
adonis(seq.mnw.tnw ~ site*zone, data=samdf.mnw.tnw, permutations=999)
#MNW & TNW not different
samdf.mse.tnw <- subset(samdf.trim,site==c("MSE","TNW"))
sams.mse.tnw <- c(rownames(samdf.mse.tnw))
seq.mse.tnw <- seq.trim.rel[(row.names(seq.trim.rel) %in% sams.mse.tnw),]
adonis(seq.mse.tnw ~ site, data=samdf.mse.tnw, permutations=999)
#TNW & MSE not sig different
samdf.mnw.mse <- subset(samdf.trim,site==c("MNW","MSE"))
sams.mnw.mse <- c(rownames(samdf.mnw.mse))
seq.mnw.mse <- seq.trim.rel[(row.names(seq.trim.rel) %in% sams.mnw.mse),]
adonis(seq.mnw.mse ~ site, data=samdf.mnw.mse, permutations=999)
#MSE & MNW not sig different

#easier method for site pairwise:
samdf.trim.rare <- data.frame(ps.rare.trim@sam_data)
row.names(samdf.trim.rare) == row.names(seq.trim)

#install.packages("remotes")
#remotes::install_github("Jtrachsel/funfuns")
library("funfuns")
pairwise.adonis(seq.trim, factors = samdf.trim.rare$site, permutations = 999)
# pairs  F.Model         R2 p.value p.adjusted
# 1 MNW vs MSE 2.909039 0.05023464   0.010      0.010
# 2 MNW vs TNW 1.592082 0.02970741   0.141      0.141
# 3 MSE vs TNW 3.450676 0.05903570   0.006      0.006

#by site
samdf.mnw <- subset(samdf.trim,site=="MNW")
sam.mnw <- rownames(samdf.mnw)
seq.mnw <- seq.trim.rel[(row.names(seq.trim.rel) %in% sam.mnw),]

dist.mnw <- vegdist(seq.mnw)
bet.mnw <- betadisper(dist.mnw,samdf.mnw$zone,bias.adjust = TRUE,type="median")
anova(bet.mnw)
permutest(bet.mnw, pairwise = FALSE, permutations = 99)
plot(bet.mnw)
#not sig
adonis(seq.mnw ~ zone, strata=samdf.mnw$site, data=samdf.mnw, permutations=999)
#0.299

samdf.mse <- subset(samdf.trim,site=="MSE")
sam.mse <- rownames(samdf.mse)
seq.mse <- seq.trim[(row.names(seq.trim) %in% sam.mse),]

dist.mse <- vegdist(seq.mse)
bet.mse <- betadisper(dist.mse,samdf.mse$zone,bias.adjust = TRUE,type="median")
anova(bet.mse)
permutest(bet.mse, pairwise = FALSE, permutations = 99)
plot(bet.mse)
#not sig
adonis(seq.mse ~ zone, strata=samdf.mse$site, data=samdf.mse, permutations=999)
#sig! 0.004 **

samdf.tnw <- subset(samdf.trim,site=="TNW")
sam.tnw <- rownames(samdf.tnw)
seq.tnw <- seq.trim[(row.names(seq.trim) %in% sam.tnw),]

dist.tnw <- vegdist(seq.tnw)
bet.tnw <- betadisper(dist.tnw,samdf.tnw$zone,bias.adjust = TRUE,type="median")
anova(bet.tnw)
permutest(bet.tnw, pairwise = FALSE, permutations = 99)
plot(bet.tnw)
#not sig
adonis(seq.tnw ~ zone, strata=samdf.tnw$site, data=samdf.tnw, permutations=999)
#sig! 0.001 **

#### stats - rarefied ####
library(vegan)
#by site
dist.seqtab <- vegdist(seq.trim)
bet.all <- betadisper(dist.seqtab,samdf.rare$site)
anova(bet.all) #not sig 
plot(bet.all)
adonis(seq.trim ~ site, data=samdf.rare, permutations=999)

#site pairs
samdf.mnw.tnw <- subset(samdf.rare,site==c("MNW","TNW"))
sams.mnw.tnw <- c(rownames(samdf.mnw.tnw))
seq.mnw.tnw <- seq.trim[(row.names(seq.trim) %in% sams.mnw.tnw),]
adonis(seq.mnw.tnw ~ site*zone, data=samdf.mnw.tnw, permutations=999)
#MNW & TNW not different
samdf.mse.tnw <- subset(samdf.rare,site==c("MSE","TNW"))
sams.mse.tnw <- c(rownames(samdf.mse.tnw))
seq.mse.tnw <- seq.trim[(row.names(seq.trim) %in% sams.mse.tnw),]
adonis(seq.mse.tnw ~ site, data=samdf.mse.tnw, permutations=999)
#TNW & MSE sig different
samdf.mnw.mse <- subset(samdf.rare,site==c("MNW","MSE"))
sams.mnw.mse <- c(rownames(samdf.mnw.mse))
seq.mnw.mse <- seq.trim[(row.names(seq.trim) %in% sams.mnw.mse),]
adonis(seq.mnw.mse ~ site*zone, data=samdf.mnw.mse, permutations=999)
#MSE & MNW not different

#by zone
dist.seqtab <- vegdist(seq.trim)
bet.all <- betadisper(dist.seqtab,samdf.rare$zone)
anova(bet.all)
#not sig
plot(bet.all)
adonis(seq.trim ~ zone, strata=samdf.rare$site, data=samdf.rare, permutations=999)
#not sig

#by site
samdf.mnw <- subset(samdf.rare,site=="MNW")
sam.mnw <- rownames(samdf.mnw)
seq.mnw <- seq.trim[(row.names(seq.trim) %in% sam.mnw),]

dist.mnw <- vegdist(seq.mnw)
bet.mnw <- betadisper(dist.mnw,samdf.mnw$zone,bias.adjust = TRUE,type="median")
bet.mnw
anova(bet.mnw)
permutest(bet.mnw, pairwise = FALSE, permutations = 99)
plot(bet.mnw)
#not sig
adonis(seq.mnw ~ zone, data=samdf.mnw, permutations=999)
#0.384

samdf.mse <- subset(samdf.rare,site=="MSE")
sam.mse <- rownames(samdf.mse)
seq.mse <- seq.trim[(row.names(seq.trim) %in% sam.mse),]

dist.mse <- vegdist(seq.mse)
bet.mse <- betadisper(dist.mse,samdf.mse$zone,bias.adjust = TRUE,type="median")
anova(bet.mse)
permutest(bet.mse, pairwise = FALSE, permutations = 99)
plot(bet.mse)
#not sig
adonis(seq.mse ~ zone, data=samdf.mse, permutations=999)
#sig! 0.001 **

samdf.tnw <- subset(samdf.rare,site=="TNW")
sam.tnw <- rownames(samdf.tnw)
seq.tnw <- seq.trim[(row.names(seq.trim) %in% sam.tnw),]

dist.tnw <- vegdist(seq.tnw)
bet.tnw <- betadisper(dist.tnw,samdf.tnw$zone,bias.adjust = TRUE,type="median")
anova(bet.tnw)
permutest(bet.tnw, pairwise = FALSE, permutations = 99)
plot(bet.tnw)
#not sig
adonis(seq.tnw ~ zone, data=samdf.tnw, permutations=999)
#sig! 0.003 **

#### adonis stats with size ####
library('vegan')
size <- read.csv("~/moorea_holobiont/mr_ITS2/mr_size.csv",header=TRUE,row.names=1)

row.names.remove <- c(109)
size2 <- size[!(row.names(size) %in% row.names.remove),]

row.names(samdf) <- samdf$sample
row.names(size2) <- size2$coral_id

samdf.size2 <- merge(samdf,size2,by=0)

row.names(samdf.size2) <- samdf.size2$id

size.rows <- row.names(samdf.size2)
size3 <- seq.trim[(row.names(seq.trim) %in% size.rows),]
# size.rows2 <- c(row.names(size3))
# samdf.size2 <- samdf.size[(row.names(samdf.size) %in% size.rows2),]

#samdf.size.sorted <- samdf.size2[nrow(samdf.size2):1, ]
#size3.sorted <- size3[nrow(size3):1, ]

#tahiti 
samdf.size.tnw <- subset(samdf.size2,place=="T")
names <- c(samdf.size.tnw$id)
counts.size.tnw <- size3[(rownames(size3) %in% names),]
samdf.size.tnw2 <- samdf.size.tnw[(rownames(samdf.size.tnw) %in% row.names(counts.size.tnw)),]

samdf.size.tnw.sorted <- samdf.size.tnw2[sort(rownames(samdf.size.tnw2)),]
counts.size.tnw.sorted <- counts.size.tnw[sort(rownames(counts.size.tnw)),]

row.names(samdf.size.tnw.sorted) == row.names(counts.size.tnw.sorted)

dist.tnw <- vegdist(counts.size.tnw.sorted)
bet.tnw <- betadisper(dist.tnw,samdf.size.tnw.sorted$zone.y)
anova(bet.tnw)
permutest(bet.tnw, pairwise = FALSE, permutations = 99)
plot(bet.tnw) #not sig

adonis(counts.size.tnw.sorted ~ zone.y*size, data=samdf.size.tnw.sorted, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# zone.y       1   0.44388 0.44388  4.4881 0.14763  0.004 **
#   size         1   0.24380 0.24380  2.4650 0.08108  0.097 . 
# zone.y:size  1   0.04432 0.04432  0.4481 0.01474  0.784   
# Residuals   23   2.27475 0.09890         0.75655          
# Total       26   3.00675                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#moorea nw 
samdf.size.mnw <- subset(samdf.size2,place=="MNW")
names <- c(samdf.size.mnw$id)
counts.size.mnw <- size3[(rownames(size3) %in% names),]
samdf.size.mnw2 <- samdf.size.mnw[(rownames(samdf.size.mnw) %in% row.names(counts.size.mnw)),]

samdf.size.mnw.sorted <- samdf.size.mnw2[sort(rownames(samdf.size.mnw2)),]
counts.size.mnw.sorted <- counts.size.mnw[sort(rownames(counts.size.mnw)),]

row.names(samdf.size.mnw.sorted) == row.names(counts.size.mnw.sorted)

dist.mnw <- vegdist(counts.size.mnw.sorted)
bet.mnw <- betadisper(dist.mnw,samdf.size.mnw.sorted$zone.y)
anova(bet.mnw)
permutest(bet.mnw, pairwise = FALSE, permutations = 99)
plot(bet.mnw) #not sig

adonis(counts.size.mnw.sorted ~ zone.y*size, data=samdf.size.mnw.sorted, permutations=999)
#not sig

#### rename ASVs ####
library(rlang)
library(stringr)

tax <- as.data.frame(ps.rare.trim@tax_table@.Data)

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "D_0__",""),
                        Phylum = str_replace(tax[,2], "D_1__",""),
                        Class = str_replace(tax[,3], "D_2__",""),
                        Order = str_replace(tax[,4], "D_3__",""),
                        Family = str_replace(tax[,5], "D_4__",""),
                        Genus = str_replace(tax[,6], "D_5__",""),
                        Species = str_replace(tax[,7], "D_6__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
####### Fill holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

tax_table(ps.rare.trim) <- as.matrix(tax.clean)

#### bar plot summed by reef zone & site ####
library(dplyr)
ps.all <- transform_sample_counts(ps.rare.trim, function(OTU) OTU/sum(OTU))
pa <- psmelt(ps.all)
tb <- psmelt(ps.all)%>%
  filter(!is.na(Abundance))%>%
  group_by(site,zone,site_zone,Class,OTU)%>%
  summarize_at("Abundance",mean)

tb$zone <- gsub("out","FR",tb$zone)
tb$zone <- gsub("in","BR",tb$zone)

tb$site <- gsub("MNW","Mo'orea NW",tb$site)
tb$site <- gsub("MSE","Mo'orea SE",tb$site)
tb$site <- gsub("TNW","Tahiti NW",tb$site)

quartz()
ggplot(tb,aes(x=zone,y=Abundance,fill=Class))+
  geom_bar(stat="identity")+
  theme_cowplot()+
  #theme(legend.position="none")+
  xlab('Reef zone')+  
  facet_wrap(~site)

#### core v accessory microbiome ####
#BiocManager::install("microbiome")
#remotes::install_github("r-lib/rlang")
library(microbiome)

pseq.core <- core(ps.rare.trim, detection = 0, prevalence = .7)
pseq.core #10 taxa

ps_glom <- tax_glom(pseq.core, "Genus")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))

plot_bar(ps2, fill="Genus")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  scale_fill_brewer(palette="BrBG")

#by site
ps.core.mnw <- subset_samples(pseq.core,site=="MNW")
ps.core.mse <- subset_samples(pseq.core,site=="MSE")
ps.core.tnw <- subset_samples(pseq.core,site=="TNW")

#alternatively,
library(dplyr)

ps.all <- transform_sample_counts(pseq.core, function(OTU) OTU/sum(OTU))
pa <- psmelt(ps.all)
tb <- psmelt(ps.all)%>%
  #filter(!is.na(Abundance))%>%
  group_by(site_zone,Genus)#%>%
summarize_at("Abundance",mean)

ggplot(tb,aes(x=site_zone,y=Abundance,fill=Genus))+
  geom_bar(stat="identity")+
  theme_cowplot()+
  #theme(legend.position="none")+
  xlab('Reef zone')#+  
#facet_wrap(~site)

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

#alternatively alternatively
ps.core.rel <- transform_sample_counts(pseq.core, function(x) x / sum(x))
psm.all <- psmelt(ps.core.rel)
#melt.mnw$Abundance[is.na(melt.mnw$Abundance)] <- 0
psm.all.sum <- summarySE(psm.all,measurevar="Abundance",groupvars = c("Genus","site_zone"))
ggplot(psm.all.sum,aes(x=zone,y=Abundance,fill=Genus))+
  geom_bar(stat="identity")+
  facet_wrap(~site)
ggplot(psm.all.sum,aes(x=site_zone,y=Abundance,fill=Genus))+
  geom_bar(stat="identity")#+
facet_wrap(~site)

out <- data.frame(pseq.core@otu_table)
# calculating core abundances #
core.sqs <- tax_table(pseq.core)
core.sqs.ids <- row.names(core.sqs)
core.sqs.ids

ps.rare.trim.rel <- transform_sample_counts(ps.rare.trim, function(x) x / sum(x))
seq.rare.rel <- data.frame(otu_table(ps.rare.trim.rel))
tax.core <- tax_table(ps.rare.trim.rel)

library('tidyverse')
seq.core <- seq.rare.rel %>% select(c(core.sqs.ids))
colMeans(seq.core)

#just checking what they were in the non-core table:
colMeans(seq.rare.rel)

ps.core.rel <- phyloseq(otu_table(seq.core, taxa_are_rows=FALSE), 
                        sample_data(samdf.rare))

tax_table(ps.core.rel) <- as.matrix(tax.clean)
ps.core.rel #10 taxa

plot_bar(ps.core.rel, x="zone",fill="Genus")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  facet_grid(site~Genus)

#DESEQ to find differentially abundant cores 
library(DESeq2)

ds.all = phyloseq_to_deseq2(pseq.core, ~ site*zone)
dds.all <- estimateSizeFactors(ds.all,type="poscounts")
stat.all = DESeq(dds.all, test="Wald", fitType="parametric")
res = results(stat.all, cooksCutoff = FALSE)
alpha = 0.05
sigtab.all = res[which(res$padj < alpha), ]
sigtab.all = cbind(as(sigtab.all, "data.frame"), as(tax_table(ps.rare.trim)[rownames(sigtab.all), ], "matrix"))
sigtab.all
dim(sigtab.all) #sq 3 & 5

#stats/plotting?
library(vegan)

#MOOREA NW CORE
mnw.core.otu <- data.frame(ps.core.mnw@otu_table)
mnw.core.sam <- data.frame(ps.core.mnw@sam_data)

dist.all <- vegdist(mnw.core.otu)
row.names(mnw.core.otu) == row.names(mnw.core.sam)
bet.all <- betadisper(dist.all,mnw.core.sam$zone,bias.adjust = TRUE,type="median")
anova(bet.all) #nope
permutest(bet.all, pairwise = FALSE, permutations = 99)
plot(bet.all) #nope
adonis(mnw.core.otu ~ zone, data=mnw.core.sam, permutations=999)
#nope nope

#MOOREA SE CORE
mse.core.otu <- data.frame(ps.core.mse@otu_table)
mse.core.sam <- data.frame(ps.core.mse@sam_data)

dist.all <- vegdist(mse.core.otu)
row.names(mse.core.otu) == row.names(mse.core.sam)
bet.all <- betadisper(dist.all,mse.core.sam$zone,bias.adjust = TRUE,type="median")
anova(bet.all) #nope
permutest(bet.all, pairwise = FALSE, permutations = 99)
plot(bet.all) #nope
adonis(mse.core.otu ~ zone, data=mse.core.sam, permutations=999)
#yes
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# zone       1    0.4215 0.42153  4.0226 0.12562  0.008 **
#   Residuals 28    2.9341 0.10479         0.87438          
# Total     29    3.3556                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#TAHITI NW
tnw.core.otu <- data.frame(ps.core.tnw@otu_table)
tnw.core.sam <- data.frame(ps.core.tnw@sam_data)

dist.all <- vegdist(tnw.core.otu)
row.names(tnw.core.otu) == row.names(tnw.core.sam)
bet.all <- betadisper(dist.all,tnw.core.sam$zone,bias.adjust = TRUE,type="median")
anova(bet.all) #nope
permutest(bet.all, pairwise = FALSE, permutations = 99)
plot(bet.all) #nope
adonis(tnw.core.otu ~ zone, data=tnw.core.sam, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# zone       1   0.36622 0.36622  4.7046 0.15838  0.008 **
#   Residuals 25   1.94610 0.07784         0.84162          
# Total     26   2.31233                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#biplot - not core
ps.mnw <- subset_samples(ps.rare.trim, site=="MNW")
plot_ordination(ps.mnw, ordinate(ps.rare.trim, "PCoA"), type="biplot", color="Class", shape="zone", title="biplot")

ps.mnw.sam <- ps.mnw@sam_data

ps.mse <- subset_samples(ps.rare.trim, site=="MSE")
plot_ordination(ps.mse, ordinate(ps.mse, "PCoA"), type="split", color="Class", shape="zone", title="biplot")+
  stat_ellipse()

#### BIPLOT - CORE ####
plot_ordination(pseq.core, ordinate(pseq.core, "PCoA"), type="biplot", color="Genus", shape="zone", title="biplot")

#MOOREA NW
ord.mnw <- ordinate(ps.core.mnw, "PCoA", "bray")
plot_ordination(ps.core.mnw, ord.mnw, type="split", color="Genus", shape="zone", title="biplot")#+
theme(legend.position="none")

plot_ordination(ps.core.mnw, ord.mnw,type="biplot",color="Genus", shape="zone")+
  stat_ellipse(group="zone")


stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("BR","FR"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("BR","FR"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  ggtitle("Mo'orea NW")+
  #annotate(geom="text", x=0.7, y=0.2, label="p < 0.01**",size=4)+ #rarefied
  #annotate(geom="text", x=-0.25, y=0.6, label="p < 0.01**",size=4)+ #not rarefied
  xlab("Axis 1 (69.4%)")+ #rarefied
  ylab("Axis 2 (25.6%)")+ #rarefied
  #xlab("Axis 1 (34.3%)")+#non-rarefied
  #ylab("Axis 2 (26.5%)")+#non-rarefied
  theme(axis.text=element_text(size=10))
gg.mnw

#MOOREA SE
ord.mse <- ordinate(ps.core.mse, "PCoA", "bray")
plot_ordination(ps.core.mse, ord.mse, type="biplot", color="Genus", shape="zone", title="biplot")#+
theme(legend.position="none")

quartz()
plot_ordination(ps.core.mse, ord.mse,type="biplot",color="Genus", shape="zone",label="Genus")+
  theme_cowplot()

gg.mse <- plot_ordination(ps.mse, ord.mse,color="zone", shape="zone")+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("BR","FR"))+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("BR","FR"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  ggtitle("Mo'orea NW")+
  #annotate(geom="text", x=0.7, y=0.2, label="p < 0.01**",size=4)+ #rarefied
  #annotate(geom="text", x=-0.25, y=0.6, label="p < 0.01**",size=4)+ #not rarefied
  xlab("Axis 1 (69.4%)")+ #rarefied
  ylab("Axis 2 (25.6%)")+ #rarefied
  #xlab("Axis 1 (34.3%)")+#non-rarefied
  #ylab("Axis 2 (26.5%)")+#non-rarefied
  theme(axis.text=element_text(size=10))
gg.mse

ord.tnw <- ordinate(ps.core.tnw, "PCoA", "bray")
plot_ordination(ps.core.tnw, ord.tnw, type="biplot", color="Genus", shape="zone", title="biplot")#+
theme(legend.position="none")

quartz()
plot_ordination(ps.core.tnw, ord.tnw,type="biplot",color="Genus", shape="zone",label="Genus")+
  theme_cowplot()

#### core + size ####
library('vegan')
size <- read.csv("~/moorea_holobiont/mr_ITS2/mr_size.csv",header=TRUE,row.names=1)

row.names.remove <- c(109)
size2 <- size[!(row.names(size) %in% row.names.remove),]

row.names(samdf) <- samdf$sample
row.names(size2) <- size2$coral_id

samdf.size2 <- merge(samdf,size2,by=0)

row.names(samdf.size2) <- samdf.size2$id

size.rows <- row.names(samdf.size2)

core.otu <- data.frame(pseq.core@otu_table)
otu.size <- core.otu[(row.names(core.otu) %in% size.rows),]

#tahiti 
samdf.size.tnw <- subset(samdf.size2,place=="T")
names <- c(samdf.size.tnw$id)
counts.size.tnw <- otu.size[(rownames(otu.size) %in% names),]
samdf.size.tnw2 <- samdf.size.tnw[(rownames(samdf.size.tnw) %in% row.names(counts.size.tnw)),]

samdf.size.tnw.sorted <- samdf.size.tnw2[sort(rownames(samdf.size.tnw2)),]
counts.size.tnw.sorted <- counts.size.tnw[sort(rownames(counts.size.tnw)),]

row.names(samdf.size.tnw.sorted) == row.names(counts.size.tnw.sorted)

dist.tnw <- vegdist(counts.size.tnw.sorted)
bet.tnw <- betadisper(dist.tnw,samdf.size.tnw.sorted$zone.y)
anova(bet.tnw)
permutest(bet.tnw, pairwise = FALSE, permutations = 99)
plot(bet.tnw) #not sig

adonis(counts.size.tnw.sorted ~ zone.y*size, data=samdf.size.tnw.sorted, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# zone.y       1   0.36622 0.36622  5.0531 0.15838  0.003 **
#   size         1   0.25267 0.25267  3.4864 0.10927  0.071 . 
# zone.y:size  1   0.02652 0.02652  0.3659 0.01147  0.764   
# Residuals   23   1.66691 0.07247         0.72088          
# Total       26   2.31233                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#moorea nw 
samdf.size.mnw <- subset(samdf.size2,place=="MNW")
names <- c(samdf.size.mnw$id)
counts.size.mnw <- otu.size[(rownames(otu.size) %in% names),]
samdf.size.mnw2 <- samdf.size.mnw[(rownames(samdf.size.mnw) %in% row.names(counts.size.mnw)),]

samdf.size.mnw.sorted <- samdf.size.mnw2[sort(rownames(samdf.size.mnw2)),]
counts.size.mnw.sorted <- counts.size.mnw[sort(rownames(counts.size.mnw)),]

row.names(samdf.size.mnw.sorted) == row.names(counts.size.mnw.sorted)

dist.mnw <- vegdist(counts.size.mnw.sorted)
bet.mnw <- betadisper(dist.mnw,samdf.size.mnw.sorted$zone.y)
anova(bet.mnw)
permutest(bet.mnw, pairwise = FALSE, permutations = 99)
plot(bet.mnw) #not sig

adonis(counts.size.mnw.sorted ~ zone.y*size, data=samdf.size.mnw.sorted, permutations=999)
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# zone.y       1   0.07666 0.076664 0.72375 0.02920  0.583
# size         1   0.08426 0.084263 0.79549 0.03210  0.516
# zone.y:size  1   0.02786 0.027859 0.26300 0.01061  0.868
# Residuals   23   2.43631 0.105926         0.92808       
# Total       26   2.62509                  1.00000       

#### accessory ####
ps.rare.trim.otu <- data.frame(ps.rare.trim@otu_table)
core.tax <- data.frame(pseq.core@tax_table)
core.ids <- c(rownames(core.tax))
ps.rare.trim.acc.otu <- ps.rare.trim.otu[,!colnames(ps.rare.trim.otu) %in% core.ids ]

#remake phyloseq object
ps.acc <- phyloseq(otu_table(ps.rare.trim.acc.otu, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa2))
ps.acc #197 taxa accessory

ps.acc.mnw <- subset_samples(ps.acc,site=="MNW")
mnw.acc.otu <- data.frame(ps.acc.mnw@otu_table)
mnw.acc.sam <- data.frame(ps.acc.mnw@sam_data)

dist.all <- vegdist(mnw.acc.otu)
row.names(mnw.acc.otu) == row.names(mnw.acc.sam)
bet.all <- betadisper(dist.all,mnw.acc.sam$zone,bias.adjust = TRUE,type="median")
anova(bet.all) #nope
permutest(bet.all, pairwise = FALSE, permutations = 99)
plot(bet.all) #nope
adonis(mnw.acc.otu ~ zone, data=mnw.acc.sam, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# zone       1    0.6966 0.69663  1.8609 0.06928  0.011 *
#   Residuals 25    9.3589 0.37436         0.93072         
# Total     26   10.0555                 1.00000         

ps.acc.mse <- subset_samples(ps.acc,site=="MSE")
mse.acc.otu <- data.frame(ps.acc.mse@otu_table)
mse.acc.sam <- data.frame(ps.acc.mse@sam_data)

dist.all <- vegdist(mse.acc.otu)
row.names(mse.acc.otu) == row.names(mse.acc.sam)
bet.all <- betadisper(dist.all,mse.acc.sam$zone,bias.adjust = TRUE,type="median")
anova(bet.all) #nope
permutest(bet.all, pairwise = FALSE, permutations = 99)
plot(bet.all) #nope
adonis(mse.acc.otu ~ zone, data=mse.acc.sam, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# zone       1    0.9654 0.96535  2.7769 0.09023  0.001 ***
#   Residuals 28    9.7338 0.34764         0.90977           
# Total     29   10.6992                 1.00000          

ps.acc.tnw <- subset_samples(ps.acc,site=="TNW")
tnw.acc.otu <- data.frame(ps.acc.tnw@otu_table)
tnw.acc.sam <- data.frame(ps.acc.tnw@sam_data)

dist.all <- vegdist(tnw.acc.otu)
row.names(tnw.acc.otu) == row.names(tnw.acc.sam)
bet.all <- betadisper(dist.all,tnw.acc.sam$zone,bias.adjust = TRUE,type="median")
anova(bet.all) #nope
permutest(bet.all, pairwise = FALSE, permutations = 99)
plot(bet.all) #nope
adonis(tnw.acc.otu ~ zone, data=tnw.acc.sam, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# zone       1    0.9984 0.99840  2.7155 0.09798  0.001 ***
#   Residuals 25    9.1915 0.36766         0.90202           
# Total     26   10.1899                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1       

#### accessory + size ####
#tahiti
samdf.size.tnw <- subset(samdf.size2,place=="T")
names <- c(samdf.size.tnw$id)
counts.size.tnw <- ps.rare.trim.acc.otu[(rownames(ps.rare.trim.acc.otu) %in% names),]
samdf.size.tnw2 <- samdf.size.tnw[(rownames(samdf.size.tnw) %in% row.names(counts.size.tnw)),]

samdf.size.tnw.sorted <- samdf.size.tnw2[sort(rownames(samdf.size.tnw2)),]
counts.size.tnw.sorted <- counts.size.tnw[sort(rownames(counts.size.tnw)),]

row.names(samdf.size.tnw.sorted) == row.names(counts.size.tnw.sorted)

dist.tnw <- vegdist(counts.size.tnw.sorted)
bet.tnw <- betadisper(dist.tnw,samdf.size.tnw.sorted$zone.y)
anova(bet.tnw)
permutest(bet.tnw, pairwise = FALSE, permutations = 99)
plot(bet.tnw) #not sig

adonis(counts.size.tnw.sorted ~ zone.y*size, data=samdf.size.tnw.sorted, permutations=999)
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# zone.y       1    0.9984 0.99840 2.70096 0.09798  0.001 ***
#   size         1    0.3185 0.31854 0.86173 0.03126  0.727    
# zone.y:size  1    0.3711 0.37114 1.00403 0.03642  0.451    
# Residuals   23    8.5019 0.36965         0.83434           
# Total       26   10.1899                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#moorea nw 
samdf.size.mnw <- subset(samdf.size2,place=="MNW")
names <- c(samdf.size.mnw$id)
counts.size.mnw <- ps.rare.trim.acc.otu[(rownames(ps.rare.trim.acc.otu) %in% names),]
samdf.size.mnw2 <- samdf.size.mnw[(rownames(samdf.size.mnw) %in% row.names(counts.size.mnw)),]

samdf.size.mnw.sorted <- samdf.size.mnw2[sort(rownames(samdf.size.mnw2)),]
counts.size.mnw.sorted <- counts.size.mnw[sort(rownames(counts.size.mnw)),]

row.names(samdf.size.mnw.sorted) == row.names(counts.size.mnw.sorted)

dist.mnw <- vegdist(counts.size.mnw.sorted)
bet.mnw <- betadisper(dist.mnw,samdf.size.mnw.sorted$zone.y)
anova(bet.mnw)
permutest(bet.mnw, pairwise = FALSE, permutations = 99)
plot(bet.mnw) #not sig

adonis(counts.size.mnw.sorted ~ zone.y*size, data=samdf.size.mnw.sorted, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# zone.y       1    0.6966 0.69663 1.85696 0.06928  0.008 **
#   size         1    0.4516 0.45160 1.20380 0.04491  0.213   
# zone.y:size  1    0.2790 0.27901 0.74375 0.02775  0.855   
# Residuals   23    8.6283 0.37514         0.85806          
# Total       26   10.0555                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### Bar-plots ####
ps_glom <- tax_glom(ps.rare.trim, "Family")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Family")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(legend.position="none")

totalsums <- colSums(seqtab.rare)
summary(totalsums)

#### indicspecies by all sites combined ####
library(indicspecies)

#normal
indval = multipatt(seq.trim, samdf.trim$zone, control = how(nperm=999))
summary(indval)
#10 go to inshore, 12 go to offshore (22 total)

#rarefied - done 3 times
rownames(seq.trim)
indval_rare_1 = multipatt(seq.trim, samdf.rare$zone, control = how(nperm=999))
summary(indval_rare_1) #8 to in, 12 to off

indval_rare_2 = multipatt(seq.trim, samdf.rare$zone, control = how(nperm=999))
summary(indval_rare_2) #7 to in, 11 to off

indval_rare_3 = multipatt(seq.trim, samdf.rare$zone, control = how(nperm=999))
summary(indval_rare_3) #9 to in, 10 to off

#plotting abundances
sqs1 <- data.frame(indval_rare_1[["sign"]])
sqs_sig1 <- subset(sqs1,p.value < 0.05)
sqs_sig1$sqs <- rownames(sqs_sig1)

sqs2 <- data.frame(indval_rare_2[["sign"]])
sqs_sig2 <- subset(sqs2,p.value < 0.05)
sqs_sig2$sqs <- rownames(sqs_sig2)

sqs3 <- data.frame(indval_rare_3[["sign"]])
sqs_sig3 <- subset(sqs3,p.value < 0.05)
sqs_sig3$sqs <- rownames(sqs_sig3)

sqs_12 <- merge(sqs_sig1,sqs_sig2,by="sqs")
sqs <- merge(sqs_sig3,sqs_12,by="sqs")
#saving
write.csv(sqs,"sqs.csv")

goodtaxa <- sqs$sqs
allTaxa <- taxa_names(ps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
indic.rare <- prune_taxa(allTaxa, ps.rare.trim)
plot_bar(indic.rare,x="site_zone",fill="Family")+
  facet_wrap(~Family,scales="free")

#now by more specific results
sqs_in <- subset(sqs,index==1)
sqs_out <- subset(sqs,index==2)

taxa_in <- sqs_in$sqs
allTaxa <- taxa_names(ps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% taxa_in)]
indic.in <- prune_taxa(allTaxa, ps.rare.trim)
indic.in.rel <- transform_sample_counts(indic.in, function(x) x / sum(x))
plot_bar(indic.in.rel,x="site_zone",fill="Family")+
  facet_wrap(~Family,scales="free")

taxa_out <- sqs_out$sqs
allTaxa <- taxa_names(ps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% taxa_out)]
indic.out <- prune_taxa(allTaxa, ps.rare.trim)
plot_bar(indic.out,x="site_zone",fill="Family")+
  facet_wrap(~Family,scales="free")

#### indicspecies by site ####
library(indicspecies)

ps.mnw = subset_samples(ps.rare.trim, site=="MNW")
samdf.mnw <- subset(samdf.rare,site=="MNW")
mnw <- c(row.names(samdf.mnw))
seq.trim.mnw <- seq.trim[row.names(seq.trim) %in% mnw, ]

#rarefied
indval.mnw = multipatt(seq.trim.mnw, samdf.mnw$zone, control = how(nperm=999))
summary(indval.mnw) #6 for in, 5 for off

#benjamini hochberg correction
sqs.mnw <- data.frame(indval.mnw[["sign"]])
sqs.mnw.nona <- sqs.mnw[complete.cases(sqs.mnw),]
sqs.mnw.nona$padj <- p.adjust(sqs.mnw.nona$p.value,method="BH")

#plotting abundances
sqs.mnw2 <- subset(sqs.mnw.nona,padj <= 0.1)
sqs.mnw2$sqs <- rownames(sqs.mnw2)

#saving
write.csv(sqs.mnw2,"sqs_mnw_padj.csv")
sqs.mnw2 <- read.csv("sqs_mnw_padj.csv",row.names=1)

goodtaxa <- sqs.mnw2$sqs
allTaxa <- taxa_names(ps.mnw)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ps.mnw.rel <- transform_sample_counts(ps.mnw, function(x) x / sum(x))
indic.mnw <- prune_taxa(allTaxa, ps.mnw.rel)
plot_bar(indic.mnw,x="zone",fill="Genus")+
  facet_wrap(~Genus,scales="free")

#now by more specific results
sqs_mnwi <- subset(sqs.mnw2,index==1)

#MSE
ps.mse = subset_samples(ps.rare.trim, site=="MSE")
samdf.mse <- subset(samdf.rare,site=="MSE")
mse <- c(row.names(samdf.mse))
seq.trim.mse <- seq.trim[row.names(seq.trim) %in% mse, ]
#rarefied
indval.mse = multipatt(seq.trim.mse, samdf.mse$zone, control = how(nperm=999))
summary(indval.mse) #3 for in, 17 for out

#multitest correction
sqs.mse <- data.frame(indval.mse[["sign"]])
sqs.mse.nona <- sqs.mse[complete.cases(sqs.mse),]
sqs.mse.nona$padj <- p.adjust(sqs.mse.nona$p.value,method="BH")

#plotting abundances
sqs.mse2 <- subset(sqs.mse.nona,padj <= 0.1)
sqs.mse2$sqs <- rownames(sqs.mse2)
#saving
write.csv(sqs.mse2,"sqs_mse_padj.csv")
sqs.mse2 <- read.csv("sqs_mse_padj.csv",row.names=1)

#plotting
goodtaxa <- sqs.mse2$sqs
allTaxa <- taxa_names(ps.mse)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ps.mse.rel <- transform_sample_counts(ps.mse, function(x) x / sum(x))
indic.mse <- prune_taxa(allTaxa, ps.mse.rel)
plot_bar(indic.mse,x="zone",fill="Genus")+
  facet_wrap(~Genus,scales="free")

#now by more specific results
sqs_msei <- subset(sqs.mse2,index==1)
sqs_mseo <- subset(sqs.mse2,index==2)

#TAHITI NW
ps.tnw = subset_samples(ps.rare.trim, site=="TNW")
samdf.tnw <- subset(samdf.rare,site=="TNW")
tnw <- c(row.names(samdf.tnw))
seq.trim.tnw <- seq.trim[row.names(seq.trim) %in% tnw, ]
#rarefied
indval.tnw = multipatt(seq.trim.tnw, samdf.tnw$zone, control = how(nperm=999))
summary(indval.tnw) #3 for in, 17 for out

#multitest correction
sqs.tnw <- data.frame(indval.tnw[["sign"]])
sqs.tnw.nona <- sqs.tnw[complete.cases(sqs.tnw),]
sqs.tnw.nona$padj <- p.adjust(sqs.tnw.nona$p.value,method="BH")

#plotting abundances
sqs.tnw2 <- subset(sqs.tnw.nona,padj <= 0.1)
sqs.tnw2$sqs <- rownames(sqs.tnw2)
#saving
#write.csv(sqs.tnw2,"sqs_tnw_padj.csv")
sqs.tnw2 <- read.csv("sqs_tnw_padj.csv",row.names=1)

#plotting
goodtaxa <- sqs.tnw2$sqs
allTaxa <- taxa_names(ps.tnw)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ps.tnw.rel <- transform_sample_counts(ps.tnw, function(x) x / sum(x))
indic.tnw <- prune_taxa(allTaxa, ps.tnw.rel)
plot_bar(indic.tnw,x="zone",fill="Genus")+
  facet_wrap(~Genus,scales="free")

#now by more specific results
sqs_tnwi <- subset(sqs.tnw2,index==1)
sqs_tnwo <- subset(sqs.tnw2,index==2)

#### bubble plot ####
library(dplyr)
library(tidyverse)
library(grid)
#install.packages("ggplotify")
library("ggplotify")
library(ggpubr)

#MOOREA NW
melt.mnw <- psmelt(indic.mnw)
melt.mnw.sum <- summarySE(melt.mnw,measurevar="Abundance",groupvars = c("OTU","Genus","zone","Species","Family","Order"))

melt.mnw.sum$bool_col <- melt.mnw.sum$OTU %in% sqs_mnwi$sqs
melt.mnw.sum$indic_res <- ifelse(melt.mnw.sum$bool_col == T, "Back reef", "Fore reef")

#have to rename the 2 Pseudomonas instances
#melt.mnw.sum[5:6,2] <- NA
#melt.mnw.sum$Genus = factor(melt.mnw.sum$Genus, levels=c(levels(melt.mnw.sum$Genus), "Pseudomonas.1"))
#melt.mnw.sum$Genus[is.na(melt.mnw.sum$Genus)] = "Pseudomonas.1"

#melt.mnw.sum[11:12,2] <- NA
#melt.mnw.sum$Genus = factor(melt.mnw.sum$Genus, levels=c(levels(melt.mnw.sum$Genus), "Pseudomonas.2"))
#melt.mnw.sum$Genus[is.na(melt.mnw.sum$Genus)] = "Pseudomonas.2"

#MOOREA SE
melt.mse <- psmelt(indic.mse)
melt.mse.sum <- summarySE(melt.mse,measurevar="Abundance",groupvars = c("OTU","Genus","zone","Species","Family","Order"))

melt.mse.sum$bool_col <- melt.mse.sum$OTU %in% sqs_msei$sqs
melt.mse.sum$indic_res <- ifelse(melt.mse.sum$bool_col == T, "Back reef", "Fore reef")

# #2 endozoicos 
# melt.mse.sum[13:14,2] <- NA
# melt.mse.sum$Genus = factor(melt.mse.sum$Genus, levels=c(levels(melt.mse.sum$Genus), "Endozoicomonas.1"))
# melt.mse.sum$Genus[is.na(melt.mse.sum$Genus)] = "Endozoicomonas.1"
# 
# melt.mse.sum[33:34,2] <- NA
# melt.mse.sum$Genus = factor(melt.mse.sum$Genus, levels=c(levels(melt.mse.sum$Genus), "Endozoicomonas.2"))
# melt.mse.sum$Genus[is.na(melt.mse.sum$Genus)] = "Endozoicomonas.2"

#TAHITI NW
melt.tnw <- psmelt(indic.tnw)
melt.tnw.sum <- summarySE(melt.tnw,measurevar="Abundance",groupvars = c("OTU","Genus","zone","Species","Family","Order"))

melt.tnw.sum$bool_col <- melt.tnw.sum$OTU %in% sqs_tnwi$sqs
melt.tnw.sum$indic_res <- ifelse(melt.tnw.sum$bool_col == T, "Back reef", "Fore reef")

#### indicspecies sites all at once ####
melt.mnw.sum$site <- rep("mnw",2)
melt.mse.sum$site <- rep("mse",14)
melt.tnw.sum$site <- rep("tnw",10)
allsites <- rbind(melt.mnw.sum,melt.mse.sum,melt.tnw.sum)

# #genera as indicators across both fore & back reef
allsites[c(3:4,19:20),2] <- NA
allsites$Genus[is.na(allsites$Genus)] = "Lacibacter*"

allsites[c(5:6,21:22),2] <- NA
allsites$Genus[is.na(allsites$Genus)] = "Phreatobacter*"
# 
# allsites[c(31:32,85:86),2] <- NA
# allsites$Genus = factor(allsites$Genus, levels=c(levels(allsites$Genus), "Methylobacterium*"))
# allsites$Genus[is.na(allsites$Genus)] = "Methylobacterium*"

allsites$Genus <- as.factor(allsites$Genus)
genera <- levels(allsites$Genus)
write.csv(genera,file="genera2_padj.csv")

##all of them:
#allsites$Genus <- factor(allsites$Genus, levels=rev(c("Pseudomonas.1", 	"Pseudomonas.2", 	"Rheinheimera", 	"Rubrobacter", 	"Acinetobacter", 	"Class_Bacteroidia", 	"Order_Entomoplasmatales", 	"Family_Simkaniaceae", 	"Methylobacterium", 	"Enhydrobacter", 	"Finegoldia", 	"Peptoniphilus", 	"Enhydrobacter*", 	"Prochlorococcus_MIT9313", 	"Pseudomonas", 	"Algoriphagus", 	"Curvibacter", 	"DSSD61", 	"Endozoicomonas", 	"Family_Neisseriaceae", 		"Haemophilus", 	"Lacibacter", 	"Lactococcus", 	"Order_Obscuribacterales", 	"Phreatobacter", 	"Rhodobacter", 	"Veillonella", 	"Endozoicomonas.1", 	"Endozoicomonas.2", 	"Acidovorax", 	"Alteromonas", 	"Rubrobacter*", "Methylobacterium*","Asinibacterium", 	"Cloacibacterium", 	"Kingdom_Bacteria", 	"Order_Rhodospirillales", 	"Pseudoalteromonas")))
##after multiple test correction:
allsites$Genus <- factor(allsites$Genus, levels=rev(c("Pseudomonas","Family_Simkaniaceae","Phreatobacter*","Lacibacter*","Rubrobacter","Curvibacter","Endozoicomonas","Enhydrobacter","Order_Obscuribacterales","Kingdom_Bacteria","Acidovorax")))

allsites$site <- gsub("mnw","Moorea NW",allsites$site)
allsites$site <- gsub("mse","Moorea SE",allsites$site)
allsites$site <- gsub("tnw","Tahiti NW",allsites$site)

allsites$indic_res <- gsub("Back reef","Back reef indicator",allsites$indic_res)
allsites$indic_res <- gsub("Fore reef","Fore reef indicator",allsites$indic_res)

quartz()
ggplot(allsites, aes(x = zone, y = Genus)) + 
  geom_point(aes(size = Abundance,fill=Genus), alpha = 0.75, shape = 21)+
  facet_grid(indic_res~site,scales="free",space="free")+
  theme_cowplot()+
  ylab("Genus")+
  xlab('Reef zone')+
  guides(fill=FALSE)+
  scale_x_discrete(labels=c("BR","FR"))

#### indicspecies by genus ####
library(indicspecies)

#creating counts by genus
ps.glom.genus.temp <- tax_glom(ps.rare.trim, "Genus")
ps.glom.genus.temp #down to 124 taxa
seq.glom.genus <- data.frame(ps.glom.genus.temp@otu_table)
tax.glom.genus <- as.data.frame(ps.glom.genus.temp@tax_table@.Data)
sqs <- rownames(tax.glom.genus)

rownames(tax.glom.genus) == sqs #check if names match up!
colnames(seq.glom.genus) <- tax.glom.genus$Genus
rownames(tax.glom.genus) <- tax.glom.genus$Genus
rownames(tax.glom.genus) == colnames(seq.glom.genus)
tax.glom.genus <- as.matrix(tax.glom.genus)

#renamed ps object
ps.glom.genus <- phyloseq(otu_table(seq.glom.genus, taxa_are_rows=FALSE), 
                          sample_data(samdf.rare), 
                          tax_table(tax.glom.genus))

#indicspecies
indval.genus <- multipatt(seq.glom.genus, samdf.rare$zone, control = how(nperm=999))
summary(indval.genus)
#cool stuff

#now by site
#moorea NW
samdf.mnw <- subset(samdf.rare,site=="MNW")
mnw <- c(row.names(samdf.mnw))
seq.glom.mnw <- seq.glom.genus[row.names(seq.glom.genus) %in% mnw, ]

row.names(samdf.mnw) == row.names(seq.glom.mnw)
indval.glom.mnw <- multipatt(seq.glom.mnw, samdf.mnw$zone, control = how(nperm=999))
summary(indval.glom.mnw) #4 in, 6 out

#plotting abundances
ps.glom.mnw = subset_samples(ps.glom.genus, site=="MNW")
genus.mnw <- data.frame(indval.glom.mnw[["sign"]])
genus.mnw <- subset(genus.mnw,p.value <= 0.05)
genus.mnw$genus <- rownames(genus.mnw)
#saving
write.csv(genus.mnw,"genus_mnw.csv")
genus.mnw <- read.csv("genus_mnw.csv",row.names=1)

#plotting abundance
goodtaxa <- genus.mnw$genus
allTaxa <- row.names(tax.glom.genus)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
indic.mnw <- prune_taxa(allTaxa, ps.glom.mnw)
indic.mnw.rel <- transform_sample_counts(indic.mnw, function(x) x / sum(x))
plot_bar(indic.mnw.rel,x="zone",fill="Genus")+
  facet_wrap(~Genus,scales="free")

#now by more specific results
genus.mnwi <- subset(genus.mnw,index==1)
genus.mnwo <- subset(genus.mnw,index==2)

#moorea SE
samdf.mse <- subset(samdf.rare,site=="MSE")
mse <- c(row.names(samdf.mse))
seq.glom.mse <- seq.glom.genus[row.names(seq.glom.genus) %in% mse, ]

row.names(samdf.mse) == row.names(seq.glom.mse)
indval.glom.mse <- multipatt(seq.glom.mse, samdf.mse$zone, control = how(nperm=999))
summary(indval.glom.mse) #2 in, 13 out

#plotting abundances
ps.glom.mse = subset_samples(ps.glom.genus, site=="MSE")
genus.mse <- data.frame(indval.glom.mse[["sign"]])
genus.mse <- subset(genus.mse,p.value <= 0.05)
genus.mse$genus <- rownames(genus.mse)
#saving
write.csv(genus.mse,"genus_mse.csv")
genus.mse <- read.csv("genus_mse.csv",row.names=1)

#plotting abundance
goodtaxa <- genus.mse$genus
allTaxa <- row.names(tax.glom.genus)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
indic.mse <- prune_taxa(allTaxa, ps.glom.mse)
indic.mse.rel <- transform_sample_counts(indic.mse, function(x) x / sum(x))
plot_bar(indic.mse.rel,x="zone",fill="Genus")+
  facet_wrap(~Genus,scales="free")

#now by more specific results
genus.msei <- subset(genus.mse,index==1)
genus.mseo <- subset(genus.mse,index==2)

#tahiti NW
samdf.tnw <- subset(samdf.rare,site=="TNW")
tnw <- c(row.names(samdf.tnw))
seq.glom.tnw <- seq.glom.genus[row.names(seq.glom.genus) %in% tnw, ]

row.names(samdf.tnw) == row.names(seq.glom.tnw)
indval.glom.tnw <- multipatt(seq.glom.tnw, samdf.tnw$zone, control = how(nperm=999))
summary(indval.glom.tnw) #8 in, 2 out

#plotting abundances
ps.glom.tnw = subset_samples(ps.glom.genus, site=="TNW")
genus.tnw <- data.frame(indval.glom.tnw[["sign"]])
genus.tnw <- subset(genus.tnw,p.value <= 0.05)
genus.tnw$genus <- rownames(genus.tnw)
#saving
write.csv(genus.tnw,"genus_tnw.csv")
genus.tnw <- read.csv("genus_tnw.csv",row.names=1)

#plotting abundance
goodtaxa <- genus.tnw$genus
allTaxa <- row.names(tax.glom.genus)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
indic.tnw <- prune_taxa(allTaxa, ps.glom.tnw)
indic.tnw.rel <- transform_sample_counts(indic.tnw, function(x) x / sum(x))
plot_bar(indic.tnw.rel,x="zone",fill="Genus")+
  facet_wrap(~Genus,scales="free")

#now by more specific results
genus.tnwi <- subset(genus.tnw,index==1)
genus.tnwo <- subset(genus.tnw,index==2)

#### bubble plot - genus ####
library(dplyr)

#MOOREA NW
melt.mnw <- psmelt(indic.mnw.rel)
melt.mnw$Abundance[is.na(melt.mnw$Abundance)] <- 0
melt.mnw.sum <- summarySE(melt.mnw,measurevar="Abundance",groupvars = c("Genus","zone"))

melt.mnw.sum$bool_col <- melt.mnw.sum$Genus %in% genus.mnwi$genus
melt.mnw.sum$indic_res <- ifelse(melt.mnw.sum$bool_col == T, "Back reef", "Fore reef")

bub.genus.mnw <- ggplot(melt.mnw.sum, aes(x = zone, y = Genus)) + 
  geom_point(aes(size = Abundance, fill = Genus), alpha = 0.75, shape = 21)+
  facet_grid(indic_res~.,scales="free",space="free")+
  theme_cowplot()+
  #scale_fill_manual(values=cols.mnw, guide=FALSE)+
  ggtitle("Moorea NW")+
  ylab("ASV")+
  xlab("")+
  scale_x_discrete(labels=c("Back","Fore"))
bub.genus.mnw

#MOOREA SE
melt.mse <- psmelt(indic.mse.rel)
melt.mse$Abundance[is.na(melt.mse$Abundance)] <- 0
melt.mse.sum <- summarySE(melt.mse,measurevar="Abundance",groupvars = c("Genus","zone"))

melt.mse.sum$bool_col <- melt.mse.sum$Genus %in% genus.msei$genus
melt.mse.sum$indic_res <- ifelse(melt.mse.sum$bool_col == T, "Back reef", "Fore reef")

bub.genus.mse <- ggplot(melt.mse.sum, aes(x = zone, y = Genus)) + 
  geom_point(aes(size = Abundance, fill = Genus), alpha = 0.75, shape = 21)+
  facet_grid(indic_res~.,scales="free",space="free")+
  theme_cowplot()+
  #scale_fill_manual(values=cols.mse, guide=FALSE)+
  ggtitle("Moorea SE")+
  ylab("ASV")+
  xlab("")+
  scale_x_discrete(labels=c("Back","Fore"))
bub.genus.mse

#TAHITI NW
melt.tnw <- psmelt(ps.glom.tnw)
melt.tnw$Abundance[is.na(melt.tnw$Abundance)] <- 0
melt.tnw.sum <- summarySE(melt.tnw,measurevar="Abundance",groupvars = c("Genus","zone"))

melt.tnw.sum$bool_col <- melt.tnw.sum$Genus %in% genus.tnwi$genus
melt.tnw.sum$indic_res <- ifelse(melt.tnw.sum$bool_col == T, "Back reef", "Fore reef")

bub.genus.tnw <- ggplot(melt.tnw.sum, aes(x = zone, y = Genus)) + 
  geom_point(aes(size = Abundance, fill = Genus), alpha = 0.75, shape = 21)+
  facet_grid(indic_res~.,scales="free",space="free")+
  theme_cowplot()+
  #scale_fill_manual(values=cols.tnw, guide=FALSE)+
  ggtitle("Tahiti NW")+
  ylab("ASV")+
  xlab("")+
  scale_x_discrete(labels=c("Back","Fore"))
bub.genus.tnw

#### presence absence #### 
#first genus
seq.glom.genus.pa <- seq.glom.genus
seq.glom.genus.pa[seq.glom.genus.pa>0] <-1

#### DESEQ to find differentially abundant ASVs ####
library(DESeq2)

ds.all = phyloseq_to_deseq2(ps.rare.trim, ~ zone)
dds.all <- estimateSizeFactors(ds.all,type="poscounts")
stat.all = DESeq(dds.all, test="Wald", fitType="parametric")
res = results(stat.all, cooksCutoff = FALSE)
alpha = 0.05
sigtab.all = res[which(res$padj < alpha), ]
sigtab.all = cbind(as(sigtab.all, "data.frame"), as(tax_table(ps.rare.trim)[rownames(sigtab.all), ], "matrix"))
sigtab.all
dim(sigtab.all) #10!  
write.csv(sigtab.all,"sigtab.all_justzone.csv")

ps.mnw = subset_samples(ps.rare.trim, site=="MNW")
#ps.mnw = subset_samples(ps.trim, site=="MNW")
ps.mnw.0 <- prune_taxa(taxa_sums(ps.mnw)>0,ps.mnw) #remove samples with 0 total

ds.mnw = phyloseq_to_deseq2(ps.mnw.0, ~ zone)
dds.mnw <- estimateSizeFactors(ds.mnw,type="poscounts")
stat.mnw = DESeq(dds.mnw, test="Wald", fitType="parametric")

res = results(stat.mnw, cooksCutoff = FALSE)
alpha = 0.05
sigtab.mnw = res[which(res$padj < alpha), ]
sigtab.mnw = cbind(as(sigtab.mnw, "data.frame"), as(tax_table(ps.mnw)[rownames(sigtab.mnw), ], "matrix"))
sigtab.mnw
dim(sigtab.mnw)
write.csv(sigtab.mnw,"sigtab.mnw.csv")
#7 for MNW after rarefying, removing low count samples & taxa, rd2 based

goodtaxa <- c(row.names(sigtab.mnw))
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex1 = prune_taxa(allTaxa, ps.rare)
ps.mnw = subset_samples(ex1, site=="MNW")

quartz()
plot_bar(ps.mnw, x="zone", fill="class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  ggtitle("Moorea NW")+
  facet_wrap(~class)

ps.mse = subset_samples(ps.rare.trim, site=="MSE")
#ps.mse = subset_samples(ps.trim, site=="MSE")
ps.mse.0 <- prune_taxa(taxa_sums(ps.mse)>0,ps.mse) #remove samples with 0 total

ds.mse = phyloseq_to_deseq2(ps.mse.0, ~ zone)
dds.mse <- estimateSizeFactors(ds.mse,type="poscounts")
stat.mse = DESeq(dds.mse, test="Wald", fitType="parametric")

res = results(stat.mse, cooksCutoff = FALSE)
alpha = 0.05
sigtab.mse = res[which(res$padj < alpha), ]
sigtab.mse = cbind(as(sigtab.mse, "data.frame"), as(tax_table(ps.mse)[rownames(sigtab.mse), ], "matrix"))
dim(sigtab.mse)
sigtab.mse
write.csv(sigtab.mse,"sigtab.mse.csv")
#5 for MSE 
#6 without rarefying

goodtaxa <- c(row.names(sigtab.mse))
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex1 = prune_taxa(allTaxa, ps.rare)
ps.mse = subset_samples(ex1, site=="MSE")

quartz()
plot_bar(ps.mse, x="zone", fill="class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  ggtitle("Moorea SE")+
  facet_wrap(~class)

ps.tnw = subset_samples(ps.rare.trim, site=="TNW")
#ps.t = subset_samples(ps.trim, site=="TNW")
ps.tnw.0 <- prune_taxa(taxa_sums(ps.tnw)>0,ps.tnw) #remove taxa with 0 total

ds.tnw = phyloseq_to_deseq2(ps.tnw.0, ~ zone)
dds.tnw <- estimateSizeFactors(ds.tnw,type="poscounts")
stat.tnw = DESeq(dds.tnw, test="Wald", fitType="parametric")

res = results(stat.tnw, cooksCutoff = FALSE)
alpha = 0.05
sigtab.tnw = res[which(res$padj < alpha), ]
sigtab.tnw = cbind(as(sigtab.tnw, "data.frame"), as(tax_table(ps.tnw)[rownames(sigtab.tnw), ], "matrix"))
dim(sigtab.tnw)
sigtab.tnw
#9 for TNW
write.csv(sigtab.tnw,"sigtab.tah.csv")

#how many are shared?
merge(sigtab.mnw,sigtab.mse,by=0) #none
merge(sigtab.mse,sigtab.tnw,by=0) #just 41
merge(sigtab.mnw,sigtab.tnw,by=0) #none

theme_set(theme_cowplot())
#scale_fill_discrete <- function(palname = "Set1", ...) {
#  scale_fill_brewer(palette = palname, ...)
#}
# Phylum orderc
#x = tapply(sigtab.mnw$log2FoldChange, sigtab.mnw$phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtab.mnw$Phylum = factor(as.character(sigtab.mnw$phylum), levels=names(x))
# class 

#not sure if I need any of this?
# x = tapply(sigtab.mnw$log2FoldChange, sigtab.mnw$Genus, function(x) max(x))
# x = sort(x, TRUE)
# sigtab.mnw$Genus = factor(as.character(sigtab.mnw$Genus), levels=names(x))
# quartz()

sigtab.mnw$site <- c(rep("Moorea NW",7))
sigtab.mse$site <- c(rep("Moorea SE",4))
sigtab.tnw$site <- c(rep("Tahiti NW",9))

sigtab.toplot <- rbind(sigtab.mnw,sigtab.tnw,sigtab.mse)

deseq.genera <- c(levels(sigtab.toplot$Genus))
sigtab.toplot$Genus <- factor(sigtab.toplot$Genus,levels=c("Finegoldia","Peptoniphilus","Pseudomonas","Acinetobacter","Class_Bacteroidia","Anoxybacillus","Order_Entomoplasmatales","Endozoicomonas","Phreatobacter","Acidovorax","Kingdom_Bacteria","Lacibacter","Methylobacterium","Ralstonia","Family_Burkholderiaceae"))

quartz()
ggplot(sigtab.toplot,aes(x=Genus,y=log2FoldChange,color=Genus))+
  geom_point(size=6)+
  facet_grid(site~.)+
  theme_bw()+
  theme(legend.position = "none",axis.text.x=element_text(hjust=1,angle=45),text=element_text(size=14))+
  ylim(-25,25)+
  ylab('Log2 fold change')+
  geom_hline(yintercept=0,linetype="dashed")

#### outdated maybe ####
goodtaxa = c("sq4","sq16","sq26")
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex1 = prune_taxa(allTaxa, ps.rare)
mnw1 <- subset_samples(ex1,site=="MNW")
quartz()
plot_bar(mnw1, x="zone", fill="class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45))+
  ggtitle("MNW")

goodtaxa = c("sq16","sq26","sq54","sq56","sq65")
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex2 = prune_taxa(allTaxa, ps.rare)
mnw <- subset_samples(ex2,site=="MNW")

goodtaxa = c("sq16","sq26","sq54","sq56","sq65")
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex2 = prune_taxa(allTaxa, ps.rare)
mnw <- subset_samples(ex2,site=="MNW")

quartz()
plot_bar(mnw, x="zone", fill="genus")+
  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  ggtitle("Moorea NW")+
  facet_wrap(~genus)

#### function ####

#### DESEQ to normalize counts ####
#BiocManager::install("DESeq2")
library("DESeq2")
#BiocManager::install("edgeR")
#BiocManager::install("DEFormats")
library(edgeR)
library(DEFormats)
library("MCMC.OTU")
library("cowplot")

#count data
seq <- read.csv("mr16s_AllOTUs_mcmc.csv")
#getting rid of crap ones
#mcmc.otu requires a blank column with sample names first, did manually in excel

seq2 <- purgeOutliers(seq,count.columns=10:2285,sampleZcut=-2.5,otu.cut=0.0001)
#bad samples: 304 (F9)
#50 OTUs pass frequency cutoff

row.names(seq2) <- seq2$id #setting my sample names as row names so I can remove the sample name column
seq2 <- seq2[,10:263] #now removing sample name columns
seqt <- t(seq2) #deseq requires count data in columns

#sample data
sam <- read.csv("mr16s_sampledata.csv")
row.names(sam) <- sam$id
row.names.remove <- c("F9")
sam <- sam[!(row.names(sam) %in% row.names.remove), ]

#### started here ####
seqt <- t(seqtab.trim) #deseq requires count data in columns

#method from here:
#https://astrobiomike.github.io/amplicon/dada2_workflow_ex#dereplication
deseq_counts <- DESeqDataSetFromMatrix(seqt, colData = samdf, design = ~site*zone) 
# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter

deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that. 

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
vst_trans_count_tab[vst_trans_count_tab < 0.0] <- 0.0 #getting rid of negatives

# back into phyloseq
vstt <- t(vst_trans_count_tab)

ps.norm <- phyloseq(otu_table(vstt, taxa_are_rows=FALSE), 
                    sample_data(samdf), 
                    tax_table(taxa2))
ps.norm

ps.mnw = subset_samples(ps.norm, site=="MNW")
ds.mnw = phyloseq_to_deseq2(ps.mnw, ~ zone)
dds.mnw <- estimateSizeFactors(ds.mnw,type="poscounts")
stat.mnw = DESeq(dds.mnw, test="Wald", fitType="parametric")

res = results(stat.mnw, cooksCutoff = FALSE)
alpha = 0.05
sigtab.mnw = res[which(res$padj < alpha), ]
sigtab.mnw = cbind(as(sigtab.mnw, "data.frame"), as(tax_table(ps.mnw)[rownames(sigtab.mnw), ], "matrix"))
head(sigtab.mnw)
dim(sigtab.mnw)

#### not sure what this was ####

# sam$new.norm <- sample_sums(ps.norm)
# plot(new.norm~site_zone,data=sam) #more normalized now
# a1 <- aov(new.norm~site_zone,data=sam)
# summary(a1) #p = 0.0646 lol
# plot(new.norm~island,data=sam) #even
# plot(new.norm~site,data=sam) #MSE still a little more
# plot(new.norm~in_off,data=sam) #outer more
# a1 <- aov(new.norm~in_off,data=sam)
# summary(a1) # p = significant

#### rarefy - outdated ####
library(vegan)
library(MCMC.OTU)
library(phyloseq)
setwd("~/moorea_holobiont/mr_16S/")
samdf<-read.csv("mr16s_sampledata.csv")
rownames(samdf) <- samdf$id

rarecurve(seqtab.nochim,step=100,label=FALSE)
total <- rowSums(seqtab.nochim)
subset(total, total <10000)
row.names.remove <- c("A5","B5","B6","B8","C2","C8","D12","F8","F9","G9","H3","H4","H5","H6")
seqtab.nochim <- seqtab.nochim[!(row.names(seqtab.nochim) %in% row.names.remove), ]
samdf <- samdf[!(row.names(samdf) %in% row.names.remove), ]

seq.rare <- rrarefy(seqtab.nochim,sample=10000)
rarecurve(seq.rare,step=100,label=FALSE)

write.csv(seq.rare,"~/moorea_holobiont/mr_16S/seq.rare10k_2.csv")
#added the sample column that mcmcotu requires
seq.rare <- read.csv("~/moorea_holobiont/mr_16S/seq.rare10k_2.csv")
str(seq.rare)

seq2 <- purgeOutliers(seq.rare,count.columns=3:2278,sampleZcut=-2.5,otu.cut=0.001)
#keeps 50 OTUs

row.names(seq2) <- seq2$sample #setting my sample names as row names so I can remove the sample name column
seq2 <- seq2[,3:52] #now removing sample name columns
row.names(samdf) <- samdf$id

ps.rare <- phyloseq(otu_table(seq2, taxa_are_rows=FALSE), 
                    sample_data(samdf), 
                    tax_table(taxid))

#making ids easier
ids <- paste0("sq", seq(1, length(colnames(seq2))))
colnames(seq2)<-ids
##replace sequences with shorter names (correspondence table output below)
taxa_names(ps.rare)<-ids

#### Pcoa #####
write.csv(vstt,"~/moorea_holobiont/mr_16S/vstt.16s.csv")
write.csv(seq.rare,"~/moorea_holobiont/mr_16s/seq.rare10k.csv")

#before normalizing
ps.ord <- ordinate(ps_no87, method="MDS", distance="euclidean")
ps.ord.eigs <- ps.ord$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
plot_ordination(ps_no87, ps.ord, color="in_off") + 
  geom_point(size=1) #+ 

#after normalizing
vst_pcoa <- ordinate(ps.norm, method="MDS", distance="euclidean")
#vst_pcoa <- ordinate(ps.rare, method="MDS", distance="euclidean")
#^this looks horrible
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
quartz()
plot_ordination(ps.norm, vst_pcoa, color="zone",shape="zone") + 
  geom_point(size=1) +
  facet_wrap(~site)+
  stat_ellipse(level=0.95)+
  theme_cowplot()+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_shape_manual(values=c(16,15),labels=c("Backreef","Forereef"))+
  labs(color="Reef zone",shape="Reef zone")
#geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
#coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA")
#scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
#theme(legend.position="none")

#after rarefying
rare_pcoa <- ordinate(ps.rare, method="MDS", distance="euclidean")
#vst_pcoa <- ordinate(ps.rare, method="MDS", distance="euclidean")
#^this looks horrible
eigen_vals <- rare_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
quartz()
plot_ordination(ps.rare, rare_pcoa, color="zone",shape="zone") + 
  geom_point(size=1) +
  facet_wrap(~site)+
  stat_ellipse(level=0.95)+
  theme_cowplot()+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_shape_manual(values=c(16,15),labels=c("Backreef","Forereef"))+
  labs(color="Reef zone",shape="Reef zone")
#geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
#coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA")
#scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
#theme(legend.position="none")

#trying to select only the variables I want
#resource for analysis here:
#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
#https://joey711.github.io/phyloseq-extensions/DESeq2.html

ds.all = phyloseq_to_deseq2(ps.rare, ~ site*zone)
dds.all <- estimateSizeFactors(ds.all,type="poscounts")
stat.all = DESeq(dds.all, test="Wald", fitType="parametric")

res = results(stat.all, cooksCutoff = FALSE)
alpha = 0.05
sigtab.all = res[which(res$padj < alpha), ]
sigtab.all= cbind(as(sigtab.all, "data.frame"), as(tax_table(ps.rare)[rownames(sigtab.all), ], "matrix"))
head(sigtab.all)
dim(sigtab.all)

goodtaxa <- c(row.names(sigtab.all))
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex1 = prune_taxa(allTaxa, ps.rare)

quartz()
plot_bar(ex1, x="zone", fill="class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45))+
  ggtitle("MNW")+
  facet_wrap(~site+genus)


#### bar plot
samdf$newname <- paste(samdf$site_zone,samdf$sample,sep="_")
row.names(samdf) <- samdf$newname
row.names(seq2) <- samdf$newname

taxid2 <- tax_table(ps.rare)

ps.rare2 <- phyloseq(otu_table(seq2, taxa_are_rows=FALSE), 
                     sample_data(samdf), 
                     tax_table(taxid2))

top3 <- names(sort(taxa_sums(ps.rare2), decreasing=TRUE))[1:118]
ps.top3 <- transform_sample_counts(ps.rare2, function(OTU) OTU/sum(OTU))
ps.top3 <- prune_taxa(top3, ps.top3)
quartz()
plot_bar(ps.top3, x="newname",fill="class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

#summed bar plot
ps <- tax_glom(ps.rare, "family")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))

ps.melt <- psmelt(ps2)

ps.melt$site_zone <- gsub(1,"MNW-B",ps.melt$site_zone)
ps.melt$site_zone <- gsub(2,"MNW-F",ps.melt$site_zone)
ps.melt$site_zone <- gsub(3,"MSE-B",ps.melt$site_zone)
ps.melt$site_zone <- gsub(4,"MSE-F",ps.melt$site_zone)
ps.melt$site_zone <- gsub(5,"TNW-B",ps.melt$site_zone)
ps.melt$site_zone <- gsub(6,"TNW-F",ps.melt$site_zone)

ggplot(ps.melt,aes(x=site_zone,y=Abundance,fill=family))+
  geom_bar(stat="identity", colour="black")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=0),text=element_text(family="Gill Sans MT"))+
  xlab('Site')+
  labs(fill="Family")

quartz()
tb2 <- tb%>%
  group_by(site_zone,family)%>%
  summarize_at("Abundance",sum)
plot_bar(ps.all, x="site_zone",fill="family") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

#############################################

# Principal coordinate analysis 

library(vegan)
#install.packages("MCMC.OTU")
library(MCMC.OTU)
#install.packages("ggfortify")
library(ggfortify)
library(cluster)
#install.packages('labdsv')
library(labdsv)

#have to add a first column with no name numbered 1-93
alldat<-read.csv("mr16s_AllOTUs_mcmc.csv")
names(alldat)
str(alldat)
head(alldat) 
#dat<-alldat[c(1:3,5:110),] #to remove rows with missing data, but I don't have any

# purging under-sequenced samples; and OTUs represented in less than 3% of all samples (~must occur at least 3 separate times)
goods2=purgeOutliers(alldat,count.columns=10:2285,sampleZcut=-2.5,otu.cut=0.032)

# creating a log-transfromed normalized dataset for PCoA:
goods.log=logLin(data=goods2,count.columns=10:length(names(goods2)))

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
goods.dist=vegdist(goods.log,method="manhattan")
goods.pcoa=pcoa(goods.dist)

scores=goods.pcoa$vectors
conditions=goods2[,1:9]
margin=0.01
quartz()
plot(scores[,1], scores[,2],type="n",
     xlim=c(min(scores[,1])-margin,max(scores[,1])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     #mgp=c(2.3,1,0),
     xlab=paste("Axis1 (",round(goods.pcoa$values$Relative_eig[1]*100,1),"%)",sep=""),
     ylab=paste("Axis2 (",round(goods.pcoa$values$Relative_eig[2]*100,1),"%)",sep=""),
     main="")
points(scores[conditions$zone=="in",1],scores[conditions$zone=="in",2],pch=1)
points(scores[conditions$zone=="out",1],scores[conditions$zone=="out",2],pch=19)
legend("bottomright", c("Inshore","Offshore"), pch=c(1, 19), cex=0.8)

#prettier PCA
scores=goods.pcoa$vectors
conditions=goods[,1:9]
zone <- conditions$zone

p1 <- ggordiplots::gg_ordiplot(ord=goods.pcoa$vectors,groups=zone,scaling=0.8,choices=c(1,2),conf=0.95,spiders=TRUE,pt.size=2)
plot <- p1$plot
plot + xlab('Axis1 (19.%)')+ 
  ylab('Axis2 (14.42%)')+
  scale_shape_manual(values=c(1,5))
names(p1)
spid <- p1$df_spiders
p1$df_spiders
gg <- p1$df_ord
quartz()
ggplot(gg,aes(x=x,y=y,group=Group,color=Group,shape=Group,fill=Group))+
  geom_point(size=2)+
  geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.95)+
  xlab('Axis 1 (19.5%)')+
  ylab('Axis 2 (14.4%)')+
  theme_classic()+
  theme(text=element_text(family="Gill Sans MT"))+
  scale_shape_manual(values=c(21,24),labels=c("Inshore","Offshore"))+
  scale_fill_manual(values=c("coral1","cyan3"),labels=c("Inshore","Offshore"))+
  scale_color_manual(values=c("coral1","cyan3"),labels=c("Inshore","Offshore"))+
  labs(shape="Reef zone",color="Reef zone",fill="Reef zone")


