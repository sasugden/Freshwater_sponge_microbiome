##################################################### DATA IMPORT AND QUALITY CONTROL ######################################
#### PREPARE WORKSPACE #####
library(phyloseq)
library(decontam)
library(ALDEx2)
library(vegan)
library(DECIPHER)
library(phangorn)
library(Biostrings)
library(ShortRead)
library(picante)
library(iNEXT)
library(iNextPD)
library(dada2)
library(randomForest)
library(ggtree)
library(ggplot2)
library(gridExtra)

sample_data <- read.csv("~/Documents/sponges_final/sponge_metadata_2020-10-8.csv")
rownames(sample_data) <- sample_data$SampleID
sample_data$Type <- factor(sample_data$Type, levels=c("sponge", "water", "biofilm"))
sample_data$Source <- factor(sample_data$Source, levels=c("Sooke", "Nanaimo", "Cowichan"))

#### DADA2 - RAW SEQUENCE PROCESSING TO GENERATE ASVs ####
# Initial sequence processing (based on the DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html)
# Note that all our sequence data from a single sequencing run was processed together using this code.
# We then separated the sponge data for this manuscript.
# Because of how DADA2 generates ASVs, running this code using only the sequences from the sponge manuscript may therefore produce slightly different results.
# Other sequencing data that was processed concurrently came from:
# Sugden et al. (2020) - An altered microbiome in urban coyotes... (Scientific Reports)
# Sugden et al. (2021) - Individual and site-specific variation... (Microbial Ecology)

# Set a seed so that analyses are reproducible.
set.seed(100)

# Name the folder where the sequence files are stored.
input_path <- "raw_data/16S_sequences/"

# Identify forward (Fs) and reverse (Rs) reads.
fnFs <- sort(list.files(input_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(input_path, pattern="_R2_001.fastq"))

# Identify sample names (should be the beginning of the filename).
sampleNames <- sapply(strsplit(fnFs, "_"), '[', 1)
fnFs <- file.path(input_path, fnFs)
fnRs <- file.path(input_path, fnRs)

# Create a new directory where filtered reads will be stored.
filt_path <- file.path(input_path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)

# Name filtered sequences with the same names as input sequences with _filt added.
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

# This filtering step allows a maximum of two expected errors per read, and zero ambiguous (N) bases.
# It also trims the forward reads to 240 bp and the reverse reads to 160 bp.
# The "out" table gives you a table with the number of input and output reads from the QC step.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# Dereplicate sequences.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Assign sample names.
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

# Calculate error probabilities.
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Use error probabilities to determine exact sequence variants.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

# Merged paired-end reads and construct an OTU table.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtabAll <- makeSequenceTable(mergers)

# Remove chimeras.
seqtabNoC <- removeBimeraDenovo(seqtabAll, multithread=TRUE)

# Examine the distribution of sequence lengths. Most sequences should be ~253bp.
table(nchar(getSequences(seqtabNoC)))

# Save a sequence table with ASVs that don't fall within the expected sequence length.
seqtab.wrong.length.asvs <- seqtabAll[,nchar(colnames(seqtabNoC)) %in% seq(239,249) |
                                        nchar(colnames(seqtabNoC)) %in% seq(257,400)]
table(nchar(getSequences(seqtab.wrong.length.asvs)))

# Trim sequence table to remove amplicons that are not within 250-256bp.
seqtabNoC <- seqtabNoC[,nchar(colnames(seqtabNoC)) %in% seq(250,256)]
table(nchar(getSequences(seqtabNoC)))

# This section produces a table that shows how many reads were lost during the processing steps.
getN <- function(x) sum(getUniques(x))
read_counts <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
                     rowSums(seqtabAll), rowSums(seqtabNoC))
colnames(read_counts) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtab", "nochimera")
rownames(read_counts) <- sampleNames

# Identify the reference database you will use for assigning taxonomy. This was downloaded from https://benjjneb.github.io/dada2/training.html.
fastaRef <- "raw_data/dada2_rdp_train_set_16.fa.gz"

# Assign taxonomy to your sequence table.
taxtabNoC <- assignTaxonomy(seqtabNoC, refFasta=fastaRef, multithread=TRUE)
taxtab.wrong.length <- assignTaxonomy(seqtab.wrong.length.asvs, refFasta=fastaRef, multithread=TRUE)

# Check to see that it worked.
unname(head(taxtabNoC))

# Import into phyloseq object and clean workspace
pseq.raw <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE),
                     sample_data(sample_data),
                     tax_table(taxtabNoC))

pseq.wrong.length <- phyloseq(otu_table(seqtab.wrong.length.asvs, taxa_are_rows=FALSE),
                              sample_data(sample_data),
                              tax_table(taxtab.wrong.length))

rm(seqtabNoC, sample_data, taxtabNoC, fastaRef, seqtabAll, mergers, out, derepRs, derepFs, errF, errR,
   dadaFs, dadaRs, filt_path, input_path, fnFs, fnRs, sampleNames, seqtab.wrong.length.asvs, taxtab.wrong.length)

# Remove chloroplasts
temp <- as.data.frame(tax_table(pseq.raw))
temp$chloroplast <- temp$Class == "Chloroplast"
temp <- subset(temp, chloroplast=="TRUE")
badTaxa <- rownames(temp)
goodTaxa <- setdiff(taxa_names(pseq.raw), badTaxa)
pseq.raw <- prune_taxa(goodTaxa, pseq.raw)
rm(temp, goodTaxa, badTaxa)

temp <- as.data.frame(tax_table(pseq.wrong.length))
temp$chloroplast <- temp$Class == "Chloroplast"
temp <- subset(temp, chloroplast=="TRUE")
badTaxa <- rownames(temp)
goodTaxa <- setdiff(taxa_names(pseq.wrong.length), badTaxa)
pseq.wrong.length <- prune_taxa(goodTaxa, pseq.raw)
rm(temp, goodTaxa, badTaxa)

#### IDENTIFY ASVS NOT WITHIN 250-256bp (PART OF SPONGE SPECIES IDENTIFICATION) ####
# Rename ASVs for clarity
dna <- Biostrings::DNAStringSet(taxa_names(pseq.wrong.length))
names(dna) <- taxa_names(pseq.wrong.length)
pseq.wrong.length <- merge_phyloseq(pseq.wrong.length, dna)
rm(dna)
taxa_names(pseq.wrong.length) <- paste0("WrongLength", seq(ntaxa(pseq.wrong.length)))

# Determine which taxa are not assigned to the phylum level for use in a BLAST search.
temp.unclass <- subset_taxa(pseq.wrong.length, Phylum !="NA")
temp.unclass <- as.data.frame(cbind(tax_table(temp.unclass)))
class.string <- rownames(temp.unclass)
temp.unclass <- subset_taxa(pseq.wrong.length, !(taxa_names(pseq.wrong.length) %in% class.string))

# Create a data frame showing wrong length ASVs by ASV sequence and total abundance.
seqtab.nochim <- as.matrix(otu_table(temp.unclass))
seqs <- refseq(temp.unclass)
ids <- taxa_names(temp.unclass)
db_out <- data.frame(ids=ids, seqs=seqs, count=colSums(seqtab.nochim))
db_out <- subset(db_out, count > 0)
db_out <- db_out[order(-db_out$count), ] # Sort by abundance
db_out$Percent <- 100*db_out$count / sum(db_out$count) # Relative abundance (only among wrong-length ASVs)

# Export sequences to a FASTA file for use in a BLAST search.
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
writeFasta(fasta, file = "wrong_length_seqs.fna")
write.csv(db_out, "wrong_length_seqs.csv")

# Two ASVs (WrongLength1 and WrongLength2) were identified as E. muelleri sequences in the BLAST search.
# These two sequences were then used in a MAFFT alignment against other freshwater sponge 16S rRNA sequences.

# Import the results from the MAFFT alignment.
speciesid.amplicon.tree <- ape::read.tree("mafft_16S_rRNA_gene_alignment/output.tree")

# Look at the abundance of these two ASVs across all samples
sp.pseq.otu <- as.data.frame(cbind(otu_table(pseq.wrong.length)))
sp.pseq.otu <- sp.pseq.otu[,c("WrongLength1", "WrongLength2")]

# Add sample data to create a table showing each sample and the abundance of E. muelleri reads
speciesid.amplicon.abundance <- transform(merge(sample_data, sp.pseq.otu, by=0, all.x=TRUE, all.y=FALSE),
                                          row.names=Row.names, Row.names=NULL)
speciesid.amplicon.abundance <- speciesid.amplicon.abundance[,c("SampleID", "Type","Source","WrongLength1", "WrongLength2")]
speciesid.amplicon.abundance$Sum <- speciesid.amplicon.abundance$WrongLength1 + speciesid.amplicon.abundance$WrongLength2

# Clean workspace
rm(sp.pseq.otu, fasta, db_out, seqs, seqtab.nochim, class.string, temp.unclass, pseq.wrong.length)

#### MOCK COMMUNITY ANALYSIS ####
# Subset to mock community
mock.comm <- subset_samples(pseq.raw, sample_names(pseq.raw) %in% c("NCMC1"))
mock.comm <- prune_taxa(taxa_sums(mock.comm) > 0, mock.comm)

# Replace empty taxa names
mock.comm.sub <- mock.comm
for(i in 1:nrow(tax_table(mock.comm.sub))){
  for(j in 2:6){
    if(is.na(tax_table(mock.comm.sub)[i,j])==TRUE){
      if(substr(tax_table(mock.comm.sub)[i,j-1], 1, 4)=="Uncl"){
        tax_table(mock.comm.sub)[i,j] <- tax_table(mock.comm.sub)[i,j-1]}
      else {
        tax_table(mock.comm.sub)[i,j] <- paste0("Uncl_", tax_table(mock.comm.sub)[i,j-1])}}
  }}

# Extract OTU table and sequences
mock.comm.data <- as.data.frame(t(otu_table(transform_sample_counts(mock.comm.sub, function(x) 100*x/sum(x)))))
mock.comm.data <- cbind(as.data.frame(cbind(tax_table(mock.comm.sub))),
                        mock.comm.data)
mock.comm.data$Kingdom <- NULL

mock.comm.data <- mock.comm.data[order(-mock.comm.data$NCMC1), ]
mock.comm.data$Taxon <- paste0("MockCom", seq(nrow(mock.comm.data)))
levels(mock.comm.data$Genus)[levels(mock.comm.data$Genus)=="Escherichia/Shigella"] <- "Escherichia"
mock.comm.data$Taxon <- paste0(mock.comm.data$Taxon, ": ", mock.comm.data$Genus)

# Export sequence data for a MAFFT alignment and BLAST search.
fasta <- ShortRead(sread = DNAStringSet(rownames(mock.comm.data)),
                   id = BStringSet(mock.comm.data$Taxon))

writeFasta(fasta, file = "mock.community.fasta")
write.csv(mock.comm.data, "mock.community.csv")
rm(fasta)

# Run MAFFT alignment between mock community sequences and downloaded reference genomes
# Import results of MAFFT alignment.
mock.comm.tree <- ape::read.tree("~/Downloads/sponge_ISME_revisions/mafft_mock_community_alignment/output.tree")

# Visualize tree
ggtree(mock.comm.tree) + 
  geom_tiplab() + 
  xlim(0,0.5) + 
  labs(tag="a") + 
  plot_theme + 
  theme(panel.border=element_blank()) #+ geom_text(aes(label=node))

rm(mock.comm, mock.comm.sub, mock.comm.tree)

#### IDENTIFY AND REMOVE CONTAMINANTS ####
# Transition to new phyloseq object.
pseq.filter <- pseq.raw

# Drop mock community from subsequent analysis.
pseq.filter <- subset_samples(pseq.filter, !sample_names(pseq.filter) %in% c("NCMC1"))
pseq.filter <- prune_taxa(taxa_sums(pseq.filter) > 0, pseq.filter)
sample_data <- subset(sample_data, SampleID !="NCMC1")

# Store taxa names as a 'refseq' object and rename taxa with temporary names.
dna <- Biostrings::DNAStringSet(taxa_names(pseq.filter))
names(dna) <- taxa_names(pseq.filter)
pseq.filter <- merge_phyloseq(pseq.filter, dna)
rm(dna)
taxa_names(pseq.filter) <- paste0("temp", seq(ntaxa(pseq.filter)))

# Determine which taxa are not assigned to the phylum level for use in a BLAST search.
temp.unclass <- subset_taxa(pseq.filter, Phylum !="NA")
temp.unclass <- as.data.frame(cbind(tax_table(temp.unclass)))
class.string <- rownames(temp.unclass)
temp.unclass <- subset_taxa(pseq.filter, !(taxa_names(pseq.filter) %in% class.string))

# Export the ASV sequences for these taxa. BLAST search to identify closest sequence identity.
seqtab.nochim <- as.matrix(otu_table(temp.unclass))
seqs <- refseq(temp.unclass)
ids <- taxa_names(temp.unclass)
db_out <- data.frame(ids=ids, seqs=seqs, count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
writeFasta(fasta, file = "ASVs_with_no_phylum.fna")
write.csv(db_out, "ASVs_with_no_phylum.csv")
rm(fasta, db_out, ids, seqs, temp.unclass, seqtab.nochim, class.string)

# Test for archaeal presence/absence and abundance
archaea <- transform_sample_counts(pseq.filter, function(x) 100*x/sum(x))
archaea <- subset_taxa(archaea, Kingdom=="Archaea")
tax_table(archaea) # Only four archaeal ASVs
rowMeans(otu_table(archaea)) # All at extremely low abundances
pseq.filter <- subset_taxa(pseq.filter, !(Kingdom %in% c("Archaea"))) # Remove archaea from analysis
rm(archaea)

# Remove taxa not assigned to the phylum level
# The BLAST results (above) suggested that most of these were noise - chloroplasts, mitochondria, etc.
pseq.filter <- subset_taxa(pseq.filter, is.na(Phylum)=="FALSE")

# Identify contaminants using decontam.
sample_data(pseq.filter)$is.control <- sample_data(pseq.filter)$Type == "control" # Create logical identifying controls.
contamdf.prev <- isContaminant(pseq.filter, method="prevalence", neg="is.control", threshold=0.25)
table(contamdf.prev$contaminant)

# Plot presence/absence in negative controls vs. true samples
ctrl.contam.data <- list()

ps.pa <- transform_sample_counts(pseq.filter, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Type == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Type != "control", ps.pa)
ctrl.contam.data$prevalence <- data.frame(pa.pos=taxa_sums(ps.pa.pos),
                                     pa.neg=taxa_sums(ps.pa.neg),
                                     contaminant=contamdf.prev$contaminant)

ggplot(data=ctrl.contam.data$prevalence, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

rm(ps.pa, ps.pa.neg, ps.pa.pos)

# Check relative abundances of contaminants
contamdf.prev <- subset(contamdf.prev, contaminant==TRUE) # Subset to true contaminants
ctrl.contam.data$abundance <- transform_sample_counts(pseq.filter, function(x) 100*x/sum(x))
ctrl.contam.data$abundance <- subset_taxa(ctrl.contam.data$abundance, 
                                     taxa_names(ctrl.contam.data$abundance) %in% as.character(rownames(contamdf.prev))) # Select only contaminants
ctrl.contam.data$abundance <- cbind(as.data.frame(cbind(tax_table(ctrl.contam.data$abundance))),
                               as.data.frame(t(otu_table(ctrl.contam.data$abundance))))
ctrl.contam.data[["abundance"]]$means <- rowMeans(ctrl.contam.data$abundance[,c(7:49)])
ctrl.contam.data[["abundance"]]$ASV <- rownames(ctrl.contam.data[["abundance"]]) 

ggplot(ctrl.contam.data$abundance, aes(x=ASV, y=means)) + geom_bar(stat="identity", position="dodge") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# How many reads in each control sample?
ctrl.nc.data <- list()
sample_sums(pseq.filter)[c("K1", "K2", "NCFILT1")]

pseq.nc <- subset_samples(pseq.filter, sample_names(pseq.filter) %in% c("K1", "K2", "NCFILT1"))
pseq.nc <- prune_taxa(taxa_sums(pseq.nc) > 0, pseq.nc)
for(i in 1:nrow(tax_table(pseq.nc))){
  for(j in 2:6){
    if(is.na(tax_table(pseq.nc)[i,j])==TRUE){
      if(substr(tax_table(pseq.nc)[i,j-1], 1, 4)=="Uncl"){
        tax_table(pseq.nc)[i,j] <- tax_table(pseq.nc)[i,j-1]}
      else {
        tax_table(pseq.nc)[i,j] <- paste0("Uncl_", tax_table(pseq.nc)[i,j-1])}}
  }}
pseq.nc <- tax_glom(pseq.nc, taxrank="Genus", NArm=FALSE)
pseq.nc <- transform_sample_counts(pseq.nc, function(x) 100*x/sum(x))
taxa_names(pseq.nc) <- tax_table(pseq.nc)[,"Genus"]

ctrl.nc.data$genus.abundances <- as.data.frame(t(otu_table(pseq.nc)))

pseq.nc <- subset_samples(pseq.filter, !sample_names(pseq.filter) %in% c("K0SKE2", "K0SKE3", "K0SKE5"))
pseq.nc <- prune_taxa(taxa_sums(pseq.nc) > 0, pseq.nc)
temp <- sample_data(pseq.nc)
levels(temp$Type)[levels(temp$Type) %in% c("sponge", "water", "biofilm")] <- "experiment"
sample_data(pseq.nc) <- temp
pseq.nc <- transform_sample_counts(pseq.nc, function(x) 100*x/sum(x))
pseq.nc.otu <- as.data.frame(cbind(temp$Type, otu_table(pseq.nc)))
pseq.nc.otu <- pseq.nc.otu %>% dplyr::group_by(V1) %>% dplyr::summarise_if(is.numeric, mean) %>% as.data.frame()
rownames(pseq.nc.otu) <- c("experiment", "negative control")
pseq.nc.otu$V1 <- NULL
pseq.nc.otu <- as.data.frame(t(pseq.nc.otu))

ggplot(pseq.nc.otu, aes(x=experiment, y=`negative control`)) + geom_point() + scale_y_log10() + scale_x_log10()

ctrl.nc.data$abundance.corr <- pseq.nc.otu
rm(pseq.nc, pseq.nc.otu, temp)

# Remove contaminants from samples
pseq.filter <- subset_taxa(pseq.filter, !(taxa_names(pseq.filter) %in% as.character(rownames(contamdf.prev)))) # Remove 30 contaminants.
sample_data(pseq.filter)$is.control <- NULL
rm(contamdf.prev)

# Remove negative controls from phyloseq object and sample data.
pseq.filter <- subset_samples(pseq.filter, Type !="control") # Remove control samples from further analysis.
pseq.filter <- prune_taxa(taxa_sums(pseq.filter) > 0, pseq.filter)

sample_data <- subset(sample_data, !SampleID %in% c("NCMC1", "K1", "K2", "NCFILT1"))

#### SEQUENCING REPLICABILITY ANALYSIS ####
# Create a phyloseq object containing only the original three samples and their replicates.
ctrl.replicates <- subset_samples(pseq.filter, sample_names(pseq.filter) %in% c("SKE2", "SKE3", "SKE5",
                                                                                     "K0SKE2", "K0SKE3", "K0SKE5"))
ctrl.replicates <- prune_taxa(taxa_sums(ctrl.replicates) > 0, ctrl.replicates)
sample_sums(ctrl.replicates)

# Create a data frame where each column is a list of the ASVs present in each sample.
Venn <- data.frame(ncol=1)
for(i in (1:length(sample_names(ctrl.replicates)))){
  a <- subset_samples(ctrl.replicates, SampleID==sample_names(ctrl.replicates)[[i]])
  a <- prune_taxa(taxa_sums(a) > 0, a)
  b <- as.data.frame(rownames(tax_table(a)))
  colnames(b) <- sample_names(ctrl.replicates)[[i]]
  Venn <- rowr::cbind.fill(Venn, b, fill=NA)
  rm(a, b)
}
Venn$ncol <- NULL

#### Calculate abundance-weighted OTU overlap for each pair of technical replicates.
temp <- subset_taxa(ctrl.replicates, taxa_names(ctrl.replicates) %in% intersect(Venn$SKE2, Venn$K0SKE2))
temp <- subset_samples(temp, SampleID=="SKE2" | SampleID=="K0SKE2")
temp <- as.data.frame(otu_table(temp))
(sum(temp[1,])+sum(temp[2,]))/sum(sample_sums(ctrl.replicates)["SKE2"]+sample_sums(ctrl.replicates)["K0SKE2"]) # 0.9291981

temp <- subset_taxa(ctrl.replicates, taxa_names(ctrl.replicates) %in% intersect(Venn$SKE3, Venn$K0SKE3))
temp <- subset_samples(temp, SampleID=="SKE3" | SampleID=="K0SKE3")
temp <- as.data.frame(otu_table(temp))
(sum(temp[1,])+sum(temp[2,]))/sum(sample_sums(ctrl.replicates)["SKE3"]+sample_sums(ctrl.replicates)["K0SKE3"]) # 0.9522301

temp <- subset_taxa(ctrl.replicates, taxa_names(ctrl.replicates) %in% intersect(Venn$SKE5, Venn$K0SKE5))
temp <- subset_samples(temp, SampleID=="SKE5" | SampleID=="K0SKE5")
temp <- as.data.frame(otu_table(temp))
(sum(temp[1,])+sum(temp[2,]))/sum(sample_sums(ctrl.replicates)["SKE5"]+sample_sums(ctrl.replicates)["K0SKE5"]) # 0.9369973

rm(Venn, temp)

#### Check for differences in alpha diversity.
temp <- cbind(sample_data(ctrl.replicates), estimate_richness(ctrl.replicates))

car::leveneTest(Observed ~ Type, temp) # p = 0.5582
t.test(subset(temp, Type=="replicate")$Observed,
       subset(temp, Type=="sponge")$Observed,
       paired=TRUE, var.equal=TRUE) # t=1.2908, d=2, p=0.3259

car::leveneTest(Shannon ~ Type, temp) # p = 0.1612
t.test(subset(temp, Type=="replicate")$Shannon,
       subset(temp, Type=="sponge")$Shannon,
       paired=TRUE, var.equal=TRUE) # t=1.4644, df=2, p=0.2807

rm(temp)

#### Check for differences in beta diversity using 1) bdiv.ordinations.all and 2) between- and among-sample distances
ctrl.replicate.data <- list()

for(i in c(1:2)){
  # Do this both using the standard Bray-Curtis approach on rarefied data (i==1) 
  # and the CLR transform approach (i==2)
  if(i==1){
    # Rarefy
    temp.otu <- as.data.frame(otu_table(rarefy_even_depth(ctrl.replicates,
                                                          sample.size=min(sample_sums(ctrl.replicates)),
                                                          rngseed=700,
                                                          replace=FALSE,
                                                          trimOTUs=TRUE)))
    # Calculate a distance matrix
    temp.dist <- vegdist(temp.otu, method="bray")
    
    # Calculate an ordination
    temp.pca <- capscale(temp.dist ~ 1)
  }
  
  if(i==2){
    # Convert to relative abundance and remove low-abundance taxa from the original phyloseq object
    x1 <- transform_sample_counts(pseq.filter, function(x) 100*x/sum(x))
    x2 <- as.data.frame(colMeans(otu_table(x1)))
    x1 <- filter_taxa(x1, function(x) mean(x) > 0.001, TRUE)
    temp.otu <- subset_taxa(ctrl.replicates, taxa_names(ctrl.replicates) %in% taxa_names(x1))
    rm(x1, x2) 
    
    # Perform the CLR transform, calculate a distance metric, and perform ordination.
    temp.otu <- as.data.frame(otu_table(microbiome::transform(ctrl.replicates, "clr")))
    temp.dist <- vegdist(temp.otu, method="euclidean")
    temp.pca <- rda(temp.otu)
  }
  
  # Extract axis scores from each ordination.
  temp.scores <- cbind(sample_data(ctrl.replicates),
                       vegan::scores(temp.pca, choices=c(1:2), display="sites"))
  colnames(temp.scores)[c(4:5)] <- c("PC1", "PC2")
  
  # Save the ordination as a plot
  plot <- ggplot(temp.scores, aes(x=PC1, y=PC2)) +
    geom_point(aes(color=Type)) +
    geom_label(aes(label=SampleID))
  
  # PERMANOVA between sample groups
  test <- adonis(temp.dist ~ Type,
                 as.data.frame(cbind(sample_data(ctrl.replicates))), permutations=999)
  
  # Compare distances betwene paired samples to distances between non-paires
  temp.dist <- as.matrix(temp.dist)
  
  # Distances betwen pairs.
  vector1 <- c(temp.dist["SKE2","K0SKE2"], temp.dist["SKE3","K0SKE3"], temp.dist["SKE5","K0SKE5"]) # Between pairs.
  
  # Distances between non-pairs.
  temp.dist[,"K0SKE2"] <- rep(NA)
  temp.dist[,"K0SKE3"] <- rep(NA)
  temp.dist[,"K0SKE5"] <- rep(NA)
  temp.dist[temp.dist==0] <- NA
  vector2 <- rowMeans(temp.dist, na.rm=TRUE)[c(1:3)]
  
  # Assemble into data frame.
  temp <- data.frame(Sample=c("SKE2", "SKE3", "SKE5"),
                     Between=vector1,
                     Among=vector2)
  
  # Save data.
  ctrl.replicate.data[[i]] <- list()
  ctrl.replicate.data[[i]][[1]] <- plot
  ctrl.replicate.data[[i]][[2]] <- test
  ctrl.replicate.data[[i]][[3]] <- temp
  
  names(ctrl.replicate.data[[i]]) <- c("plot", "adonis", "distance")
  
  rm(temp.otu, temp.dist, temp.pca, temp.scores,
     plot, test, vector1, vector2, temp)
}

#### Check for differences in general taxon abundances between individual samples and their replicates.
temp.replicate <- tax_glom(ctrl.replicates, taxrank="Genus") # Do this at the genus level.

# Create separate phyloseq objects for each pair of samples.
replicate1 <- subset_samples(temp.replicate, sample_names(temp.replicate) %in% c("SKE2", "K0SKE2"))
replicate2 <- subset_samples(temp.replicate, sample_names(temp.replicate) %in% c("SKE3", "K0SKE3"))
replicate3 <- subset_samples(temp.replicate, sample_names(temp.replicate) %in% c("SKE5", "K0SKE5"))

replicate1 <- prune_taxa(taxa_sums(replicate1) > 0, replicate1)
replicate2 <- prune_taxa(taxa_sums(replicate2) > 0, replicate2)
replicate3 <- prune_taxa(taxa_sums(replicate3) > 0, replicate3)

# Convert to data frame showing relative abundance.
replicate1 <- as.data.frame(t(otu_table(transform_sample_counts(replicate1, function(x) 100*x/sum(x)))))
replicate2 <- as.data.frame(t(otu_table(transform_sample_counts(replicate2, function(x) 100*x/sum(x)))))
replicate3 <- as.data.frame(t(otu_table(transform_sample_counts(replicate3, function(x) 100*x/sum(x)))))

# Plot relative abundance of taxa in original sample (x-axis) vs. relative abundance in replicate (y-axis).
grid.arrange(
  ggplot(replicate1, aes(x=SKE2, y=K0SKE2)) + geom_point() + scale_y_log10() + scale_x_log10(),
  ggplot(replicate2, aes(x=SKE3, y=K0SKE3)) + geom_point() + scale_y_log10() + scale_x_log10(),
  ggplot(replicate3, aes(x=SKE5, y=K0SKE5)) + geom_point() + scale_y_log10() + scale_x_log10(),
  ncol=3)

# Test correlations between originals and replicates.
cor.test(replicate1$SKE2, replicate1$K0SKE2, method = "spearman") # p < 0.001
cor.test(replicate2$SKE3, replicate2$K0SKE3, method = "spearman") # p < 0.001
cor.test(replicate3$SKE5, replicate3$K0SKE5, method = "spearman") # p < 0.001

rm(replicate1, replicate2, replicate3, temp.replicate)

# Remove replicates from analysis and move on.
pseq.filter <- subset_samples(pseq.filter, Type !="replicate")
pseq.filter <- prune_taxa(taxa_sums(pseq.filter) > 0, pseq.filter)

sample_data <- subset(sample_data, !SampleID %in% c("K0SKE2", "K0SKE3", "K0SKE5"))

#### DATA PRE-PROCESSING (RAREFACTION, TAXON FILTERING, ETC.) ####
sample_data <- subset(sample_data, Type !="replicate" & Type !="control")

# Examine read count distribution to determine appropriate cutoff threshold.
temp <- as.data.frame(sample_sums(pseq.filter))
ggplot(temp, aes(x=temp[,1])) + geom_histogram() # No need to remove any samples.
rm(temp)

# Calculate read profiles for each original sample
df <- as.data.frame(sample_sums(pseq.new))
df <- subset(df, !(rownames(df) %in% c("NCFILT1", "K1", "K2", "K0SKE2", "K0SKE3", "K0SKE5", "NCMC1")))
min(df[,1])
max(df[,1])
mean(df[,1])
sd(df[,1])
rm(df)

# Rarefy samples to the minimum remaining library size (8000 reads). Average across 1,000 rarefactions.
taxa_names(pseq.filter) <- paste0("ASV", seq(ntaxa(pseq.filter)))

pseq.clr <- pseq.filter
pseq.rarefied <- pseq.filter

rarefaction.average <- list()
for(i in 1:1000){
  temp.rarefy <- rarefy_even_depth(pseq.filter, 
                                   sample.size = min(sample_sums(pseq.filter)),
                                   replace = FALSE,
                                   trimOTUs = FALSE,
                                   rngseed=i)
  rarefaction.average[[i]] <- as.data.frame(otu_table(temp.rarefy))
}

dfAvg <- Reduce("+", rarefaction.average)/length(rarefaction.average)
dfAvg <- round(dfAvg, 0)
dfAvg <- otu_table(dfAvg, taxa_are_rows=FALSE)

# Replace feature table in phyloseq object with new, rarefied feature table.
otu_table(pseq.rarefied) <- dfAvg
pseq.rarefied <- prune_taxa(taxa_sums(pseq.rarefied) > 0, pseq.rarefied)
rm(temp.rarefy, rarefaction.average, dfAvg)

# Construct a phylogenetic tree for both filtered and rarefied data, for use in a phyloseq object.
pseq.list <- list(pseq.filter, pseq.rarefied)

for(i in c(1:2)){
  vector <- taxa_names(pseq.list[[i]])
  taxa_names(pseq.list[[i]]) <- refseq(pseq.list[[i]])
  seqtabNoC <- as.matrix(as.data.frame(otu_table(pseq.list[[i]])))
  seqs <- getSequences(seqtabNoC)
  
  names(seqs) <- seqs
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose=FALSE)
  phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phangAlign)
  treeNJ <- NJ(dm)
  fit = pml(treeNJ, data=phangAlign)
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement="stochastic", control=pml.control(trace=0))
  
  phy_tree(pseq.list[[i]]) <- fitGTR$tree
  taxa_names(pseq.list[[i]]) <- vector
  
  rm(fitGTR, fit, treeNJ, dm, phangAlign, alignment, seqs, seqtabNoC, vector)
}

pseq.clr <- pseq.list[[1]]
pseq.rarefied <- pseq.list[[2]]
rm(pseq.list)

# Replace 'NA' taxa with 'Uncl' to prevent errors during agglomeration.
pseq.rarefied.sub <- pseq.rarefied
for(i in 1:nrow(tax_table(pseq.rarefied.sub))){
  for(j in 2:6){
    if(is.na(tax_table(pseq.rarefied.sub)[i,j])==TRUE){
      if(substr(tax_table(pseq.rarefied.sub)[i,j-1], 1, 4)=="Uncl"){
        tax_table(pseq.rarefied.sub)[i,j] <- tax_table(pseq.rarefied.sub)[i,j-1]}
      else {
        tax_table(pseq.rarefied.sub)[i,j] <- paste0("Uncl_", tax_table(pseq.rarefied.sub)[i,j-1])}}
  }}

pseq.clr.sub <- pseq.clr
for(i in 1:nrow(tax_table(pseq.clr.sub))){
  for(j in 2:6){
    if(is.na(tax_table(pseq.clr.sub)[i,j])==TRUE){
      if(substr(tax_table(pseq.clr.sub)[i,j-1], 1, 4)=="Uncl"){
        tax_table(pseq.clr.sub)[i,j] <- tax_table(pseq.clr.sub)[i,j-1]}
      else {
        tax_table(pseq.clr.sub)[i,j] <- paste0("Uncl_", tax_table(pseq.clr.sub)[i,j-1])}}
  }}

# Calculate alpha diversity metrics from rarefied data.
alpha.data <- estimate_richness(pseq.rarefied) # Alpha diversity.
alpha.data$SampleID <- rownames(alpha.data)
sample_data <- merge(sample_data, alpha.data, by="SampleID", all=TRUE)
rm(alpha.data)

# Calculate phylogenetic diversity using picante
temp.phylo.class = phy_tree(pseq.rarefied)
temp.vegan.otu.table = as.data.frame(otu_table(pseq.rarefied))
picante.pd.result <- pd(temp.vegan.otu.table, temp.phylo.class, include.root=FALSE)

picante.results <- data.frame(SampleID=rownames(picante.pd.result),
                              PD = picante.pd.result$PD)

sample_data <- merge(sample_data, picante.results, by="SampleID", all=TRUE)
rm(picante.results, temp.phylo.class, temp.vegan.otu.table, picante.pd.result)


# Calculate richness and diversity using rarefaction and extrapolation from unrarefied data (iNext)
temp <- as.data.frame(t(otu_table(pseq.clr)))
iNext.result <- iNEXT(temp, q=0, datatype="abundance") # All samples, except the five low-read samples, have completion >98%
rm(temp)

temp <- iNext.result$AsyEst
temp.richness <- subset(temp, Diversity=="Species richness")
temp.richness$SampleID <- temp.richness$Site
temp.richness$Observed_Extrap <- temp.richness$Estimator
temp.richness$Observed_Extrap <- round(temp.richness$Observed_Extrap, digits = 0)

temp.diversity <- subset(temp, Diversity=="Shannon diversity")
temp.diversity$SampleID <- temp.diversity$Site
temp.diversity$Shannon_Extrap <- log(temp.diversity$Estimator)

sample_data <- merge(sample_data, temp.richness[,c("SampleID", "Observed_Extrap")], by="SampleID", all=TRUE)
sample_data <- merge(sample_data, temp.diversity[,c("SampleID", "Shannon_Extrap")], by="SampleID", all=TRUE)

rm(temp.diversity, temp.richness, temp)

# Calculate distance matrices, as appropriate.
dist.sp.list <- list()
dist.sp.list[["Bray-Curtis"]] <- vegdist(as.data.frame(otu_table(pseq.rarefied)), method="bray")
dist.sp.list[["Jaccard"]] <- vegdist(as.data.frame(otu_table(pseq.rarefied)), method="jaccard")
dist.sp.list[["Aitchison"]] <- vegdist(as.data.frame(otu_table(microbiome::transform(pseq.clr, "clr"))), method="euclidean")
dist.sp.list[["wUF"]] <- UniFrac(pseq.rarefied, weighted=TRUE, normalized=TRUE)
dist.sp.list[["uwUF"]] <- UniFrac(pseq.rarefied, weighted=FALSE, normalized=TRUE)

temp.rarefied <- subset_samples(pseq.rarefied, Type=="sponge")
temp.clr <- subset_samples(pseq.clr, Type=="sponge")
temp.rarefied <- prune_taxa(taxa_sums(temp.rarefied) > 0, temp.rarefied)
temp.clr <- prune_taxa(taxa_sums(temp.clr) > 0, temp.clr)

dist.sp.list.sponge <- list()
dist.sp.list.sponge[["Bray-Curtis"]] <- vegdist(as.data.frame(otu_table(temp.rarefied)), method="bray")
dist.sp.list.sponge[["Jaccard"]] <- vegdist(as.data.frame(otu_table(temp.rarefied)), method="jaccard")
dist.sp.list.sponge[["Aitchison"]] <- vegdist(as.data.frame(otu_table(microbiome::transform(temp.clr, "clr"))), method="euclidean")
dist.sp.list.sponge[["wUF"]] <- UniFrac(temp.rarefied, weighted=TRUE, normalized=TRUE)
dist.sp.list.sponge[["uwUF"]] <- UniFrac(temp.rarefied, weighted=FALSE, normalized=TRUE)

rm(temp.rarefied, temp.clr)

# Create a summary of all taxon abundances and prevalences, for subsetting purposes later.
prev.data <- list()
taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")
for(i in c(1:6)){
  temp.physeq <- pseq.clr.sub
  if(i<6){temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])}
  if(i<6){taxa_names(temp.physeq) <- as.data.frame(cbind(tax_table(temp.physeq)))[,i+1]}
  prev <- as.data.frame(microbiome::prevalence(temp.physeq, detection=0, count=TRUE))
  colnames(prev) <- "Prevalence"
  prev$Taxon <- rownames(prev)
  temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp <- as.data.frame(otu_table(temp.physeq))
  prev$Abundance <- colMeans(temp)
  
  # Mean relative abundances by type
  temp$Type <- sample_data$Type
  temp2 <- temp %>% dplyr::group_by(Type) %>% dplyr::summarise_if(is.numeric, mean) %>% as.data.frame()
  rownames(temp2) <- temp2$Type
  temp2$Type <- NULL
  temp2 <- as.data.frame(t(temp2))
  prev <- cbind(prev, temp2)
  
  # Mean relative abundance by source/type
  temp$Merge <- paste0(sample_data$Source, "_", sample_data$Type)
  temp$Type <- NULL
  temp$Source <- NULL
  temp2 <- temp %>% dplyr::group_by(Merge) %>% dplyr::summarise_if(is.numeric, mean) %>% as.data.frame()
  rownames(temp2) <- temp2$Merge
  temp2$Merge <- NULL
  temp2 <- as.data.frame(t(temp2))
  prev <- cbind(prev, temp2)
  
  prev.data[[i]] <- prev
  rm(prev, temp.physeq, temp, temp2)
}
names(prev.data) <- c(taxranks, "ASV")

##################################################### 16S rRNA AMPLICON ANALYSIS ###########################################
#### CALCULATE RAREFACTION CURVES ##########
# Define rarefaction curve calculation function (obtained from https://github.com/joey711/phyloseq/issues/143)
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

# Calculate rarefaction curves.
rarefaction_curve_data <- calculate_rarefaction_curves(pseq.clr, c('Observed', 'Shannon'), rep(c(1, 10, 100, seq(1000, 37000, by=1000)), each = 10))

# Summarize alpha diversity and add sample data.
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary$Sample <- gsub("X", "", rarefaction_curve_data_summary$Sample)
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(pseq.rarefied)), by.x = 'Sample', by.y = 'row.names')
rm(rarefaction_curve_data, rarefaction_curve_data_summary)

rarefaction_curve_data <- rarefaction_curve_data_summary_verbose
rm(rarefaction_curve_data_summary_verbose)

ggplot(subset(rarefaction_curve_data, Measure=="Observed"), aes(x=Depth, y=Alpha_diversity_mean)) + 
  geom_point() + geom_line(aes(group=SampleID)) + 
  facet_wrap(~Type)

#### ALPHA DIVERSITY ANALYSES: SIGNIFICANT DIFFERENCES AMONG GROUPINGS #########
# Among sample types (all rivers pooled)
adiv.signif.all <- data.frame(matrix(nrow=0, ncol=8))
adiv.signif.sponge <- data.frame(matrix(nrow=0, ncol=8))
adiv.signif.water <- data.frame(matrix(nrow=0, ncol=8))
adiv.signif.biofilm <- data.frame(matrix(nrow=0, ncol=8))
adiv.signif.Sooke <- data.frame(matrix(nrow=0, ncol=8))
adiv.signif.Nanaimo <- data.frame(matrix(nrow=0, ncol=8))
adiv.signif.Cowichan <- data.frame(matrix(nrow=0, ncol=8))

for(i in c("Observed", "Observed_Extrap", "Chao1", "Shannon", "Shannon_Extrap", "PD")){
  
  # Comparisons among sample types (all rivers pooled)
  levene <- car::leveneTest(sample_data[,i] ~ Type, sample_data)
  anova <- summary(aov(sample_data[,i] ~ Type, sample_data))
  tukey <- TukeyHSD(aov(sample_data[,i] ~ Type, sample_data))
  
  adiv.signif.all[i,1] <- i # variable name
  adiv.signif.all[i,2] <- levene$`Pr(>F)`[1] # Levene's p
  adiv.signif.all[i,3] <- anova[[1]][1,4] # ANOVA F
  adiv.signif.all[i,4] <- anova[[1]][1,1] # ANOVA df
  adiv.signif.all[i,5] <- anova[[1]][1,5] # ANOVA p
  adiv.signif.all[i,6] <- tukey[[1]][1,4] # Tukey 1 water-sponge
  adiv.signif.all[i,7] <- tukey[[1]][2,4] # Tukey 2 biofilm-sponge
  adiv.signif.all[i,8] <- tukey[[1]][3,4] # Tukey 3 biofilm-water
  
  rm(levene, anova, tukey)
  
  # Comparisons among sponges (separated by river)
  test_data <- subset(sample_data, Type=="sponge")
  
  levene <- car::leveneTest(test_data[,i] ~ Source, test_data)
  anova <- summary(aov(test_data[,i] ~ Source, test_data))
  tukey <- TukeyHSD(aov(test_data[,i] ~ Source, test_data))
  
  adiv.signif.sponge[i,1] <- i # variable name
  adiv.signif.sponge[i,2] <- levene$`Pr(>F)`[1] # Levene's p
  adiv.signif.sponge[i,3] <- anova[[1]][1,4] # ANOVA F
  adiv.signif.sponge[i,4] <- anova[[1]][1,1] # ANOVA df
  adiv.signif.sponge[i,5] <- anova[[1]][1,5] # ANOVA p
  adiv.signif.sponge[i,6] <- tukey[[1]][1,4] # Tukey 1
  adiv.signif.sponge[i,7] <- tukey[[1]][2,4] # Tukey 2
  adiv.signif.sponge[i,8] <- tukey[[1]][3,4] # Tukey 3
  
  rm(levene, anova, tukey, test_data)
  
  # Comparisons across water (separated by river)
  test_data <- subset(sample_data, Type=="water")
  
  levene <- car::leveneTest(test_data[,i] ~ Source, test_data)
  anova <- summary(aov(test_data[,i] ~ Source, test_data))
  tukey <- TukeyHSD(aov(test_data[,i] ~ Source, test_data))
  
  adiv.signif.water[i,1] <- i # variable name
  adiv.signif.water[i,2] <- levene$`Pr(>F)`[1] # Levene's p
  adiv.signif.water[i,3] <- anova[[1]][1,4] # ANOVA F
  adiv.signif.water[i,4] <- anova[[1]][1,1] # ANOVA df
  adiv.signif.water[i,5] <- anova[[1]][1,5] # ANOVA p
  adiv.signif.water[i,6] <- tukey[[1]][1,4] # Tukey 1 Nanaimo-Sooke
  adiv.signif.water[i,7] <- tukey[[1]][2,4] # Tukey 2 Cowichan-Sooke
  adiv.signif.water[i,8] <- tukey[[1]][3,4] # Tukey 3 Cowichan-Nanaimo
  
  rm(levene, anova, tukey, test_data)
  
  # Comparisons across biofilms (separated by river). The Cowichan river must be excluded due to only one sample.
  test_data <- subset(sample_data, Type=="biofilm" & Source != "Cowichan")
  
  levene <- car::leveneTest(test_data[,i] ~ Source, test_data)
  anova <- t.test(test_data[,i] ~ Source, test_data, var.equal=TRUE)
  
  adiv.signif.biofilm[i,1] <- i # variable name
  adiv.signif.biofilm[i,2] <- levene$`Pr(>F)`[1] # Levene's p
  adiv.signif.biofilm[i,3] <- anova$statistic[1]
  adiv.signif.biofilm[i,4] <- anova$parameter[1]
  adiv.signif.biofilm[i,5] <- anova$p.value[1]
  
  rm(levene, anova, test_data)
  
  # Comparisons across sample types within the Sooke River.
  test_data <- subset(sample_data, Source=="Sooke")
  
  levene <- car::leveneTest(test_data[,i] ~ Type, test_data)
  anova <- summary(aov(test_data[,i] ~ Type, test_data))
  tukey <- TukeyHSD(aov(test_data[,i] ~ Type, test_data))
  
  adiv.signif.Sooke[i,1] <- i # variable name
  adiv.signif.Sooke[i,2] <- levene$`Pr(>F)`[1] # Levene's p
  adiv.signif.Sooke[i,3] <- anova[[1]][1,4] # ANOVA F
  adiv.signif.Sooke[i,4] <- anova[[1]][1,1] # ANOVA df
  adiv.signif.Sooke[i,5] <- anova[[1]][1,5] # ANOVA p
  adiv.signif.Sooke[i,6] <- tukey[[1]][1,4] # Tukey 1
  adiv.signif.Sooke[i,7] <- tukey[[1]][2,4] # Tukey 2
  adiv.signif.Sooke[i,8] <- tukey[[1]][3,4] # Tukey 3
  
  rm(levene, anova, tukey, test_data)
  
  # Comparisons across sample types within the Nanaimo River.
  test_data <- subset(sample_data, Source=="Nanaimo")
  
  levene <- car::leveneTest(test_data[,i] ~ Type, test_data)
  anova <- summary(aov(test_data[,i] ~ Type, test_data))
  tukey <- TukeyHSD(aov(test_data[,i] ~ Type, test_data))
  
  adiv.signif.Nanaimo[i,1] <- i # variable name
  adiv.signif.Nanaimo[i,2] <- levene$`Pr(>F)`[1] # Levene's p
  adiv.signif.Nanaimo[i,3] <- anova[[1]][1,4] # ANOVA F
  adiv.signif.Nanaimo[i,4] <- anova[[1]][1,1] # ANOVA df
  adiv.signif.Nanaimo[i,5] <- anova[[1]][1,5] # ANOVA p
  adiv.signif.Nanaimo[i,6] <- tukey[[1]][1,4] # Tukey 1 # water-sponge
  adiv.signif.Nanaimo[i,7] <- tukey[[1]][2,4] # Tukey 2 # biofilm-sponge
  adiv.signif.Nanaimo[i,8] <- tukey[[1]][3,4] # Tukey 3 # biofilm-water
  
  rm(levene, anova, tukey, test_data)
  
  # Comparisons across sample types within the Cowichan River (biofilms excluded)
  test_data <- subset(sample_data, Type != "biofilm" & Source=="Cowichan")
  
  levene <- car::leveneTest(test_data[,i] ~ Type, test_data)
  anova <- t.test(test_data[,i] ~ Type, test_data, var.equal=TRUE)
  
  adiv.signif.Cowichan[i,1] <- i # variable name
  adiv.signif.Cowichan[i,2] <- levene$`Pr(>F)`[1] # Levene's p
  adiv.signif.Cowichan[i,3] <- anova$statistic[1]
  adiv.signif.Cowichan[i,4] <- anova$parameter[1]
  adiv.signif.Cowichan[i,5] <- anova$p.value[1]
  
  rm(levene, anova, test_data)
  
}

# Calculate pairwise comparisons
temp_data <- sample_data
temp_data$var <- paste0(temp_data$Type, temp_data$Source)
adiv.tukey.HSD <- as.data.frame(as.matrix(TukeyHSD(aov(Observed_Extrap ~ var, temp_data))[[1]]))
rm(temp_data)

# Combine the data into a single list.
colnames(adiv.signif.all) <- c("Variable", "Levene.p", "Anova.F", "Anova.df", "Anova.p", "sp-wat", "sp-bio", "wat-bio")

colnames(adiv.signif.sponge) <- c("Variable", "Levene.p", "Anova.F", "Anova.df", "Anova.p", "SKE-NMO", "SKE-COW", "NMO-COW")
colnames(adiv.signif.water) <- c("Variable", "Levene.p", "Anova.F", "Anova.df", "Anova.p", "SKE-NMO", "SKE-COW", "NMO-COW")
colnames(adiv.signif.biofilm) <- c("Variable", "Levene.p", "Student.t", "Student.df", "Student.p")

colnames(adiv.signif.Sooke) <- c("Variable", "Levene.p", "Anova.F", "Anova.df", "Anova.p", "sp-wat", "sp-bio", "wat-bio")
colnames(adiv.signif.Nanaimo) <- c("Variable", "Levene.p", "Anova.F", "Anova.df", "Anova.p", "sp-wat", "sp-bio", "wat-bio")
colnames(adiv.signif.Cowichan) <- c("Variable", "Levene.p", "Student.t", "Student.df", "Student.p")

adiv.signif.tests <- list(all = adiv.signif.all,
                          sponge = adiv.signif.sponge,
                          water = adiv.signif.water,
                          biofilm = adiv.signif.biofilm,
                          Sooke = adiv.signif.Sooke,
                          Nanaimo = adiv.signif.Nanaimo,
                          Cowichan = adiv.signif.Cowichan)
rm(adiv.signif.all, adiv.signif.sponge, adiv.signif.water, adiv.signif.biofilm,
   adiv.signif.Sooke, adiv.signif.Nanaimo, adiv.signif.Cowichan)

adiv.summary <- sample_data %>% dplyr::group_by(Type) %>% dplyr::summarise_if(is.numeric, .funs=list(mean)) %>% as.data.frame()
temp <- sample_data %>% dplyr::group_by(Source, Type) %>% dplyr::summarise_if(is.numeric, .funs=list(mean)) %>% as.data.frame()
temp$Type <- c("SKE-sp", "SKE-wat", "SKE-bio", "NMO-sp", "NMO-wat", "NMO-bio", "COW-sp", "COW-wat", "COW-bio")
temp$Source <- NULL
adiv.summary <- rbind(adiv.summary, temp)
rm(temp)
rownames(adiv.summary) <- adiv.summary$Type
adiv.summary$Type <- NULL
adiv.summary <- as.data.frame(t(adiv.summary))

# Quick visual comparison
temp.melt <- melt(sample_data)
ggplot(subset(temp.melt, variable %in% c("Observed", "Observed_Extrap", "Chao1", "Shannon", "Shannon_Extrap", "PD")),
       aes(x=Source, y=value, fill=Type)) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free_y") +
  scale_fill_manual(values=c("darkorange1","cornflowerblue","chartreuse3"))

rm(temp.melt)

#### BETA DIVERSITY: CLUSTERING ANALYSIS ##########
bdiv.ordinations.all <- list()
bdiv.ordinations.all[["models"]] <- list()
bdiv.ordinations.all[["points"]] <- list()
bdiv.ordinations.all[["ellipses"]] <- list()

bdiv.ordinations.sp <- list()
bdiv.ordinations.sp[["models"]] <- list()
bdiv.ordinations.sp[["points"]] <- list()
bdiv.ordinations.sp[["ellipses"]] <- list()

veganCovEllipse <-function (cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Calculate PCoA bdiv.ordinations.all for the Bray-Curtis, Jaccard, weighted and unweighted UniFrac distances.
for(i in c(1,2,4,5)){
  # Run the PCoA model for each of the non-Aitchison distances.
  bdiv.ordinations.all[["models"]][[i]] <- capscale(dist.sp.list[[i]]~1)
  bdiv.ordinations.sp[["models"]][[i]] <- capscale(dist.sp.list.sponge[[i]]~1)
  
  # Extract the first two axes loadings.
  bdiv.ordinations.all[["points"]][[i]] = merge(scores(bdiv.ordinations.all[["models"]][[i]],
                                                       display="sites", choices=c(1:2)),
                                                sample_data[,c("SampleID", "Type","Source")],
                                                by=0, all=TRUE)
                                                
  bdiv.ordinations.sp[["points"]][[i]] = merge(scores(bdiv.ordinations.sp[["models"]][[i]],
                                                      display="sites", choices=c(1:2)),
                                               subset(sample_data, Type=="sponge")[,c("SampleID", "Source")],
                                               by=0, all=TRUE)
  
  # Calculate confidence ellipses for the all-types comparison.
  plot.new()
  temp <- ordiellipse(bdiv.ordinations.all[["models"]][[i]], 
                      bdiv.ordinations.all[["points"]][[i]]$Type, 
                      display="sites", 
                      kind="sd", 
                      conf=0.90, 
                      label=T)
  
  bdiv.ordinations.all[["ellipses"]][[i]] <- data.frame()
  
  for(g in levels(bdiv.ordinations.all[["points"]][[i]]$Type)){
    bdiv.ordinations.all[["ellipses"]][[i]] <- rbind(bdiv.ordinations.all[["ellipses"]][[i]],
                                            cbind(as.data.frame(with(bdiv.ordinations.all[["points"]][[i]][bdiv.ordinations.all[["points"]][[i]]$Type==g,],
                                                                     veganCovEllipse(temp[[g]]$cov,
                                                                                     temp[[g]]$center,
                                                                                     temp[[g]]$scale)))
                                                  ,Type=g))
  }
  
  # Calculate confidence ellipses for the only-sponges comparison.
  plot.new()
  temp <- ordiellipse(bdiv.ordinations.sp[["models"]][[i]], 
                      bdiv.ordinations.sp[["points"]][[i]]$Source, 
                      display="sites", 
                      kind="sd", 
                      conf=0.90, 
                      label=T)
  
  bdiv.ordinations.sp[["ellipses"]][[i]] <- data.frame()
  
  for(g in levels(bdiv.ordinations.sp[["points"]][[i]]$Source)){
    bdiv.ordinations.sp[["ellipses"]][[i]] <- rbind(bdiv.ordinations.sp[["ellipses"]][[i]],
                                               cbind(as.data.frame(with(bdiv.ordinations.sp[["points"]][[i]][bdiv.ordinations.sp[["points"]][[i]]$Source==g,],
                                                                        veganCovEllipse(temp[[g]]$cov,
                                                                                        temp[[g]]$center,
                                                                                        temp[[g]]$scale)))
                                                     ,Source=g))
  }
  
  rm(temp)
}

# Calculate PCA ordination for the Aitchison distance.
rownames(sample_data) <- sample_data$SampleID

for(i in c(3)){
  bdiv.ordinations.all[["models"]][[i]] <- rda(as.data.frame(otu_table(microbiome::transform(pseq.clr, "clr"))))
  bdiv.ordinations.sp[["models"]][[i]] <- rda(as.data.frame(otu_table(microbiome::transform(subset_samples(pseq.clr, Type=="sponge"), "clr"))))
  
  bdiv.ordinations.all[["points"]][[i]] <- merge(scores(bdiv.ordinations.all[["models"]][[i]], display="sites", choices=c(1:3)),
                                        sample_data[,c(1:3)],
                                        by=0, all=FALSE)
  bdiv.ordinations.sp[["points"]][[i]] <- merge(scores(bdiv.ordinations.sp[["models"]][[i]], display="sites", choices=c(1:3)),
                                           subset(sample_data, Type=="sponge")[,c(1:3)],
                                           by=0, all=FALSE)
  
  plot.new()
  temp <- ordiellipse(bdiv.ordinations.all[["models"]][[i]], 
                      bdiv.ordinations.all[["points"]][[i]]$Type, 
                      display="sites", 
                      kind="sd", 
                      conf=0.90, 
                      label=T)
  
  bdiv.ordinations.all[["ellipses"]][[i]] <- data.frame()
  
  for(g in levels(bdiv.ordinations.all[["points"]][[i]]$Type)){
    bdiv.ordinations.all[["ellipses"]][[i]] <- rbind(bdiv.ordinations.all[["ellipses"]][[i]],
                                            cbind(as.data.frame(with(bdiv.ordinations.all[["points"]][[i]][bdiv.ordinations.all[["points"]][[i]]$Type==g,],
                                                                     veganCovEllipse(temp[[g]]$cov,
                                                                                     temp[[g]]$center,
                                                                                     temp[[g]]$scale)))
                                                  ,Type=g))
  }
  
  # Calculate confidence ellipses for the only-sponges comparison.
  plot.new()
  temp <- ordiellipse(bdiv.ordinations.sp[["models"]][[i]], 
                      bdiv.ordinations.sp[["points"]][[i]]$Source, 
                      display="sites", 
                      kind="sd", 
                      conf=0.90, 
                      label=T)
  
  bdiv.ordinations.sp[["ellipses"]][[i]] <- data.frame()
  
  for(g in levels(bdiv.ordinations.sp[["points"]][[i]]$Source)){
    bdiv.ordinations.sp[["ellipses"]][[i]] <- rbind(bdiv.ordinations.sp[["ellipses"]][[i]],
                                               cbind(as.data.frame(with(bdiv.ordinations.sp[["points"]][[i]][bdiv.ordinations.sp[["points"]][[i]]$Source==g,],
                                                                        veganCovEllipse(temp[[g]]$cov,
                                                                                        temp[[g]]$center,
                                                                                        temp[[g]]$scale)))
                                                     ,Source=g))
  }
  
  rm(temp)
}

# Add location-based sub-ellipses to the bigger plots
bdiv.ordinations.all[["subellipses"]] <- list()

for(i in c(1:5)){
  bdiv.ordinations.all[["points"]][[i]]$Merge <- paste0(bdiv.ordinations.all[["points"]][[i]]$Source, "_",
                                                        bdiv.ordinations.all[["points"]][[i]]$Type)
  bdiv.ordinations.all[["points"]][[i]]$Merge[bdiv.ordinations.all[["points"]][[i]]$Merge %in% c("Sooke_biofilm", "Nanaimo_biofilm", "Cowichan_biofilm")] <- "biofilm"
  bdiv.ordinations.all[["points"]][[i]]$Merge <- factor(bdiv.ordinations.all[["points"]][[i]]$Merge)
  
  plot.new()
  temp <- ordiellipse(bdiv.ordinations.all[["models"]][[i]], 
                      bdiv.ordinations.all[["points"]][[i]]$Merge, 
                      display="sites", 
                      kind="sd", 
                      conf=0.90, 
                      label=T)
  
  bdiv.ordinations.all[["subellipses"]][[i]] <- data.frame()
  
  for(g in levels(bdiv.ordinations.all[["points"]][[i]]$Merge)){
    bdiv.ordinations.all[["subellipses"]][[i]] <- rbind(bdiv.ordinations.all[["subellipses"]][[i]],
                                                        cbind(as.data.frame(with(bdiv.ordinations.all[["points"]][[i]][bdiv.ordinations.all[["points"]][[i]]$Merge==g,],
                                                                                 veganCovEllipse(temp[[g]]$cov,
                                                                                                 temp[[g]]$center,
                                                                                                 temp[[g]]$scale)))
                                                              ,Merge=g))
  }
  
  bdiv.ordinations.all[["subellipses"]][[i]] <- subset(bdiv.ordinations.all[["subellipses"]][[i]], Merge !="biofilm")
  
  rm(temp)
}

# Calculate PERMANOVA results for each distance measure.
bdiv.permanova <- data.frame(matrix(ncol=6, nrow=10))
bdiv.permanova[,1] <- rep(c("Bray-Curtis", "Jaccard", "Aitchison", "wUF", "uwUF"), 2)
bdiv.permanova[,2] <- c(rep("all samples", 5), rep("sponges only", 5))
colnames(bdiv.permanova) <- c("Metric", "Comparison", "F", "df", "R2", "p")

for(i in c(1:10)){
  if(i %in% c(1:5)){ # By sample type
    temp <- adonis(dist.sp.list[[i]] ~ Type, sample_data, permutations=1000)
    bdiv.permanova[i,3] <- temp[["aov.tab"]]$F.Model[1]
    bdiv.permanova[i,4] <- temp[["aov.tab"]]$Df[1]
    bdiv.permanova[i,5] <- temp[["aov.tab"]]$R2[1]
    bdiv.permanova[i,6] <- temp[["aov.tab"]]$`Pr(>F)`[1]
  }
  if(i %in% c(6:10)){ # By sample source (sponges only)
    temp <- adonis(dist.sp.list.sponge[[i-5]] ~ Source, subset(sample_data, Type=="sponge"), permuations=1000)
    bdiv.permanova[i,3] <- temp[["aov.tab"]]$F.Model[1]
    bdiv.permanova[i,4] <- temp[["aov.tab"]]$Df[1]
    bdiv.permanova[i,5] <- temp[["aov.tab"]]$R2[1]
    bdiv.permanova[i,6] <- temp[["aov.tab"]]$`Pr(>F)`[1]
  }
  rm(temp)
}

bdiv.betadisper <- data.frame(matrix(ncol=5, nrow=10))
bdiv.betadisper[,1] <- rep(c("Bray-Curtis", "Jaccard", "Aitchison", "wUF", "uwUF"), 2)
bdiv.betadisper[,2] <- c(rep("all samples", 5), rep("sponges only", 5))
colnames(bdiv.betadisper) <- c("Metric", "Comparison", "F", "df", "p")

for(i in c(1:10)){
  if(i %in% c(1:5)){
    temp <- permutest(betadisper(dist.sp.list[[i]], sample_data$Type), permuations=1000, pairwise=TRUE)
    bdiv.betadisper[i,3] <- temp[["tab"]][1,4]
    bdiv.betadisper[i,4] <- temp[["tab"]][1,1]
    bdiv.betadisper[i,5] <- temp[["tab"]][1,6]
  }
  if(i %in% c(6:10)){
    temp <- permutest(betadisper(dist.sp.list.sponge[[i-5]], subset(sample_data, Type=="sponge")$Source), permuations=1000)
    bdiv.betadisper[i,3] <- temp[["tab"]][1,4]
    bdiv.betadisper[i,4] <- temp[["tab"]][1,1]
    bdiv.betadisper[i,5] <- temp[["tab"]][1,6]
  }
  rm(temp)
}

#### RANDOM FOREST MODELS ###########
rf.results <- list()
rf.results[["Models"]] <- list()
rf.results[["Gini_scores"]] <- list()

# This code is adapted from https://rpubs.com/michberr/randomforestmicrobe. 
# It is repeated twice: 
# first to classify samples based on type [sponge/water/biofllm], and then
# second to classify samples by location AND type
# third to classify only sponges by location

for(i in c(1:3)){
  prunescale = 0.0001
  
  if(i<3){temp.physeq <- pseq.rarefied}
  if(i==3){
    temp.physeq <- subset_samples(pseq.rarefied, Type=="sponge")
    temp.physeq <- prune_taxa(taxa_sums(temp.physeq) > 0, temp.physeq)
  }
  
  minlib <- min(sample_sums(temp.physeq))
  tax.mean <- taxa_sums(temp.physeq)/nsamples(temp.physeq)
  sites.prune <- prune_taxa(tax.mean > prunescale*minlib, temp.physeq)
  
  sites.prune <- subset_taxa(pseq.clr.sub, taxa_names(pseq.clr.sub) %in% taxa_names(sites.prune))
  if(i==3){sites.prune <- subset_samples(sites.prune, Type=="sponge")}
  
  sites.prune <- microbiome::transform(sites.prune, "clr")
  
  # Define the model predictors as ASV abundances.
  predictors <- otu_table(sites.prune)
  
  # Define the model response as the variable of interest (segment, individual identity, or intestinal site).
  if(i==1){response <- sample_data$Type}
  if(i==2){
    x <- as.data.frame(sample_data(sites.prune))
    x$Factor <- paste0(x$Source, x$Type)
    response <- x$Factor
    rm(x)
  }
  if(i==3){response <- sample_data(temp.physeq)$Source}
  
  # Combine the predictors and response into a single data frame.
  rf.data <- data.frame(response, predictors)
  
  # Perform random forest model.
  set.seed(2)
  rf.raw.results <- randomForest(response~., data = rf.data, ntree = 1000)
  
  # Extract the 40 most important ASVs, based on their mean decrease in the Gini coefficient.
  rf.importance <- randomForest::importance(rf.raw.results)
  rf.importance <- data.frame(predictors = rownames(rf.importance), rf.importance)
  rf.importance <- arrange(rf.importance, desc(MeanDecreaseGini))
  rf.importance$predictors <- factor(rf.importance$predictors,
                                     levels = rf.importance$predictors)
  rf.importance <- rf.importance[1:40, ]
  rownames(rf.importance) <- rf.importance$predictors
  rf.importance$predictors <- NULL
  
  # Add taxonomy information to the predictors.
  temp <- tax_table(sites.prune)
  temp <- subset(temp, rownames(temp) %in% rownames(rf.importance))
  rf.importance <- merge(rf.importance, temp, by="row.names", all=TRUE)
  colnames(rf.importance)[1] <- "Taxon"
  
  if(i==1){rf.importance <- merge(rf.importance, prev.data$ASV[,c("Taxon","Abundance","sponge","water","biofilm")],
                                  by="Taxon", all=FALSE)}
  if(i==2){rf.importance <- merge(rf.importance, prev.data$ASV[,c("Taxon","Abundance", "Sooke_sponge","Sooke_water","Sooke_biofilm",
                                                                  "Nanaimo_sponge", "Nanaimo_water", "Nanaimo_biofilm",
                                                                  "Cowichan_sponge", "Cowichan_water", "Cowichan_biofilm")],
                                  by="Taxon", all=FALSE)}
  if(i==3){rf.importance <- merge(rf.importance, prev.data$ASV[,c("Taxon","Abundance", "Sooke_sponge","Nanaimo_sponge","Cowichan_sponge")],
                                  by="Taxon", all=FALSE)}
  
  # Organize by Gini score.  
  rf.importance <- arrange(rf.importance, desc(MeanDecreaseGini))
  
  # Save results and clean workspace.
  rf.results[["Models"]][[i]] <- rf.raw.results
  rf.results[["Gini_scores"]][[i]] <- rf.importance
  rm(prunescale, minlib, tax.mean, sites.prune, predictors, response, rf.data, temp, rf.raw.results, rf.importance, temp.physeq)
}
names(rf.results[["Models"]]) <- c("type", "location_and_type", "sponges_by_location")
names(rf.results[["Gini_scores"]]) <- c("type", "location_and_type", "sponges_by_location")

#### NUMBER OF ASVs SHARED B/W EVERY PAIR OF SAMPLES [not included in final manuscript] ####
temp.physeq.all <- transform_sample_counts(pseq.clr, function(x) 100*x/sum(x)) # convert to relative abundance
temp.physeq.otu <- as.data.frame(cbind(otu_table(temp.physeq.all))) # extract OTU table
temp.physeq.otu <- as.data.frame(t(temp.physeq.otu))
temp.physeq.otu[temp.physeq.otu > 0] <- 1 # convert abundance to presence/absence

# Make a list of ASVs present in every sample
vector.list.of.unique.per.sample <- list()

for(i in c(1:ncol(temp.physeq.otu))){
  temp <- subset(temp.physeq.otu, temp.physeq.otu[,i] > 0)
  vector.list.of.unique.per.sample[[i]] <- rownames(temp)
  rm(temp)
}
names(vector.list.of.unique.per.sample) <- colnames(temp.physeq.otu)

# Find the number of ASVs shared between every pair of samples.
# This involves calculate the cumulative species richness of each pair of samples (df.shared.sp.richness)
# and the number of those ASVs common to both samples (vector.intersect.bw.samples)
vector.intersect.bw.samples <- data.frame(matrix(nrow=37, ncol=37))
df.shared.sp.richness <- c("temp", "fill", "fill")

for(i in c(1:length(vector.list.of.unique.per.sample))){
  for(j in c(1:37)){
    vector.intersect.bw.samples[i, j] <- length(intersect(vector.list.of.unique.per.sample[[i]],
                                                          vector.list.of.unique.per.sample[[j]]))
    
    temp <- temp.physeq.otu
    temp$TEMP <- rowSums(temp.physeq.otu[,c(i,j)])
    temp <- subset(temp, TEMP > 0)
    
    string <- c(names(vector.list.of.unique.per.sample)[i],
                names(vector.list.of.unique.per.sample)[j],
                nrow(temp))
    df.shared.sp.richness <- rbind(df.shared.sp.richness, string)
    rm(temp, string)
  }
}

# Clean the data frames.
df.shared.sp.richness <- as.data.frame(df.shared.sp.richness)
colnames(df.shared.sp.richness) <- c("Sample1", "Sample2", "Shared.Richness")
df.shared.sp.richness <-subset(df.shared.sp.richness, Sample2 !="fill")

colnames(vector.intersect.bw.samples) <- names(vector.list.of.unique.per.sample)
rownames(vector.intersect.bw.samples) <- names(vector.list.of.unique.per.sample)
vector.intersect.bw.samples$Sample2 <- rownames(vector.intersect.bw.samples)
vector.intersect.bw.samples <- reshape2::melt(vector.intersect.bw.samples)
colnames(vector.intersect.bw.samples) <- c("Sample1", "Sample2", "Intersect")

# Combine into a single data frame and ensure data is numeric.
vector.intersect.bw.samples <- merge(vector.intersect.bw.samples,
                                     df.shared.sp.richness,
                                     by=c("Sample1", "Sample2"),
                                     all=TRUE)

vector.intersect.bw.samples$Intersect  <- as.numeric(as.character(vector.intersect.bw.samples$Intersect))
vector.intersect.bw.samples$Shared.Richness  <- as.numeric(as.character(vector.intersect.bw.samples$Shared.Richness))

# Calculate percentage overlap between sample pairs as a function of the total number of ASVs.
vector.intersect.bw.samples$Percent <- 100*(vector.intersect.bw.samples$Intersect/
                                              vector.intersect.bw.samples$Shared.Richness)
vector.intersect.bw.samples <- subset(vector.intersect.bw.samples,
                                      Percent < 100) # Removes identical sample comparisons (e.g., COW1 to COW1)

# Add type and location information for each sample
merge_data <- sample_data[,c("SampleID", "Type", "Source")]
merge_data$Sample1 <- merge_data$SampleID
merge_data$Sample2 <- merge_data$SampleID
merge_data$Type1 <- merge_data$Type
merge_data$Type2 <- merge_data$Type
merge_data$Source1 <- merge_data$Source
merge_data$Source2 <- merge_data$Source

vector.intersect.bw.samples <- merge(vector.intersect.bw.samples,
                                     merge_data[,c("Sample1","Type1","Source1")],
                                     by="Sample1", all=TRUE)
vector.intersect.bw.samples <- merge(vector.intersect.bw.samples,
                                     merge_data[,c("Sample2","Type2","Source2")],
                                     by="Sample2", all=TRUE)

rm(merge_data, df.shared.sp.richness)

# Visualize sample overlap for all samples and just sponge samples.
ggplot(vector.intersect.bw.samples, aes(x=Sample1, y=Sample2, fill=Percent)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=90))

ggplot(subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge"),
       aes(x=Sample1, y=Sample2, fill=Percent)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_bw() + plot_theme + theme(axis.text.x=element_text(angle=90))

# Calculate the average percent of ASVs shared between sponge samples vs. the average percent shared between sponges and other samples.
all.samples1 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge") # sponge-sponge
all.samples2 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2 !="sponge") # sponge - not sponge
all.samples3 <- subset(vector.intersect.bw.samples, Type1 !="sponge" & Type2=="sponge") # sponge - not sponge
all.samples1$Group <- rep("sponge-sponge")
all.samples2$Group <- rep("sponge-other")
all.samples3$Group <- rep("sponge-other")
all.samples <- rbind(all.samples1[,c("Group","Percent")],
                     all.samples2[,c("Group","Percent")],
                     all.samples3[,c("Group","Percent")])
rm(all.samples1, all.samples2, all.samples3)
all.samples <- dplyr::distinct(all.samples)

t.test(Percent ~ Group, all.samples) # sponges share significantly more ASVs with each other than w/ other sample types
all.samples %>% dplyr::group_by(Group) %>% dplyr::summarise(mean=mean(Percent))

# Repeat this analysis considering only sponges from the same river.
sponge.samples1 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge" & Source1=="Sooke" & Source2=="Sooke")
sponge.samples2 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge" & Source1=="Sooke" & Source2!="Sooke")
sponge.samples3 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge" & Source1!="Sooke" & Source2=="Sooke")

sponge.samples4 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge" & Source1=="Nanaimo" & Source2=="Nanaimo")
sponge.samples5 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge" & Source1=="Nanaimo" & Source2!="Nanaimo")
sponge.samples6 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge" & Source1!="Nanaimo" & Source2=="Nanaimo")

sponge.samples7 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge" & Source1=="Cowichan" & Source2=="Cowichan")
sponge.samples8 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge" & Source1=="Cowichan" & Source2!="Cowichan")
sponge.samples9 <- subset(vector.intersect.bw.samples, Type1=="sponge" & Type2=="sponge" & Source1!="Cowichan" & Source2=="Cowichan")

sponge.samples1$Group <- rep("within-river")
sponge.samples2$Group <- rep("among-river")
sponge.samples3$Group <- rep("among-river")
sponge.samples4$Group <- rep("within-river")
sponge.samples5$Group <- rep("among-river")
sponge.samples6$Group <- rep("among-river")
sponge.samples7$Group <- rep("within-river")
sponge.samples8$Group <- rep("among-river")
sponge.samples9$Group <- rep("among-river")

sponge.samples <- rbind(sponge.samples1[,c("Group","Percent")],
                        sponge.samples2[,c("Group","Percent")],
                        sponge.samples3[,c("Group","Percent")],
                        sponge.samples4[,c("Group","Percent")],
                        sponge.samples5[,c("Group","Percent")],
                        sponge.samples6[,c("Group","Percent")],
                        sponge.samples7[,c("Group","Percent")],
                        sponge.samples8[,c("Group","Percent")],
                        sponge.samples9[,c("Group","Percent")])
rm(sponge.samples1, sponge.samples2, sponge.samples3, sponge.samples4, sponge.samples5,
   sponge.samples6, sponge.samples7, sponge.samples8, sponge.samples9)
sponge.samples <- dplyr::distinct(sponge.samples)

t.test(Percent ~ Group, sponge.samples) # sponges share significantly more ASVs within each river than with among rivers

rm(all.samples, sponge.samples, vector.list.of.unique.per.sample)

#### DIFFERENTIAL ABUNDANCE #####
taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")
library(ALDEx2)

# (a) By type, among all samples generally ####
diff.abund.TYPE <- list()
for(i in 1:6){
  temp.physeq <- pseq.clr.sub
  if(i < 6){
    temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])  
  }
  if(i > 5){
    temp.physeq <- subset_taxa(temp.physeq, Phylum != "NA")
  }
  temp.otu.table <- as.data.frame(otu_table(temp.physeq))
  temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
  
  if(i < 6){
    colnames(temp.otu.table) <- temp.tax[,i+1] }
  temp.otu.table <- as.data.frame(t(temp.otu.table))
  
  covariates <- as.character(sample_data(temp.physeq)$Type)
  p <- aldex(temp.otu.table, covariates, mc.samples=128, test="kw", effect=FALSE, denom="all")
  
  temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.transform.otu <- as.data.frame(otu_table(temp.transform))
  
  if(i < 6){
    colnames(temp.transform.otu) <- temp.tax[,i+1] }
  
  temp.transform.otu <- cbind(sample_data(temp.physeq)$Type, temp.transform.otu)
  colnames(temp.transform.otu)[1] <- "Type"
  
  vector.means <- data.frame(matrix(nrow=9))
  for(j in 2:(ncol(temp.transform.otu))){
    poo1 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), mean)
    poo1$Group.1 <- NULL
    rownames(poo1) <- levels(temp.transform.otu$Type)
    rownames(poo1) <- paste("mean", rownames(poo1))
    
    poo2 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), median)
    poo2$Group.1 <- NULL
    rownames(poo2) <- levels(temp.transform.otu$Type)
    rownames(poo2) <- paste("median", rownames(poo2))
    
    poo3 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), sd)
    poo3$Group.1 <- NULL
    rownames(poo3) <- levels(temp.transform.otu$Type)
    rownames(poo3) <- paste("stdev", rownames(poo3))
    
    poo <- rbind(poo1, poo2, poo3)
    colnames(poo) <- colnames(temp.transform.otu)[j]
    
    vector.means <- cbind(vector.means, poo)
  }
  
  colnames(vector.means)[1] <- "Type"
  vector.means$Type <- NULL
  vector.means <- as.data.frame(t(vector.means))
  vector.means <- cbind(temp.tax, vector.means)
  if(i < 6){
    rownames(vector.means) <- vector.means[,i+1]}
  final <- merge(p, vector.means, by=0, all=TRUE)
  rownames(final) <- final$Row.names
  final$Row.names <- NULL
  
  diff.abund.TYPE[[i]] <- final
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq,
     temp.tax, temp.otu.table, covariates)
}

names(diff.abund.TYPE) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# (b) Among sponges, by location ####
taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")
diff.abund.sp.location <- list()
for(i in 1:6){
  temp.physeq <- subset_samples(pseq.clr.sub, Type=="sponge")
  if(i < 6){
    temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])  
  }
  if(i > 5){
    temp.physeq <- subset_taxa(temp.physeq, Phylum != "NA")
  }
  temp.otu.table <- as.data.frame(otu_table(temp.physeq))
  temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
  
  if(i < 6){
    colnames(temp.otu.table) <- temp.tax[,i+1] }
  temp.otu.table <- as.data.frame(t(temp.otu.table))
  
  covariates <- as.character(sample_data(temp.physeq)$Source)
  p <- aldex(temp.otu.table, covariates, mc.samples=128, test="kw", effect=FALSE, denom="all")
  
  temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.transform.otu <- as.data.frame(otu_table(temp.transform))
  
  if(i < 6){
    colnames(temp.transform.otu) <- temp.tax[,i+1] }
  
  temp.transform.otu <- cbind(sample_data(temp.physeq)$Source, temp.transform.otu)
  colnames(temp.transform.otu)[1] <- "Source"
  
  vector.means <- data.frame(matrix(nrow=9))
  for(j in 2:(ncol(temp.transform.otu))){
    poo1 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Source), mean)
    poo1$Group.1 <- NULL
    rownames(poo1) <- levels(temp.transform.otu$Source)
    rownames(poo1) <- paste("mean", rownames(poo1))
    
    poo2 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Source), median)
    poo2$Group.1 <- NULL
    rownames(poo2) <- levels(temp.transform.otu$Source)
    rownames(poo2) <- paste("median", rownames(poo2))
    
    poo3 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Source), sd)
    poo3$Group.1 <- NULL
    rownames(poo3) <- levels(temp.transform.otu$Source)
    rownames(poo3) <- paste("stdev", rownames(poo3))
    
    poo <- rbind(poo1, poo2, poo3)
    colnames(poo) <- colnames(temp.transform.otu)[j]
    
    vector.means <- cbind(vector.means, poo)
  }
  
  colnames(vector.means)[1] <- "Source"
  vector.means$Source <- NULL
  vector.means <- as.data.frame(t(vector.means))
  vector.means <- cbind(temp.tax, vector.means)
  if(i < 6){
    rownames(vector.means) <- vector.means[,i+1]}
  final <- merge(p, vector.means, by=0, all=TRUE)
  rownames(final) <- final$Row.names
  final$Row.names <- NULL
  
  diff.abund.sp.location[[i]] <- final
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq,
     temp.tax, temp.otu.table, covariates)
}

names(diff.abund.sp.location) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# (c) Among water, by location ####
taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")
diff.abund.wat.location <- list()
for(i in 1:6){
  temp.physeq <- subset_samples(pseq.clr.sub, Type=="water")
  if(i < 6){
    temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])  
  }
  if(i > 5){
    temp.physeq <- subset_taxa(temp.physeq, Phylum != "NA")
  }
  temp.otu.table <- as.data.frame(otu_table(temp.physeq))
  temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
  
  if(i < 6){
    colnames(temp.otu.table) <- temp.tax[,i+1] }
  temp.otu.table <- as.data.frame(t(temp.otu.table))
  
  covariates <- as.character(sample_data(temp.physeq)$Source)
  p <- aldex(temp.otu.table, covariates, mc.samples=128, test="kw", effect=FALSE, denom="all")
  
  temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.transform.otu <- as.data.frame(otu_table(temp.transform))
  
  if(i < 6){
    colnames(temp.transform.otu) <- temp.tax[,i+1] }
  
  temp.transform.otu <- cbind(sample_data(temp.physeq)$Source, temp.transform.otu)
  colnames(temp.transform.otu)[1] <- "Source"
  
  vector.means <- data.frame(matrix(nrow=9))
  for(j in 2:(ncol(temp.transform.otu))){
    poo1 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Source), mean)
    poo1$Group.1 <- NULL
    rownames(poo1) <- levels(temp.transform.otu$Source)
    rownames(poo1) <- paste("mean", rownames(poo1))
    
    poo2 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Source), median)
    poo2$Group.1 <- NULL
    rownames(poo2) <- levels(temp.transform.otu$Source)
    rownames(poo2) <- paste("median", rownames(poo2))
    
    poo3 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Source), sd)
    poo3$Group.1 <- NULL
    rownames(poo3) <- levels(temp.transform.otu$Source)
    rownames(poo3) <- paste("stdev", rownames(poo3))
    
    poo <- rbind(poo1, poo2, poo3)
    colnames(poo) <- colnames(temp.transform.otu)[j]
    
    vector.means <- cbind(vector.means, poo)
  }
  
  colnames(vector.means)[1] <- "Source"
  vector.means$Source <- NULL
  vector.means <- as.data.frame(t(vector.means))
  vector.means <- cbind(temp.tax, vector.means)
  if(i < 6){
    rownames(vector.means) <- vector.means[,i+1]}
  final <- merge(p, vector.means, by=0, all=TRUE)
  rownames(final) <- final$Row.names
  final$Row.names <- NULL
  
  diff.abund.wat.location[[i]] <- final
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq,
     temp.tax, temp.otu.table, covariates)
}

names(diff.abund.wat.location) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# (d) Pairwise: sponge-water ####
diff.abund.SPWAT <- list()
for(i in 1:6){
  temp.physeq <- subset_samples(pseq.clr.sub, Type !="biofilm")
  if(i < 6){
    temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])  
  }
  if(i > 5){
    temp.physeq <- subset_taxa(temp.physeq, Phylum != "NA")
  }
  temp.otu.table <- as.data.frame(otu_table(temp.physeq))
  temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
  
  if(i < 6){
    colnames(temp.otu.table) <- temp.tax[,i+1] }
  temp.otu.table <- as.data.frame(t(temp.otu.table))
  
  covariates <- as.character(sample_data(temp.physeq)$Type)
  p <- aldex(temp.otu.table, covariates, mc.samples=128, test="t", effect=FALSE, denom="all")
  
  temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.transform.otu <- as.data.frame(otu_table(temp.transform))
  
  if(i < 6){
    colnames(temp.transform.otu) <- temp.tax[,i+1] }
  
  temp.transform.otu <- cbind(sample_data(temp.physeq)$Type, temp.transform.otu)
  colnames(temp.transform.otu)[1] <- "Type"
  
  vector.means <- data.frame(matrix(nrow=6))
  for(j in 2:(ncol(temp.transform.otu))){
    poo1 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), mean)
    poo1$Group.1 <- NULL
    rownames(poo1) <- levels(temp.transform.otu$Type)
    rownames(poo1) <- paste("mean", rownames(poo1))
    
    poo2 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), median)
    poo2$Group.1 <- NULL
    rownames(poo2) <- levels(temp.transform.otu$Type)
    rownames(poo2) <- paste("median", rownames(poo2))
    
    poo3 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), sd)
    poo3$Group.1 <- NULL
    rownames(poo3) <- levels(temp.transform.otu$Type)
    rownames(poo3) <- paste("stdev", rownames(poo3))
    
    poo <- rbind(poo1, poo2, poo3)
    colnames(poo) <- colnames(temp.transform.otu)[j]
    
    vector.means <- cbind(vector.means, poo)
  }
  
  colnames(vector.means)[1] <- "Type"
  vector.means$Type <- NULL
  vector.means <- as.data.frame(t(vector.means))
  vector.means <- cbind(temp.tax, vector.means)
  if(i < 6){
    rownames(vector.means) <- vector.means[,i+1]}
  final <- merge(p, vector.means, by=0, all=TRUE)
  rownames(final) <- final$Row.names
  final$Row.names <- NULL
  
  diff.abund.SPWAT[[i]] <- final
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq,
     temp.tax, temp.otu.table, covariates)
}

names(diff.abund.SPWAT) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# (e) Pairwise: sponge-biofilm #####
diff.abund.SPBIO <- list()
for(i in 1:6){
  temp.physeq <- subset_samples(pseq.clr.sub, Type !="water")
  if(i < 6){
    temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])  
  }
  if(i > 5){
    temp.physeq <- subset_taxa(temp.physeq, Phylum != "NA")
  }
  temp.otu.table <- as.data.frame(otu_table(temp.physeq))
  temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
  
  if(i < 6){
    colnames(temp.otu.table) <- temp.tax[,i+1] }
  temp.otu.table <- as.data.frame(t(temp.otu.table))
  
  covariates <- as.character(sample_data(temp.physeq)$Type)
  p <- aldex(temp.otu.table, covariates, mc.samples=128, test="t", effect=FALSE, denom="all")
  
  temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.transform.otu <- as.data.frame(otu_table(temp.transform))
  
  if(i < 6){
    colnames(temp.transform.otu) <- temp.tax[,i+1] }
  
  temp.transform.otu <- cbind(sample_data(temp.physeq)$Type, temp.transform.otu)
  colnames(temp.transform.otu)[1] <- "Type"
  
  vector.means <- data.frame(matrix(nrow=6))
  for(j in 2:(ncol(temp.transform.otu))){
    poo1 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), mean)
    poo1$Group.1 <- NULL
    rownames(poo1) <- levels(temp.transform.otu$Type)
    rownames(poo1) <- paste("mean", rownames(poo1))
    
    poo2 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), median)
    poo2$Group.1 <- NULL
    rownames(poo2) <- levels(temp.transform.otu$Type)
    rownames(poo2) <- paste("median", rownames(poo2))
    
    poo3 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), sd)
    poo3$Group.1 <- NULL
    rownames(poo3) <- levels(temp.transform.otu$Type)
    rownames(poo3) <- paste("stdev", rownames(poo3))
    
    poo <- rbind(poo1, poo2, poo3)
    colnames(poo) <- colnames(temp.transform.otu)[j]
    
    vector.means <- cbind(vector.means, poo)
  }
  
  colnames(vector.means)[1] <- "Type"
  vector.means$Type <- NULL
  vector.means <- as.data.frame(t(vector.means))
  vector.means <- cbind(temp.tax, vector.means)
  if(i < 6){
    rownames(vector.means) <- vector.means[,i+1]}
  final <- merge(p, vector.means, by=0, all=TRUE)
  rownames(final) <- final$Row.names
  final$Row.names <- NULL
  
  diff.abund.SPBIO[[i]] <- final
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq,
     temp.tax, temp.otu.table, covariates)
}

names(diff.abund.SPBIO) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# (f) Pairwise: water-biofilm ####
diff.abund.WATBIO <- list()
for(i in 1:6){
  temp.physeq <- subset_samples(pseq.clr.sub, Type !="sponge")
  if(i < 6){
    temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])  
  }
  if(i > 5){
    temp.physeq <- subset_taxa(temp.physeq, Phylum != "NA")
  }
  temp.otu.table <- as.data.frame(otu_table(temp.physeq))
  temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
  
  if(i < 6){
    colnames(temp.otu.table) <- temp.tax[,i+1] }
  temp.otu.table <- as.data.frame(t(temp.otu.table))
  
  covariates <- as.character(sample_data(temp.physeq)$Type)
  p <- aldex(temp.otu.table, covariates, mc.samples=128, test="t", effect=FALSE, denom="all")
  
  temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.transform.otu <- as.data.frame(otu_table(temp.transform))
  
  if(i < 6){
    colnames(temp.transform.otu) <- temp.tax[,i+1] }
  
  temp.transform.otu <- cbind(sample_data(temp.physeq)$Type, temp.transform.otu)
  colnames(temp.transform.otu)[1] <- "Type"
  
  vector.means <- data.frame(matrix(nrow=6))
  for(j in 2:(ncol(temp.transform.otu))){
    poo1 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), mean)
    poo1$Group.1 <- NULL
    rownames(poo1) <- levels(temp.transform.otu$Type)
    rownames(poo1) <- paste("mean", rownames(poo1))
    
    poo2 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), median)
    poo2$Group.1 <- NULL
    rownames(poo2) <- levels(temp.transform.otu$Type)
    rownames(poo2) <- paste("median", rownames(poo2))
    
    poo3 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Type), sd)
    poo3$Group.1 <- NULL
    rownames(poo3) <- levels(temp.transform.otu$Type)
    rownames(poo3) <- paste("stdev", rownames(poo3))
    
    poo <- rbind(poo1, poo2, poo3)
    colnames(poo) <- colnames(temp.transform.otu)[j]
    
    vector.means <- cbind(vector.means, poo)
  }
  
  colnames(vector.means)[1] <- "Type"
  vector.means$Type <- NULL
  vector.means <- as.data.frame(t(vector.means))
  vector.means <- cbind(temp.tax, vector.means)
  if(i < 6){
    rownames(vector.means) <- vector.means[,i+1]}
  final <- merge(p, vector.means, by=0, all=TRUE)
  rownames(final) <- final$Row.names
  final$Row.names <- NULL
  
  diff.abund.WATBIO[[i]] <- final
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq,
     temp.tax, temp.otu.table, covariates)
}

names(diff.abund.WATBIO) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# Put pairwise results together ####
diff.abund.PAIRWISE <- list()
for(i in 1:6){
    df1 <- diff.abund.SPWAT[[i]]
    colnames(df1) <- c("SPWAT.we.ep", "SPWAT.we.eBH", "SPWAT.wi.ep", "SPWAT.wi.eBH", colnames(df1)[c(5:ncol(df1))])
    df2 <- diff.abund.SPBIO[[i]][,c(!colnames(diff.abund.SPBIO[[i]]) %in% c("Kingdom", "Phylum", "Class", "Order",
                                                                            "Family", "Genus", "mean sponge", "median sponge",
                                                                            "stdev sponge"))]
    colnames(df2) <- c("SPBIO.we.ep", "SPBIO.we.eBH", "SPBIO.wi.ep", "SPBIO.wi.eBH", colnames(df2)[c(5:ncol(df2))])
    df3 <- diff.abund.WATBIO[[i]][,c(1:4)]
    colnames(df3) <- c("WATBIO.we.ep", "WATBIO.we.eBH", "WATBIO.wi.ep", "WATBIO.wi.eBH")
    df <- merge(df1, df2, by=0, all=TRUE)
    rownames(df) <- df$Row.names
    df$Row.names <- NULL
    df <- merge(df, df3, by=0, all=TRUE)
    rownames(df) <- df$Row.names
    df$Row.names <- NULL
    
    df <- df[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "mean sponge", "stdev sponge", "mean water",
                "stdev water", "mean biofilm", "stdev biofilm", "SPWAT.we.eBH", "SPBIO.we.eBH", "WATBIO.we.eBH",
                "SPWAT.wi.eBH", "SPBIO.wi.eBH", "WATBIO.wi.eBH", "SPWAT.we.ep", "SPBIO.we.ep", "WATBIO.we.ep",
                "SPWAT.wi.ep", "SPBIO.wi.ep", "WATBIO.wi.ep")]
    diff.abund.PAIRWISE[[i]] <- df
    rm(df, df1, df2, df3)
}
names(diff.abund.PAIRWISE) <- c(taxranks, "ASV")
rm(diff.abund.SPBIO, diff.abund.SPWAT, diff.abund.WATBIO)

#### SPONGE CORE TAXA ####
## Core microbiome: ASVs present in at least 12 sponge samples
# Convert to relative abundances
temp.physeq.all <- transform_sample_counts(pseq.clr.sub, function(x) 100*x/sum(x))

# Extract OTU table and taxonomy table
temp.physeq.otu <- as.data.frame(cbind(otu_table(temp.physeq.all)))
temp.physeq.otu <- as.data.frame(t(temp.physeq.otu))
temp.physeq.tax <- as.data.frame(cbind())

# Count total number of sponge samples in which each ASV appears
temp.physeq.otu$ALL <- rowSums(temp.physeq.otu[,c("SKE1", "SKE2", "SKE3", "SKE4", "SKE5",
                                                  "NMO1", "NMO2", "NMO3", "NMO4", "NMO5",
                                                  "COW1", "COW2", "COW3", "COW4", "COW5")] != 0)

# Count number of Sooke, Nanaimo, and Cowichan sponges in which each ASV appears
temp.physeq.otu$SKE <- rowSums(temp.physeq.otu[,c("SKE1", "SKE2", "SKE3", "SKE4", "SKE5")] != 0)
temp.physeq.otu$NMO <- rowSums(temp.physeq.otu[,c("NMO1", "NMO2", "NMO3", "NMO4", "NMO5")] != 0)
temp.physeq.otu$COW <- rowSums(temp.physeq.otu[,c("COW1", "COW2", "COW3", "COW4", "COW5")] != 0)

# Calculate the mean relative abundance of each ASV across all sponges
temp.physeq.otu$MEAN <- rowMeans(temp.physeq.otu[,c("SKE1", "SKE2", "SKE3", "SKE4", "SKE5",
                                                    "NMO1", "NMO2", "NMO3", "NMO4", "NMO5",
                                                    "COW1", "COW2", "COW3", "COW4", "COW5")])

# Calculate the number of water and biofilm samples in which each ASV appears
temp.physeq.otu$WATER <- rowSums(temp.physeq.otu[,c("SKEW1", "SKEW2", "SKEW3", "SKEW4", "SKEW5",
                                                    "NMOW1", "NMOW2", "NMOW3", "NMOW4", "NMOW5",
                                                    "COW1", "COW2", "COW4", "COW5")] != 0)
temp.physeq.otu$BIOFILM <- rowSums(temp.physeq.otu[,c("SKEBio1", "SKEBio23", "SKEBio4", "SKEBio5",
                                                      "NMOBio1", "NMOBio2", "NMOBio3",
                                                      "COWBio1")] != 0)

# Remove ASVs not present at all in sponges
temp.physeq.otu <- subset(temp.physeq.otu, ALL > 0)

# Find the size of the core microbiome
nrow(subset(temp.physeq.otu, ALL==15)) # 22 ASVs present in every sponge sample
nrow(subset(temp.physeq.otu, ALL >= 12)) # 92 ASVs present in at least 80% of samples

# Identify the taxa in the core microbiome
temp.physeq.otu <- merge(as.data.frame(cbind(tax_table(temp.physeq.all))),
                         temp.physeq.otu,
                         by=0, all.x=FALSE, all.y=TRUE)
core.taxa <- subset(temp.physeq.otu, ALL>=12)

rownames(core.taxa) <- core.taxa$Row.names
core.taxa$Row.names <- NULL

# Sort by relative abundance
core.taxa <- core.taxa[order(-core.taxa$MEAN), ]

# Save the core sponge microbiome as a separate table
core.table <- core.taxa[,c("Phylum","Family","Genus","ALL","SKE","NMO","COW","MEAN","WATER","BIOFILM")]
write.csv(core.table, "sponge.core.microbiome.csv")

#### Explore differentially abundant taxa in relation to the core microbiome
# See how many ASVs are more abundant in sponges
nrow(subset(diff.abund.PAIRWISE[[6]], SPWAT.wi.eBH < 0.05 & SPBIO.wi.eBH < 0.05)) # 35 taxa
nrow(subset(diff.abund.PAIRWISE[[6]], SPWAT.we.eBH < 0.05 & SPBIO.we.eBH < 0.05)) # 19 taxa

temp <- subset(diff.abund.TYPE[[6]], glm.eBH < 0.05) # 134 differentially abundant ASVs
temp$spwat <- temp$`mean sponge` / temp$`mean water`
temp$spbio <- temp$`mean sponge` / temp$`mean biofilm`
temp <- subset(temp, spwat > 1 & spbio > 1) # 39 are more abundant in sponges

temp <- subset(temp, rownames(temp) %in% rownames(core.taxa)) # 34 differentially abundant ASVs in the core microbiome

temp <- subset(diff.abund.TYPE[[5]], glm.eBH < 0.05) # 80 differentially abundant genera
temp$spwat <- temp$`mean sponge` / temp$`mean water`
temp$spbio <- temp$`mean sponge` / temp$`mean biofilm`
temp <- subset(temp, spwat > 1 & spbio > 1) # 10 are more abundant in sponges

temp <- subset(diff.abund.TYPE[[4]], glm.eBH < 0.05) # 60 differentially abundant families
temp$spwat <- temp$`mean sponge` / temp$`mean water`
temp$spbio <- temp$`mean sponge` / temp$`mean biofilm`
temp <- subset(temp, spwat > 1 & spbio > 1) # 7 are more abundant in sponges

rm(temp.physeq.all, temp.physeq.otu, temp.physeq.tax, temp)

#### ASV OVERLAP ANALYSIS #########
# Create separate phyloseq objects for (1) sponge, (2) water, (3) biofilm, (4) SKE sponge, (5) NMO sponge, (6) COW sponge
temp.physeq <- pseq.clr

temp.physeq.sponge <- subset_samples(temp.physeq, Type=="sponge")
temp.physeq.sponge <- prune_taxa(taxa_sums(temp.physeq.sponge) > 0, temp.physeq.sponge)

temp.physeq.water <- subset_samples(temp.physeq, Type=="water")
temp.physeq.water <- prune_taxa(taxa_sums(temp.physeq.water) > 0, temp.physeq.water)

temp.physeq.biofilm <- subset_samples(temp.physeq, Type=="biofilm")
temp.physeq.biofilm <- prune_taxa(taxa_sums(temp.physeq.biofilm) > 0, temp.physeq.biofilm)

temp.physeq.sponge.SKE <- subset_samples(temp.physeq.sponge, Source=="Sooke")
temp.physeq.sponge.SKE <- prune_taxa(taxa_sums(temp.physeq.sponge.SKE) > 0, temp.physeq.sponge.SKE)

temp.physeq.sponge.NMO <- subset_samples(temp.physeq.sponge, Source=="Nanaimo")
temp.physeq.sponge.NMO <- prune_taxa(taxa_sums(temp.physeq.sponge.NMO) > 0, temp.physeq.sponge.NMO)

temp.physeq.sponge.COW <- subset_samples(temp.physeq.sponge, Source=="Cowichan")
temp.physeq.sponge.COW <- prune_taxa(taxa_sums(temp.physeq.sponge.COW) > 0, temp.physeq.sponge.COW)

temp.physeq.SKE <- subset_samples(temp.physeq, Source=="Sooke")
temp.physeq.SKE <- prune_taxa(taxa_sums(temp.physeq.SKE) > 0, temp.physeq.SKE)

temp.physeq.NMO <- subset_samples(temp.physeq, Source=="Nanaimo")
temp.physeq.NMO <- prune_taxa(taxa_sums(temp.physeq.NMO) > 0, temp.physeq.NMO)

temp.physeq.COW <- subset_samples(temp.physeq, Source=="Cowichan")
temp.physeq.COW <- prune_taxa(taxa_sums(temp.physeq.COW) > 0, temp.physeq.COW)

# Calculate PERMANOVA results (Aitchison distance) for each comparison
# (among sample types within a river, and among rivers for each sample type)
temp.list <- list(temp.physeq.sponge, temp.physeq.water, temp.physeq.biofilm,
                  temp.physeq.SKE, temp.physeq.NMO, temp.physeq.COW)
temp.dist <- list()
bdiv.pairwise.permanova <- data.frame(matrix(ncol=7))
  
for(i in 1:length(temp.list)){
  temp.list[[i]] <- microbiome::transform(temp.list[[i]], "clr")
  temp.dist[[i]] <- as.data.frame(otu_table(temp.list[[i]]))
  temp.dist[[i]] <- vegdist(temp.dist[[i]], method="euclidean")
  if(i %in% c(1:3)){
    temp <- adonis(temp.dist[[i]] ~ sample_data(temp.list[[i]])$Source)
    tempd <- permutest(betadisper(temp.dist[[i]], sample_data(temp.list[[i]])$Source))
  }
  if(i %in% c(4:6)){
    temp <- adonis(temp.dist[[i]] ~ sample_data(temp.list[[i]])$Type)
    tempd <- permutest(betadisper(temp.dist[[i]], sample_data(temp.list[[i]])$Type))
  }
  bdiv.pairwise.permanova[i,2] <- temp[["aov.tab"]]$F.Model[1]
  bdiv.pairwise.permanova[i,3] <- temp[["aov.tab"]]$Df[1]
  bdiv.pairwise.permanova[i,4] <- temp[["aov.tab"]]$R2[1]
  bdiv.pairwise.permanova[i,5] <- temp[["aov.tab"]]$`Pr(>F)`[1]
  bdiv.pairwise.permanova[i,6] <- tempd$tab[1,4]
  bdiv.pairwise.permanova[i,7] <- tempd$tab[1,6]
  rm(temp, tempd)
}
bdiv.pairwise.permanova[,1] <- c("sponge", "water", "biofilm", "SKE", "NMO", "COW")
colnames(bdiv.pairwise.permanova) <- c("ordination", "F", "Df", "R2", "p", "Disp_F", "Disp_p")

rm(temp.list, temp.dist, temp.perm)

# Create character strings with the names of taxa in each object.
sponge <- taxa_names(temp.physeq.sponge)
water <- taxa_names(temp.physeq.water)
biofilm <- taxa_names(temp.physeq.biofilm)
SKE_sponge <- taxa_names(temp.physeq.sponge.SKE)
NMO_sponge <- taxa_names(temp.physeq.sponge.NMO)
COW_sponge <- taxa_names(temp.physeq.sponge.COW)

# Calculate numbers for the Venn diagram (overlap in eahc section) 
length(sponge) # Total ASVs in sponges = 1959 / 812
length(water) # Total ASVs in water = 2073 / 783
length(biofilm) # Total ASVs in biofilm = 2042 / 824

length(setdiff(sponge, c(water, biofilm))) # Unique to sponges = 134 / 14
length(setdiff(water, c(sponge, biofilm))) # Unique to water = 107 / 1
length(setdiff(biofilm, c(sponge, water))) # Unique to biofilm = 182 / 45

length(intersect(sponge, water)) # Sponge + water = 1669 / 739
length(intersect(sponge, biofilm)) # Sponge + biofilm = 1563 / 736
length(intersect(water, biofilm)) # Water + biofilm = 1704 / 720

length(intersect(sponge, intersect(water, biofilm))) # Shared by all three = 1407 / 677

length(SKE_sponge) # Total ASVs in SKE_sponges = 1246 / 615
length(NMO_sponge) # Total ASVs in NMO_sponge = 1143 / 565
length(COW_sponge) # Total ASVs in COW_sponge = 1238 / 615

length(setdiff(SKE_sponge, c(NMO_sponge, COW_sponge))) # Unique to SKE_sponges = 276 / 63
length(setdiff(NMO_sponge, c(SKE_sponge, COW_sponge))) # Unique to NMO_sponge = 232 / 43
length(setdiff(COW_sponge, c(SKE_sponge, NMO_sponge))) # Unique to COW_sponge = 315 / 84

length(intersect(SKE_sponge, NMO_sponge)) # SKE_sponge + NMO_sponge = 745 / 452
length(intersect(SKE_sponge, COW_sponge)) # SKE_sponge + COW_sponge = 757 / 461
length(intersect(NMO_sponge, COW_sponge)) # NMO_sponge + COW_sponge = 698 / 431

length(intersect(SKE_sponge, intersect(NMO_sponge, COW_sponge))) # Shared by all three = 532 / 361

# Calculate the prevalence and mean relative abundance of each unique taxon.
Venn.histograms <- list()

temp.physeq.sponge <- transform_sample_counts(temp.physeq.sponge, function(x) 100*x/sum(x))
temp.physeq.sponge <- subset_taxa(temp.physeq.sponge,
                                  taxa_names(temp.physeq.sponge) %in% setdiff(sponge, c(water, biofilm)))
Venn.histograms[["sponge"]] <- as.data.frame(microbiome::prevalence(temp.physeq.sponge, detection=0, count=TRUE))
temp.otu <- as.data.frame(otu_table(temp.physeq.sponge))
Venn.histograms[["sponge"]]$percent <- Venn.histograms[["sponge"]][,1]/15
Venn.histograms[["sponge"]]$mean <- colMeans(temp.otu)
colnames(Venn.histograms[["sponge"]]) <- c("prevalence", "percent", "mean")

temp.physeq.water <- transform_sample_counts(temp.physeq.water, function(x) 100*x/sum(x))
temp.physeq.water <- subset_taxa(temp.physeq.water,
                                  taxa_names(temp.physeq.water) %in% setdiff(water, c(sponge, biofilm)))
Venn.histograms[["water"]] <- as.data.frame(microbiome::prevalence(temp.physeq.water, detection=0, count=TRUE))
temp.otu <- as.data.frame(otu_table(temp.physeq.water))
Venn.histograms[["water"]]$percent <- Venn.histograms[["water"]][,1]/14
Venn.histograms[["water"]]$mean <- colMeans(temp.otu)
colnames(Venn.histograms[["water"]]) <- c("prevalence", "percent", "mean")

temp.physeq.biofilm <- transform_sample_counts(temp.physeq.biofilm, function(x) 100*x/sum(x))
temp.physeq.biofilm <- subset_taxa(temp.physeq.biofilm,
                                 taxa_names(temp.physeq.biofilm) %in% setdiff(biofilm, c(sponge, water)))
Venn.histograms[["biofilm"]] <- as.data.frame(microbiome::prevalence(temp.physeq.biofilm, detection=0, count=TRUE))
temp.otu <- as.data.frame(otu_table(temp.physeq.biofilm))
Venn.histograms[["biofilm"]]$percent <- Venn.histograms[["biofilm"]][,1]/8
Venn.histograms[["biofilm"]]$mean <- colMeans(temp.otu)
colnames(Venn.histograms[["biofilm"]]) <- c("prevalence", "percent", "mean")

temp.physeq.sponge.SKE <- transform_sample_counts(temp.physeq.sponge.SKE, function(x) 100*x/sum(x))
temp.physeq.sponge.SKE <- subset_taxa(temp.physeq.sponge.SKE,
                                   taxa_names(temp.physeq.sponge.SKE) %in% setdiff(SKE_sponge, c(NMO_sponge, COW_sponge)))
Venn.histograms[["SKE_sponge"]] <- as.data.frame(microbiome::prevalence(temp.physeq.sponge.SKE, detection=0, count=TRUE))
temp.otu <- as.data.frame(otu_table(temp.physeq.sponge.SKE))
Venn.histograms[["SKE_sponge"]]$percent <- Venn.histograms[["SKE_sponge"]][,1]/5
Venn.histograms[["SKE_sponge"]]$mean <- colMeans(temp.otu)
colnames(Venn.histograms[["SKE_sponge"]]) <- c("prevalence", "percent", "mean")

temp.physeq.sponge.NMO <- transform_sample_counts(temp.physeq.sponge.NMO, function(x) 100*x/sum(x))
temp.physeq.sponge.NMO <- subset_taxa(temp.physeq.sponge.NMO,
                                      taxa_names(temp.physeq.sponge.NMO) %in% setdiff(NMO_sponge, c(SKE_sponge, COW_sponge)))
Venn.histograms[["NMO_sponge"]] <- as.data.frame(microbiome::prevalence(temp.physeq.sponge.NMO, detection=0, count=TRUE))
temp.otu <- as.data.frame(otu_table(temp.physeq.sponge.NMO))
Venn.histograms[["NMO_sponge"]]$percent <- Venn.histograms[["NMO_sponge"]][,1]/5
Venn.histograms[["NMO_sponge"]]$mean <- colMeans(temp.otu)
colnames(Venn.histograms[["NMO_sponge"]]) <- c("prevalence", "percent", "mean")

temp.physeq.sponge.COW <- transform_sample_counts(temp.physeq.sponge.COW, function(x) 100*x/sum(x))
temp.physeq.sponge.COW <- subset_taxa(temp.physeq.sponge.COW,
                                      taxa_names(temp.physeq.sponge.COW) %in% setdiff(COW_sponge, c(SKE_sponge, NMO_sponge)))
Venn.histograms[["COW_sponge"]] <- as.data.frame(microbiome::prevalence(temp.physeq.sponge.COW, detection=0, count=TRUE))
temp.otu <- as.data.frame(otu_table(temp.physeq.sponge.COW))
Venn.histograms[["COW_sponge"]]$percent <- Venn.histograms[["COW_sponge"]][,1]/5
Venn.histograms[["COW_sponge"]]$mean <- colMeans(temp.otu)
colnames(Venn.histograms[["COW_sponge"]]) <- c("prevalence", "percent", "mean")

# River-specific ASVs?
nrow(subset(Venn.histograms$SKE_sponge, prevalence >=4 & mean >= 0.01))
nrow(subset(Venn.histograms$NMO_sponge, prevalence >=4 & mean >= 0.01))
nrow(subset(Venn.histograms$COW_sponge, prevalence >=4 & mean >= 0.01))

r1 <- subset(Venn.histograms$SKE_sponge, prevalence >=4 & mean >= 0.01)
r1$river <- rep("Sooke")
r2 <- subset(Venn.histograms$NMO_sponge, prevalence >=4 & mean >= 0.01)
r2$river <- rep("Nanaimo")
r3 <- subset(Venn.histograms$COW_sponge, prevalence >=4 & mean >= 0.01)
r3$river <- rep("Cowichan")

Venn.river.specific <- rbind(r1, r2, r3)

Venn.river.specific <- merge(Venn.river.specific,
                        as.data.frame(cbind(tax_table(pseq.clr.sub))),
                        by=0, all=FALSE)

Venn.river.specific$water <- (Venn.river.specific$Row.names %in% water)
Venn.river.specific$biofilm <- (Venn.river.specific$Row.names %in% biofilm)
  
rm(sponge, water, biofilm, SKE_sponge, NMO_sponge, COW_sponge,
   temp.physeq.sponge, temp.physeq.water, temp.physeq.biofilm,
   temp.physeq.sponge.SKE, temp.physeq.sponge.NMO, temp.physeq.sponge.COW,
   temp.physeq.COW, temp.physeq.NMO, temp.physeq.SKE,
   temp.physeq, temp.otu, r1, r2, r3)

#### ASV ABUNDANCE CORRELATIONS ######
# Calculate average abundance of each ASV in each sample type.
temp.physeq <- pseq.clr.sub
temp.physeq <- subset_taxa(temp.physeq, Phylum !="NA")
temp.physeq.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))

temp.otu <- as.data.frame(t(otu_table(temp.physeq.transform)))
temp.tax <- as.data.frame(cbind(tax_table(temp.physeq.transform)))
rownames(temp.otu) <- rownames(temp.tax)
temp.otu <- as.data.frame(t(temp.otu))
temp.otu$Type <- sample_data(temp.physeq.transform)$Type

correlation.summary <- temp.otu %>% dplyr::group_by(Type) %>% dplyr::summarise_if(is.numeric, mean)
rownames(correlation.summary) <- correlation.summary$Type
correlation.summary$Type <- NULL
correlation.summary <- as.data.frame(t(correlation.summary))

for(i in 1:ncol(correlation.summary)){
  correlation.summary[,c(i)] <- as.numeric(as.character(correlation.summary[,c(i)]))
}

correlation.summary <- cbind(temp.tax, correlation.summary)
rm(temp.physeq.transform, temp.otu, temp.tax)

# Calculate pairwise differential abundance at the ASV level (sponge-water and sponge-biofilm)
temp.no.bio <- subset_samples(temp.physeq, Type !="biofilm")
temp.no.bio <- prune_taxa(taxa_sums(temp.no.bio) > 0, temp.no.bio)

temp.no.wat <- subset_samples(temp.physeq, Type !="water")
temp.no.wat <- prune_taxa(taxa_sums(temp.no.wat) > 0, temp.no.wat)

temp.otu.table <- as.data.frame(otu_table(temp.no.bio))
temp.otu.table <- as.data.frame(t(temp.otu.table))
temp.tax <- as.data.frame(cbind(tax_table(temp.no.bio)))

covariates <- as.character(sample_data(temp.no.bio)$Type)
p1 <- aldex(temp.otu.table, covariates, mc.samples=128, test="t", effect=FALSE, denom="all")
rownames(p1) <- taxa_names(temp.no.bio)
p1 <- subset(p1, we.eBH < 0.05)
colnames(p1) <- paste0("spwat.", colnames(p1))

temp.otu.table <- as.data.frame(otu_table(temp.no.wat))
temp.otu.table <- as.data.frame(t(temp.otu.table))
temp.tax <- as.data.frame(cbind(tax_table(temp.no.wat)))

covariates <- as.character(sample_data(temp.no.wat)$Type)
p2 <- aldex(temp.otu.table, covariates, mc.samples=128, test="t", effect=FALSE, denom="all")
rownames(p2) <- taxa_names(temp.no.wat)
p2 <- subset(p2, we.eBH < 0.05)
colnames(p2) <- paste0("spbio.", colnames(p2))

p.final <- merge(p1, p2, by=0, all=TRUE)
rownames(p.final) <- p.final$Row.names
p.final$Row.names <- NULL

correlation.summary <- merge(correlation.summary, p.final, by=0, all=TRUE)
rownames(correlation.summary) <- correlation.summary$Row.names
correlation.summary$Row.names <- NULL
rm(p1, p2, p.final, temp.no.bio, temp.no.wat, temp.otu.table, temp.tax, temp.physeq)

cor.test(correlation.summary$sponge, correlation.summary$water, method="spearman") # R = 0.384, p < 0.001
cor.test(correlation.summary$sponge, correlation.summary$biofilm, method="spearman") # R = 0.119, p < 0.001

#### BLAST SEQUENCING OF THREE HIGHLY ABUNDANT ASVs ######
# Subset phyloseq object to Sediminibacterium
sediminibacterium <- subset_taxa(transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x)), Genus=="Sediminibacterium")

# Create a data frame with each ASV's mean abundance and DNA sequence.
temp <- as.data.frame(t(otu_table(sediminibacterium)))
sediminibacterium <- data.frame(
  mean_abund = rowMeans(temp),
  refseq = refseq(sediminibacterium))

temp <- as.data.frame(t(temp))
temp$Group <- paste0(sample_data$Source, "_", sample_data$Type)
temp2 <- temp %>% dplyr::group_by(Group) %>% dplyr::summarise_if(is.numeric, mean)
rownames(temp2) <- temp2$Group
temp2$Group <- NULL
sediminibacterium <- cbind(sediminibacterium, as.data.frame(t(temp2)))
rm(temp, temp2)

# Create a FASTA file containing the sequences, named by ASV#.
fasta <- ShortRead(sread = DNAStringSet(sediminibacterium$refseq), id = BStringSet(rownames(sediminibacterium)))

# Export fasta file for a BLAST search and abundance information for later reference.
writeFasta(fasta, file = "sediminibacterium.fasta")
write.csv(sediminibacterium, "sediminibacterium.csv")
rm(sediminibacterium, fasta)

# Repeat for other taxa of interest.
rhodospirillales <- subset_taxa(transform_sample_counts(pseq.rarefied.sub, function(x) 100*x/sum(x)), Family=="Uncl_Rhodospirillales")
temp <- as.data.frame(t(otu_table(rhodospirillales)))
rhodospirillales <- data.frame(
  mean_abund = rowMeans(temp),
  refseq = refseq(rhodospirillales))
temp <- as.data.frame(t(temp))
temp$Group <- paste0(sample_data$Source, "_", sample_data$Type)
temp2 <- temp %>% dplyr::group_by(Group) %>% dplyr::summarise_if(is.numeric, mean)
rownames(temp2) <- temp2$Group
temp2$Group <- NULL
rhodospirillales <- cbind(rhodospirillales, as.data.frame(t(temp2)))
rm(temp, temp2)
fasta <- ShortRead(sread = DNAStringSet(rhodospirillales$refseq), id = BStringSet(rownames(rhodospirillales)))
writeFasta(fasta, file = "rhodospirillales.fasta")
write.csv(rhodospirillales, "rhodospirillales.csv")
rm(rhodospirillales, fasta)

comamonas <- subset_taxa(transform_sample_counts(pseq.rarefied, function(x) 100*x/sum(x)), Genus=="Comamonas")
temp <- as.data.frame(t(otu_table(comamonas)))
comamonas <- data.frame(
  mean_abund = rowMeans(temp),
  refseq = refseq(comamonas))
temp <- as.data.frame(t(temp))
temp$Group <- paste0(sample_data$Source, "_", sample_data$Type)
temp2 <- temp %>% dplyr::group_by(Group) %>% dplyr::summarise_if(is.numeric, mean)
rownames(temp2) <- temp2$Group
temp2$Group <- NULL
comamonas <- cbind(comamonas, as.data.frame(t(temp2)))
rm(temp, temp2)
fasta <- ShortRead(sread = DNAStringSet(comamonas$refseq), id = BStringSet(rownames(comamonas)))
writeFasta(fasta, file = "comamonas.fasta")
write.csv(comamonas, "comamonas.csv")
rm(comamonas, fasta)

##################################################### COMPARISON TO KENNY ET AL. ############################################
#### IMPORT AND PROCESS DATA ####
# Import sequence data from QIIME2 into a phyloseq object
q2.biom <- biomformat::read_biom("raw_data/qiime2_outputs/feature-table.biom")
q2.otus <- as.data.frame(as.matrix(biomformat::biom_data(q2.biom)))
q2.all.seqs <- Biostrings::readDNAStringSet("raw_data/qiime2_outputs/dna-sequences.fasta")

q2.all.data <- read.csv("raw_data/qiime2_outputs/q2.sample_data.csv")
rownames(q2.all.data) <- q2.all.data$SampleID

q2.all.pseq <- phyloseq(otu_table(q2.otus, taxa_are_rows=TRUE),
                        sample_data(q2.all.data),
                        refseq(q2.all.seqs))

rm(q2.biom, q2.otus, q2.all.seqs)

# Create table of OTU IDs with their nucleotide sequences
seqtab.all <- as.matrix(t(otu_table(q2.all.pseq)))
seqs <- refseq(q2.all.pseq)
ids <- taxa_names(q2.all.pseq)
db_all <- data.frame(ids=ids, sequence=seqs, abundance=colSums(seqtab.all))
colnames(seqtab.all) <- db_all$sequence
rm(seqs, ids)

# Examine the distribution of sequence lengths. Most sequences should be ~253bp.
table(nchar(getSequences(db_all)))

# Trim sequence tables to remove amplicons that are not within 250-256bp.
seqtab.all.trim <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(249,256)]

# Subset the sequence "database" to only include the apropriate length-ASVs
db_all_subset <- subset(db_all, sequence %in% colnames(seqtab.all.trim))
nrow(db_all_subset) # 37,438 taxa

# There were a bunch of 378-bp ASVs in the data. Export the sequences for BLAST search to see what they are.
seqtab.long.reads <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(338,338)]
db_longreads <- subset(db_all, sequence %in% colnames(seqtab.long.reads))
write.csv(db_longreads, "db_longreads.csv")

# Identify the reference database you will use for assigning taxonomy.
# This was downloaded from https://benjjneb.github.io/dada2/training.html.
fastaRef <- "~/Documents/AnalysisFolders/FASTQ_INPUT/databases/dada2_rdp_train_set_16.fa.gz"

# Assign taxonomy to your sequence table.
taxtab.all <- assignTaxonomy(db_all_subset, refFasta=fastaRef, multithread=TRUE)

# Put the phyloseq object back together, with taxonomy information and only 250-256bp ASVs 
q2.all.pseq <- phyloseq(otu_table(seqtab.all.trim, taxa_are_rows=FALSE),
                        tax_table(taxtab.all),
                        sample_data(q2.all.data))

rm(db_all, db_all_subset, db_longreads, taxtab.all, fastaRef,
   seqtab.all, seqtab.all.trim, seqtab.long.reads)

# Remove chloroplasts
temp <- as.data.frame(cbind(tax_table(q2.all.pseq)))
temp$chloroplast <- temp$Class == "Chloroplast"
temp <- subset(temp, chloroplast=="TRUE")
badTaxa <- rownames(temp)
goodTaxa <- setdiff(taxa_names(q2.all.pseq), badTaxa)
q2.all.pseq.filter <- prune_taxa(goodTaxa, q2.all.pseq)
rm(temp, goodTaxa, badTaxa) # Removes 441 OTUs

# Remove taxa that are not assigned to a Kingdom or are assigned to archaea.
q2.all.pseq.filter <- subset_taxa(q2.all.pseq.filter, is.na(Kingdom)=="FALSE")
q2.all.pseq.filter <- subset_taxa(q2.all.pseq.filter, !(Kingdom %in% c("Archaea")))

# Store sequences as a refseq object
dna <- Biostrings::DNAStringSet(taxa_names(q2.all.pseq.filter))
names(dna) <- taxa_names(q2.all.pseq.filter)
q2.all.pseq.filter <- merge_phyloseq(q2.all.pseq.filter, dna)
rm(dna)

# Determine which taxa are not assigned to the phylum level for use in a BLAST search.
temp.unclass <- subset_taxa(q2.all.pseq.filter, Phylum !="NA")
temp.unclass <- as.data.frame(cbind(tax_table(temp.unclass)))
class.string <- rownames(temp.unclass)
temp.unclass <- subset_taxa(q2.all.pseq.filter, !(taxa_names(q2.all.pseq.filter) %in% class.string))

# Export the ASV sequences for these taxa. BLAST search to identify closest sequence identity.
seqtab.nochim <- as.matrix(otu_table(temp.unclass))
seqs <- refseq(temp.unclass)
ids <- taxa_names(temp.unclass)
db_out <- data.frame(ids=ids, seqs=seqs, count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
writeFasta(fasta, file = "ASVs_with_no_phylum.fna")
write.csv(db_out, "ASVs_with_no_phylum.csv") # Note 240 sequences with no phylum assignment
rm(fasta, db_out, ids, seqs, temp.unclass, seqtab.nochim, class.string)

# Remove taxa not assigned to the phylum level
# The BLAST results (above) suggested that most of these were noise - chloroplasts, mitochondria, etc.
q2.all.pseq.filter <- subset_taxa(q2.all.pseq.filter, is.na(Phylum)=="FALSE")

# Remove ASVs that have a relative abundance < 0.001%
temp <- transform_sample_counts(q2.all.pseq.filter, function(x) 100*x/sum(x))
x2 <- as.data.frame(colMeans(otu_table(temp)))
temp <- filter_taxa(temp, function(x) mean(x) > 0.001, TRUE)
q2.all.pseq.filter <- subset_taxa(q2.all.pseq.filter, taxa_names(q2.all.pseq.filter) %in% taxa_names(temp))
rm(temp, x2)

# Rename OTUs
taxa_names(q2.all.pseq.filter) <- paste0("OTU", seq(ntaxa(q2.all.pseq.filter))) 

#### ALPHA AND BETA DIVERSITY, ALL SAMPLES ####
# Calculate alpha diversity measures using iNext
temp <- as.data.frame(t(otu_table(q2.all.pseq.filter)))
iNext.result <- iNEXT::iNEXT(temp, q=0, datatype="abundance") # All samples, except the five low-read samples, have completion >98%

temp <- iNext.result$AsyEst
temp.richness <- subset(temp, Diversity=="Species richness")
temp.richness$SampleID <- temp.richness$Site
temp.richness$Observed_Extrap <- temp.richness$Estimator
temp.richness$Observed_Extrap <- round(temp.richness$Observed_Extrap, digits = 0)

temp.diversity <- subset(temp, Diversity=="Shannon diversity")
temp.diversity$SampleID <- temp.diversity$Site
temp.diversity$Shannon_Extrap <- log(temp.diversity$Estimator)

q2.all.data <- merge(q2.all.data, temp.richness[,c("SampleID", "Observed_Extrap")], by="SampleID", all=TRUE)
q2.all.data <- merge(q2.all.data, temp.diversity[,c("SampleID", "Shannon_Extrap")], by="SampleID", all=TRUE)

rm(temp.richness, temp.diversity, temp)

# Create a dummy variable to compare sample type and experiment
q2.all.data$Type <- factor(q2.all.data$Type, levels=c("sponge", "water", "biofilm"))
q2.all.data$TestVariable <- paste0(q2.all.data$Experiment, " ", q2.all.data$Type)

# Test for significant differences
summary(aov(Observed_Extrap ~ TestVariable, q2.all.data))
TukeyHSD(aov(Observed_Extrap ~ TestVariable, q2.all.data))

summary(aov(Shannon_Extrap ~ TestVariable, q2.all.data))
TukeyHSD(aov(Shannon_Extrap ~ TestVariable, q2.all.data))

## Beta diversity (clustering) analysis
# Calculate Aitchison distance between samples
temp.otu <- as.data.frame(t(otu_table(microbiome::transform(q2.all.pseq.filter, "clr"))))
temp.otu <- as.data.frame(t(temp.otu))
temp.dist <- vegdist(temp.otu, method="euclidean")

# Run PCA
temp.pca <- rda(temp.otu)

# Extract species location in the PCA for plotting later
Kenny.all.pca.scores <- as.data.frame(vegan::scores(temp.pca, display="sites", choices=c(1:2)))
Kenny.all.pca.scores$SampleID <- rownames(Kenny.all.pca.scores)

# Add metadata
Kenny.all.pca.scores <- merge(Kenny.all.pca.scores, q2.all.data, by="SampleID", all=FALSE)

Kenny.all.pca.scores$Sample <- paste0(Kenny.all.pca.scores$Stage, " - ", Kenny.all.pca.scores$Experiment)
Kenny.all.pca.scores$Sample <- factor(Kenny.all.pca.scores$Sample,
                             levels=c("sponge - Kenny et al.",
                                      "gemmule - Kenny et al.",
                                      "sponge - This study",
                                      "water - This study",
                                      "biofilm - This study"))
Kenny.all.pca.scores$Source <- Kenny.all.pca.scores$SimpSource

#### ALPHA AND BETA DIVERSITY, SOOKE SAMPLES ONLY ####
# Subset to only Sooke samples for pairwise alpha diversity comparisons
test.data <- subset(q2.all.data, SimpSource=="Sooke")
TukeyHSD(aov(Observed_Extrap ~ TestVariable, test.data))
test.data$SampleID <- factor(test.data$SampleID)

# Subset to only Sooke samples for clustering analysis
temp.pseq <- subset_samples(q2.all.pseq.filter,
                            sample_names(q2.all.pseq.filter) %in% 
                              test.data$SampleID)
temp.pseq <- prune_taxa(taxa_sums(temp.pseq) > 0, temp.pseq)

# Aitchison distance
temp.otu <- as.data.frame(t(otu_table(microbiome::transform(temp.pseq, "clr"))))
temp.otu <- as.data.frame(t(temp.otu))
temp.dist <- vegdist(temp.otu, method="euclidean")

# Principal components analysis
temp.pca <- rda(temp.otu)

# Extract scores
Kenny.SKE.pca.scores <- as.data.frame(vegan::scores(temp.pca, display="sites", choices=c(1:2)))
Kenny.SKE.pca.scores$SampleID <- rownames(Kenny.SKE.pca.scores)

Kenny.SKE.pca.scores <- merge(Kenny.SKE.pca.scores, q2.all.data, by="SampleID", all=FALSE)

Kenny.SKE.pca.scores$Source <- Kenny.SKE.pca.scores$SimpSource
Kenny.SKE.pca.scores$Experiment <- factor(Kenny.SKE.pca.scores$Experiment)

rm(temp.pca, temp.otu, temp.dist, temp.pseq, test.data)

#### DIFFERENTIAL ABUNDANCE ####
# Fill in unnamed taxa levels to prevent them from getting lost during agglomeration
q2.all.pseq.filter.sub <- q2.all.pseq.filter
for(i in 1:nrow(tax_table(q2.all.pseq.filter.sub))){
  for(j in 2:6){
    if(is.na(tax_table(q2.all.pseq.filter.sub)[i,j])==TRUE){
      if(substr(tax_table(q2.all.pseq.filter.sub)[i,j-1], 1, 4)=="Uncl"){
        tax_table(q2.all.pseq.filter.sub)[i,j] <- tax_table(q2.all.pseq.filter.sub)[i,j-1]}
      else {
        tax_table(q2.all.pseq.filter.sub)[i,j] <- paste0("Uncl_", tax_table(q2.all.pseq.filter.sub)[i,j-1])}}
  }}

# Calculate differential abundance for every taxonomic rank
taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")

# Differential abundance between all sponge samples
diff.abund.sally.scott.sponges <- list()
for(i in 1:6){
  temp.physeq <- subset_samples(q2.all.pseq.filter.sub, Type =="sponge")
  if(i < 6){
    temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])  
  }
  if(i > 5){
    temp.physeq <- subset_taxa(temp.physeq, Phylum != "NA")
  }
  temp.otu.table <- as.data.frame(otu_table(temp.physeq))
  temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
  
  if(i < 6){
    colnames(temp.otu.table) <- temp.tax[,i+1] }
  temp.otu.table <- as.data.frame(t(temp.otu.table))
  
  covariates <- as.character(sample_data(temp.physeq)$Experiment)
  p <- aldex(temp.otu.table, covariates, mc.samples=128, test="t", effect=FALSE, denom="all")
  
  temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.transform.otu <- as.data.frame(otu_table(temp.transform))
  
  if(i < 6){
    colnames(temp.transform.otu) <- temp.tax[,i+1] }
  
  temp.transform.otu <- cbind(sample_data(temp.physeq)$Experiment, temp.transform.otu)
  colnames(temp.transform.otu)[1] <- "Experiment"
  
  vector.means <- data.frame(matrix(nrow=7))
  for(j in 2:(ncol(temp.transform.otu))){
    poo1 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Experiment), mean)
    poo1$Group.1 <- NULL
    rownames(poo1) <- levels(temp.transform.otu$Experiment)
    rownames(poo1) <- paste("mean", rownames(poo1))
    
    poo2 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Experiment), median)
    poo2$Group.1 <- NULL
    rownames(poo2) <- levels(temp.transform.otu$Experiment)
    rownames(poo2) <- paste("median", rownames(poo2))
    
    poo3 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Experiment), sd)
    poo3$Group.1 <- NULL
    rownames(poo3) <- levels(temp.transform.otu$Experiment)
    rownames(poo3) <- paste("stdev", rownames(poo3))
    
    # Effect size
    eff <- effsize::cohen.d(temp.transform.otu[,j] ~ temp.transform.otu$Experiment,
                            hedges.correction=TRUE)$estimate
    
    poo <- rbind(poo1, poo2, poo3, eff)
    colnames(poo) <- colnames(temp.transform.otu)[j]
    
    vector.means <- cbind(vector.means, poo)
  }
  
  vector.means[,1] <- NULL
  vector.means <- as.data.frame(t(vector.means))
  vector.means <- cbind(temp.tax, vector.means)
  if(i < 6){
    vector.means <- subset(vector.means, Order !="Uncl_Actinobacteria" | is.na(Order)=="TRUE")
    rownames(vector.means) <- vector.means[,i+1]}
  final <- merge(p, vector.means, by=0, all=TRUE)
  rownames(final) <- final$Row.names
  final$Row.names <- NULL
  final$logFC.sally.scott <- log(final$`mean Sally`/final$`mean Scott`)
  
  diff.abund.sally.scott.sponges[[i]] <- final
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq)
}

names(diff.abund.sally.scott.sponges) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# Differential abundance only for Sooke River samples #
diff.abund.sally.scott.SKE.sponges <- list()
for(i in 1:6){
  temp.physeq <- subset_samples(q2.all.pseq.filter.sub, Type=="sponge" & SimpSource=="Sooke")
  if(i < 6){
    temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])  
  }
  if(i > 5){
    temp.physeq <- subset_taxa(temp.physeq, Phylum != "NA")
  }
  temp.otu.table <- as.data.frame(otu_table(temp.physeq))
  temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
  
  if(i < 6){
    colnames(temp.otu.table) <- temp.tax[,i+1] }
  temp.otu.table <- as.data.frame(t(temp.otu.table))
  
  covariates <- as.character(sample_data(temp.physeq)$Experiment)
  p <- aldex(temp.otu.table, covariates, mc.samples=128, test="t", effect=FALSE, denom="all")
  
  temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp.transform.otu <- as.data.frame(otu_table(temp.transform))
  
  if(i < 6){
    colnames(temp.transform.otu) <- temp.tax[,i+1] }
  
  temp.transform.otu <- cbind(sample_data(temp.physeq)$Experiment, temp.transform.otu)
  colnames(temp.transform.otu)[1] <- "Experiment"
  
  vector.means <- data.frame(matrix(nrow=7))
  for(j in 2:(ncol(temp.transform.otu))){
    poo1 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Experiment), mean)
    poo1$Group.1 <- NULL
    rownames(poo1) <- levels(temp.transform.otu$Experiment)
    rownames(poo1) <- paste("mean", rownames(poo1))
    
    poo2 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Experiment), median)
    poo2$Group.1 <- NULL
    rownames(poo2) <- levels(temp.transform.otu$Experiment)
    rownames(poo2) <- paste("median", rownames(poo2))
    
    poo3 <- aggregate(temp.transform.otu[,j], list(temp.transform.otu$Experiment), sd)
    poo3$Group.1 <- NULL
    rownames(poo3) <- levels(temp.transform.otu$Experiment)
    rownames(poo3) <- paste("stdev", rownames(poo3))
    
    # Effect size
    eff <- effsize::cohen.d(temp.transform.otu[,j] ~ temp.transform.otu$Experiment,
                            hedges.correction=TRUE)$estimate
    
    poo <- rbind(poo1, poo2, poo3, eff)
    colnames(poo) <- colnames(temp.transform.otu)[j]
    
    vector.means <- cbind(vector.means, poo)
  }
  
  vector.means[,1] <- NULL
  vector.means <- as.data.frame(t(vector.means))
  vector.means <- cbind(temp.tax, vector.means)
  if(i < 6){
    vector.means <- subset(vector.means, Order !="Uncl_Actinobacteria" | is.na(Order)=="TRUE")
    rownames(vector.means) <- vector.means[,i+1]}
  final <- merge(p, vector.means, by=0, all=TRUE)
  rownames(final) <- final$Row.names
  final$Row.names <- NULL
  final$logFC.sally.scott <- log(final$`mean Sally`/final$`mean Scott`)
  
  diff.abund.sally.scott.SKE.sponges[[i]] <- final
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq,
     covariates, temp.otu.table, temp.tax, eff)
}

names(diff.abund.sally.scott.SKE.sponges) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

##################################################### METAGENOME ANALYSIS I: COGS ###########################################
#### EXPLORE TAXONOMIC INFORMATION IN RELATION TO AMPLICONS ####
# Import taxonomic data. #
metagenome.otu <- read.csv("raw_data/metagenome_otu_counts_table.csv", row.names="OTU_ID")
metagenome.taxa <- read.csv("raw_data/metagenome_taxonomy_table.csv", row.names="OTU_ID")
metagenome.data <- read.csv("raw_data/metagenome_sample_data.csv")
rownames(metagenome.data) <- metagenome.data$Sample

metagenome.otu <- otu_table(metagenome.otu, taxa_are_rows=TRUE)
metagenome.taxa <- tax_table(as.matrix(metagenome.taxa))

metagenome.pseq <- phyloseq(otu_table(metagenome.otu, taxa_are_rows=TRUE),
                            sample_data(metagenome.data),
                            tax_table(metagenome.taxa))

rm(metagenome.otu, metagenome.taxa)

# Check kingdom-level abundances to see how much archaea, viruses, and viroids contribute to the sample
otu_table(transform_sample_counts(tax_glom(metagenome.pseq, taxrank="Kingdom"), 
                                  function(x) 100*x/sum(x)))
tax_table(transform_sample_counts(tax_glom(metagenome.pseq, taxrank="Kingdom"), 
                                  function(x) 100*x/sum(x)))

# Archaea, viruses, and viroids are a small fraction of the sample, so remove theme from subsequent analyses
metagenome.pseq <- subset_taxa(metagenome.pseq, Kingdom == "Bacteria")

#### IMPORT AND PROCESS ABUNDANCE DATA (RAW COUNTS, RPKM) ####
# a) Import contig abundances extracted from BAM file ####
read_counts_SKE1 <- read.table("raw_data/metagenome_cog_counts/EmuEnvSum1_S7.txt", sep="\t")
read_counts_SKE4 <- read.table("raw_data/metagenome_cog_counts/EmuEnvSum4_S8.txt", sep="\t")
read_counts_SKE5 <- read.table("raw_data/metagenome_cog_counts/EmuEnvSum5_S9.txt", sep="\t")

read_counts_SKEW1 <- read.table("raw_data/metagenome_cog_counts/WatSum1_S10.txt", sep="\t")
read_counts_SKEW4 <- read.table("raw_data/metagenome_cog_counts/WatSum4_S11.txt", sep="\t")
read_counts_SKEW5 <- read.table("raw_data/metagenome_cog_counts/WatSum5_S12.txt", sep="\t")

colnames(read_counts_SKE1) <- c("contig", "length", "SKE1", "misses")
colnames(read_counts_SKE4) <- c("contig", "length", "SKE4", "misses")
colnames(read_counts_SKE5) <- c("contig", "length", "SKE5", "misses")

colnames(read_counts_SKEW1) <- c("contig", "length", "SKEW1", "misses")
colnames(read_counts_SKEW4) <- c("contig", "length", "SKEW4", "misses")
colnames(read_counts_SKEW5) <- c("contig", "length", "SKEW5", "misses")

metagenome.contig_counts <- cbind(read_counts_SKE1[,c("contig","length","SKE1")],
                                  read_counts_SKE4[,c("SKE4")],
                                  read_counts_SKE5[,c("SKE5")],
                                  read_counts_SKEW1[,c("SKEW1")],
                                  read_counts_SKEW4[,c("SKEW4")],
                                  read_counts_SKEW5[,c("SKEW5")])

colnames(metagenome.contig_counts) <- c("Contig", "Length", "SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")
metagenome.contig_counts <- subset(metagenome.contig_counts, Contig !="*")
head(metagenome.contig_counts)

rm(read_counts_SKE1, read_counts_SKE4, read_counts_SKE5,
   read_counts_SKEW1, read_counts_SKEW4, read_counts_SKEW5)

# Convert raw counts to RPKM values for each contig.
metagenome.contig_rpkm <- metagenome.contig_counts
for(i in c("SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")){
  df <- metagenome.contig_rpkm[,c(i,"Length")]
  df[,1] <- df[,1] / ((df[,2]/1000) * (sum(df[,1]/1000000)))
  metagenome.contig_rpkm[,i] <- df[,1]
  rm(df)
}

# b) Convert contig abundances to gene abundances (gene abundance = abundance of the contig on which it was found) ####
# Import and process BLAST data for gene identities.
metagenome.gene_blast <- read.csv("raw_data/metagenome_cog_counts/contigs_to_cogs_webmga.txt", sep="\t", header=FALSE)
colnames(metagenome.gene_blast) <- c("QueryID", "SubjectID", "Percent.Identity", "Length", "Mismatches", "GapOpenings", 
                                     "QueryStart", "QueryEnd", "SubjectStart", "SubjectEnd", "EValue", "Bitscore",
                                     "COG", "Gene", "Name")
metagenome.gene_blast <- tidyr::separate(metagenome.gene_blast, col=QueryID, into=c("ContigA", "ContigB", "GeneNumber"),
                                         sep="_", remove=TRUE, convert=FALSE)
metagenome.gene_blast$Contig <- paste0(metagenome.gene_blast$ContigA, "_", metagenome.gene_blast$ContigB)
metagenome.gene_blast$ContigA <- NULL
metagenome.gene_blast$ContigB <- NULL
head(metagenome.gene_blast)

# Use contig counts as a proxy for gene counts.
metagenome.gene_counts <- merge(metagenome.gene_blast[,c("Contig","GeneNumber","Percent.Identity","Length","COG","Gene","Name")],
                                metagenome.contig_counts[,c("Contig","SKE1","SKE4","SKE5","SKEW1","SKEW4","SKEW5")],
                                by="Contig", all=FALSE)
colnames(metagenome.gene_counts) <- c("Contig", "GeneNumber", "Perc.Identity", "Length", "COG", "Gene", "Name",
                                      "SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")
head(metagenome.gene_counts, 10)

# Use contig RPKM as a proxy for gene RPKM
metagenome.gene_rpkm <- metagenome.gene_counts
metagenome.gene_rpkm[,c("SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")] <- NULL
metagenome.gene_rpkm <- merge(metagenome.gene_rpkm,
                              metagenome.contig_rpkm,
                              by="Contig", all=FALSE)

# c) Calculate COG abundances by summing all genes that map to the same COG ####
# Obtain COG abundances by summing abundances of each COG.
metagenome.cog_counts <- metagenome.gene_counts %>% dplyr::group_by(COG) %>% dplyr::summarise(SKE1=sum(SKE1, na.rm=TRUE),
                                                                                              SKE4=sum(SKE4, na.rm=TRUE),
                                                                                              SKE5=sum(SKE5, na.rm=TRUE),
                                                                                              SKEW1=sum(SKEW1, na.rm=TRUE),
                                                                                              SKEW4=sum(SKEW4, na.rm=TRUE),
                                                                                              SKEW5=sum(SKEW5, na.rm=TRUE))
head(metagenome.cog_counts)

# Obtain COG RPKM by summing RPKMs of each gene assigned to the same COG. Calculate means for sponges and water.
metagenome.cog_rpkm <- metagenome.gene_rpkm %>% dplyr::group_by(COG) %>% dplyr::summarise(SKE1=sum(SKE1, na.rm=TRUE),
                                                                                          SKE4=sum(SKE4, na.rm=TRUE),
                                                                                          SKE5=sum(SKE5, na.rm=TRUE),
                                                                                          SKEW1=sum(SKEW1, na.rm=TRUE),
                                                                                          SKEW4=sum(SKEW4, na.rm=TRUE),
                                                                                          SKEW5=sum(SKEW5, na.rm=TRUE))

metagenome.cog_rpkm$Sponge <- rowMeans(metagenome.cog_rpkm[,c("SKE1", "SKE4", "SKE5")])
metagenome.cog_rpkm$Water <- rowMeans(metagenome.cog_rpkm[,c("SKEW1", "SKEW4", "SKEW5")])

# Create annotations for the COGs.
metagenome.cog_annotations <- dplyr::distinct(metagenome.gene_counts[,c("COG","Gene","Name")])

# This file lists the COG class (A, B, C, D, etc.) that each COG (COG0001, COG0002, COG0003, etc.) is in
cog_to_class <- read.csv("raw_data/cogid_contigs_class_key.txt", sep="\t")
cog_to_class <- dplyr::distinct(cog_to_class[,c("COG","Class")])
metagenome.cog_annotations <- merge(metagenome.cog_annotations,
                                    cog_to_class[,c("COG","Class")],
                                    by="COG", all=TRUE)
rownames(metagenome.cog_annotations) <- metagenome.cog_annotations$COG
rm(cog_to_class)

# Remove eukaryotic proteins!
metagenome.cog_annotations <- subset(metagenome.cog_annotations, Class !="#N/A")
metagenome.cog_annotations$Class <- factor(metagenome.cog_annotations$Class)
metagenome.cog_annotations$COG <- factor(metagenome.cog_annotations$COG)

metagenome.cog_counts <- subset(metagenome.cog_counts, COG %in% metagenome.cog_annotations$COG)
metagenome.cog_counts$COG <- factor(metagenome.cog_counts$COG)

metagenome.cog_rpkm <- subset(metagenome.cog_rpkm, COG %in% metagenome.cog_annotations$COG)
metagenome.cog_rpkm$COG <- factor(metagenome.cog_rpkm$COG)

# d) Calculate COG class abundances (RPKM only) ####
temp <- merge(metagenome.cog_rpkm,
              metagenome.cog_annotations,
              by="COG", all=TRUE)
temp <- temp %>% dplyr::group_by(Class) %>% dplyr::summarise_if(is.numeric, sum)

# A bunch of COGs are assigned to multiple classes
# To ensure they get counted in each class, separate classes into different columns
temp <- tidyr::separate(temp, col=Class, into=c("Blank", "Class1", "Class2","Class3","Class4"),
                        sep="", remove=TRUE, convert=FALSE)
temp$Blank <- NULL

temp1 <- temp %>% dplyr::group_by(Class1) %>% dplyr::summarise_if(is.numeric, sum) %>% as.data.frame()
temp2 <- temp %>% dplyr::group_by(Class2) %>% dplyr::summarise_if(is.numeric, sum) %>% as.data.frame()
temp3 <- temp %>% dplyr::group_by(Class3) %>% dplyr::summarise_if(is.numeric, sum) %>% as.data.frame()
temp4 <- temp %>% dplyr::group_by(Class4) %>% dplyr::summarise_if(is.numeric, sum) %>% as.data.frame()

temp1$Class1 <- factor(temp1$Class1)
temp2$Class2 <- factor(temp2$Class2)
temp3$Class3 <- factor(temp3$Class3)
temp4$Class4 <- factor(temp4$Class4)

temp1 <- subset(temp1, is.na(Class1)=="FALSE")
temp2 <- subset(temp2, is.na(Class2)=="FALSE")
temp3 <- subset(temp3, is.na(Class3)=="FALSE")
temp4 <- subset(temp4, is.na(Class4)=="FALSE")

rownames(temp1) <- temp1$Class1
rownames(temp2) <- temp2$Class2
rownames(temp3) <- temp3$Class3
rownames(temp4) <- temp4$Class4

# Count the total RPKM per each COG class
temp <- data.frame(matrix(nrow=25, ncol=7))

for(i in 1:nlevels(temp1$Class1)){
  temp[i,1] <- levels(temp1$Class1)[i]
  
  df <- rbind(temp1[levels(temp1$Class1)[i],c(2:7)],
              temp2[levels(temp1$Class1)[i],c(2:7)],
              temp3[levels(temp1$Class1)[i],c(2:7)],
              temp4[levels(temp1$Class1)[i],c(2:7)])
  
  temp[i,c(2:7)] <- colSums(df, na.rm=TRUE)
  rm(df)
}

metagenome.class_rpkm <- temp
colnames(metagenome.class_rpkm) <- c("Class", "SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")
rm(temp, temp1, temp2, temp3, temp4)

# e) Ordination based on COG RPKM abundances to check for visual differences ####
metagenome.sample_data <- read.csv("raw_data/metagenome_sample_data.csv")

# Log-transform RPKM values before ordination.
cog_comm_matrix <- as.data.frame(t(metagenome.cog_rpkm[,c("SKE1","SKE4","SKE5", "SKEW1", "SKEW4", "SKEW5")]))
cog_comm_matrix <- log(cog_comm_matrix+0.01)
cog_ord_PCA <- rda(cog_comm_matrix)
summary(cog_ord_PCA) # 79.48% on Axis 1, 8.99% on Axis 2

metagenome.sample_data <- cbind(metagenome.sample_data,
                                scores(cog_ord_PCA, display="sites", choices = c(1:3)))
colnames(metagenome.sample_data)[c(7:9)] <- c("COG_logPC1", "COG_logPC2", "COG_logPC3")

ggplot(data = metagenome.sample_data, aes(COG_logPC1, COG_logPC2)) + 
  geom_point(aes(color = Type, shape=Type), size=1.75) +
  scale_shape_manual(values=c(17,16)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.5),
        axis.text=element_text(color="black", size = 9, vjust=1),
        axis.title=element_text(size=9),
        strip.background=element_blank(), 
        strip.text=element_text(size=9),
        legend.title=element_blank(),
        plot.tag=element_text(face="bold"),
        legend.position="right",
        legend.justification=c(0.5,0.5),
        plot.title=element_text(size=11, color="black")) +
  scale_color_manual(values=c("maroon2", "orange2")) +
  labs(x="PC1 (79.5%)", y="PC2 (9.00%)", title="PCA - COG abundances (log RPKM)")

rm(cog_comm_matrix, cog_ord_PCA)

#### DIFFERENTIAL ABUNDANCE TESTING BY COG ####
# a) Perform generic differential abundance tests ####
library(edgeR)

dgList <- DGEList(counts=metagenome.cog_counts[,c("SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")],
                  genes=metagenome.cog_counts$COG)
designMat <- model.matrix(~Type, data=data.frame(Type=c("sponge","sponge","sponge","water","water","water")))

dgList <- calcNormFactors(dgList, method="TMM")
plotMDS(dgList)
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
plotBCV(dgList)

fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=2)
edgeR_result <- topTags(lrt)

metagenome.cog_diff_abund <- lrt[["table"]]
metagenome.cog_diff_abund$p.adj <- p.adjust(metagenome.cog_diff_abund$PValue, method="BH")
metagenome.cog_diff_abund$log10q <- -1*log(metagenome.cog_diff_abund$p.adj, base=10)

rm(fit, lrt, edgeR_result, dgList, designMat)

rownames(metagenome.cog_diff_abund) <- metagenome.cog_annotations$COG
metagenome.cog_diff_abund$COG <- rownames(metagenome.cog_diff_abund)

ggplot(metagenome.cog_diff_abund, aes(x=logFC, y=log10q)) + geom_point() # Explore a quick volcano plot

# Add annotations and RPKM data to the differential abundance data
metagenome.cog_diff_abund <- merge(metagenome.cog_diff_abund,
                                   metagenome.cog_rpkm,
                                   by="COG", all=TRUE)
metagenome.cog_diff_abund <- merge(metagenome.cog_diff_abund,
                                   metagenome.cog_annotations,
                                   by="COG", all=TRUE)

metagenome.cog_diff_abund <- metagenome.cog_diff_abund[,c("COG","Class","Gene","Name","logFC","p.adj","log10q",
                                                          "Sponge","Water","SKE1","SKE4","SKE5","SKEW1","SKEW4","SKEW5")]

nrow(subset(metagenome.cog_diff_abund, p.adj < 0.01)) # 2,479 differentially abundant genes at p < 0.01 (60%!)

nrow(subset(metagenome.cog_diff_abund, p.adj < 0.01 & logFC < -1)) # 924 genes enriched in sponges
nrow(subset(metagenome.cog_diff_abund, p.adj < 0.01 & logFC > 1)) # 493 genes enriched in water

# b) Visually compare COG RPKM b/w sponges and water, colored by differential abundance ####
metagenome.cog_compare <- metagenome.cog_diff_abund
metagenome.cog_compare$SpongePseudo <- metagenome.cog_compare$Sponge+0.01
metagenome.cog_compare$WaterPseudo <- metagenome.cog_compare$Water+0.01

metagenome.cog_compare$Enrich <- rep("none")
for(i in 1:nrow(metagenome.cog_compare)){
  if(metagenome.cog_compare[i,"logFC"] < -1 & metagenome.cog_compare[i,"p.adj"] < 0.01){
    metagenome.cog_compare[i,"Enrich"] <- "sponge"
  }
  if(metagenome.cog_compare[i,"logFC"] > 1 & metagenome.cog_compare[i,"p.adj"] < 0.01){
    metagenome.cog_compare[i,"Enrich"] <- "water"
  }
}
metagenome.cog_compare$Enrich <- factor(metagenome.cog_compare$Enrich)
levels(metagenome.cog_compare$Enrich)

ggplot(metagenome.cog_compare, aes(x=WaterPseudo, y=SpongePseudo, color=Enrich)) + geom_point() +
  scale_color_manual(values=c("grey","orange","blue")) +
  scale_y_log10() +
  scale_x_log10()

rm(metagenome.cog_compare)

#### DIFFERENTIAL ABUNDANCE TESTING: HYPERGEOMETRIC OVER-REPRESENTATION ANALYSIS ####
# a) Define a function to count the number of COGs in each class (including for COGs assigned to multiple classes) ####
cog_counter <- function(table){
  temp <- tidyr::separate(table, col=Class, into=c("Blank", "Class1", "Class2","Class3","Class4"),
                          sep="", remove=TRUE, convert=FALSE)
  temp$Blank <- NULL
  temp1 <- temp %>% dplyr::group_by(Class1) %>% dplyr::count()
  temp2 <- temp %>% dplyr::group_by(Class2) %>% dplyr::count()
  temp3 <- temp %>% dplyr::group_by(Class3) %>% dplyr::count()
  temp4 <- temp %>% dplyr::group_by(Class4) %>% dplyr::count()
  colnames(temp1) <- c("Class", "Count1")
  colnames(temp2) <- c("Class", "Count2")
  colnames(temp3) <- c("Class", "Count3")
  colnames(temp4) <- c("Class", "Count4")
  temp <- merge(temp1, temp2, by=c("Class"), all=TRUE)
  temp <- merge(temp, temp3, by=c("Class"), all=TRUE)
  temp <- merge(temp, temp4, by=c("Class"), all=TRUE)
  rm(temp1, temp2, temp3, temp4)
  temp$Count <- rowSums(temp[,c(2:5)], na.rm=TRUE)
  temp[,c(2:5)] <- NULL
  temp <- subset(temp, Class !="NA")
  colnames(temp) <- c("Class", "Total_Count")
  output <- temp
  rm(temp)
  output
}

# b) Summarize the number of COGs in each class, for each sample type (total, diff abund, more abund in sponge, more abund in water) ####
# Calculate the total number of COGs in each class, in the entire data set.
metagenome.cog_counts <- merge(metagenome.cog_counts,
                               metagenome.cog_annotations,
                               by="COG", all=TRUE)
metagenome.class_sums <- cog_counter(metagenome.cog_counts)
metagenome.cog_counts[,c("Gene","Name","Class")] <- NULL

# Calculate the total number of differentially abundant COGs in each class (p < 0.05)
tempA <- subset(metagenome.cog_diff_abund, abs(logFC) > 1 & p.adj < 0.01)
tempA <- cog_counter(tempA)
colnames(tempA) <- c("Class", "All")

# Calculate the total number of differentially abundant COGs in each class that are more abundant in sponges
tempB <- subset(metagenome.cog_diff_abund, logFC < -1 & p.adj < 0.01)
tempB <- cog_counter(tempB)
colnames(tempB) <- c("Class", "Sponge")

# Calculate the total number of differentially abundant COGs in each class that are more abundant in water
tempC <- subset(metagenome.cog_diff_abund, logFC > 1 & p.adj < 0.01)
tempC <- cog_counter(tempC)
colnames(tempC) <- c("Class", "Water")

# Create a summary of these data
metagenome.class_sums <- merge(metagenome.class_sums, tempA, by="Class", all=TRUE)
metagenome.class_sums <- merge(metagenome.class_sums, tempB, by="Class", all=TRUE)
metagenome.class_sums <- merge(metagenome.class_sums, tempC, by="Class", all=TRUE)
rm(tempA, tempB, tempC)

# c) Run the Fisher hypergeometric over-representation test. ####
metagenome.class_fisher_test <- data.frame(matrix(nrow=25, ncol=4))
metagenome.class_fisher_test[,1] <- metagenome.class_sums$Class

metagenome.class_sums[is.na(metagenome.class_sums)] <- 0

totalM <- nrow(subset(metagenome.cog_diff_abund, is.na(Class)=="FALSE")) # Total genes with a COG term (4125)

# First, for over-representation in the entire pool of differentially abundant genes
for(p in 1:nrow(metagenome.class_sums)){
  total_sp <- metagenome.class_sums[p,"Total_Count"] # Total number of genes within that COG
  
  totalDE <- sum(metagenome.class_sums$All)
  total_sp_DE <- metagenome.class_sums[p,"All"] # Total number of differentially expressed genes within that COG
  
  total_sp_notDE <- total_sp - total_sp_DE
  total_DE_notsp <- totalDE - total_sp_DE
  total_notsp_notDE <- totalM - total_sp - total_DE_notsp
  
  metagenome.class_fisher_test[p,2] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), nrow=2), alternative = "greater")$p.value
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE, totalDE)
}

# Second, for over-representation in the pool of genes more abundant in sponges
for(p in 1:nrow(metagenome.class_sums)){
  total_sp <- metagenome.class_sums[p,"Total_Count"] # Total number of genes within that COG
  
  totalDE <- sum(metagenome.class_sums$Sponge)
  total_sp_DE <- metagenome.class_sums[p,"Sponge"] # Total number of differentially expressed genes within that COG
  
  total_sp_notDE <- total_sp - total_sp_DE
  total_DE_notsp <- totalDE - total_sp_DE
  total_notsp_notDE <- totalM - total_sp - total_DE_notsp
  
  metagenome.class_fisher_test[p,3] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), nrow=2), alternative = "greater")$p.value
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE, totalDE)
}

# Third, for over-representation in the pool of genes more abundant in water
for(p in 1:nrow(metagenome.class_sums)){
  total_sp <- metagenome.class_sums[p,"Total_Count"] # Total number of genes within that COG
  
  totalDE <- sum(metagenome.class_sums$Water)
  total_sp_DE <- metagenome.class_sums[p,"Water"] # Total number of differentially expressed genes within that COG
  
  total_sp_notDE <- total_sp - total_sp_DE
  total_DE_notsp <- totalDE - total_sp_DE
  total_notsp_notDE <- totalM - total_sp - total_DE_notsp
  
  metagenome.class_fisher_test[p,4] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), nrow=2), alternative = "greater")$p.value
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE, totalDE)
}

colnames(metagenome.class_fisher_test) <- c("Class", "All", "Sponge", "Water")

metagenome.class_fisher_test$All.adj <- p.adjust(metagenome.class_fisher_test$All, method="BH")
metagenome.class_fisher_test$Sponge.adj <- p.adjust(metagenome.class_fisher_test$Sponge, method="BH")
metagenome.class_fisher_test$Water.adj <- p.adjust(metagenome.class_fisher_test$Water, method="BH")

rm(totalM)

#### EXPLORE DIFFERENTIAL ABUNDANCE BASED ON COGS OF FUNCTIONAL INTEREST ####
# a) Create heatmap/dendrogram combinations for easy visual comparison ####
# Import a data frame that lists sets of COGs based on their function.
# These functions and COGs were identified from an evaluation of previous sponge microbiome papers in the literature.
sorting_data <- read.csv("raw_data/cogs_from_other_papers.csv")

# Create a list of data frames, where each data frame includes the data for COGs with a different function
sorting_list <- list()
for(i in 1:ncol(sorting_data)){
  sorting_list[[i]] <- subset(metagenome.cog_diff_abund, COG %in% sorting_data[,i])
}
names(sorting_list) <- colnames(sorting_data)

# Create a data frame of COG RPKMs with no other data
cog_abund_matr <- as.data.frame(metagenome.cog_rpkm[,c("COG", "SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")])
rownames(cog_abund_matr) <- cog_abund_matr$COG
cog_abund_matr$COG <- NULL

# Create heat maps showing the abundance of each COG in sponge and water samples
# One heat map per functional group
sorting_list_maps <- list()
temp_dendros <- list()

for(i in 1:ncol(sorting_data)){
  cog_abund_matr_subset <- subset(cog_abund_matr, rownames(cog_abund_matr) %in% sorting_data[,i])
  cog_abund_matr_subset <- as.data.frame(t(subset(cog_abund_matr_subset)))
  x <- hclust(vegdist(cog_abund_matr_subset, method="bray"))
  temp_dendros[[i]] <- x
  y <- x$order
  y[y==1] <- "SKE1"
  y[y==2] <- "SKE4"
  y[y==3] <- "SKE5"
  y[y==4] <- "SKEW1"
  y[y==5] <- "SKEW4"
  y[y==6] <- "SKEW5"
  rm(x)
  
  temp_dendros[[i]] <- as.dendrogram(temp_dendros[[i]])
  temp_dendros[[i]] <- ggdendro::dendro_data(temp_dendros[[i]])
  
  cog_abund_matr_subset$SampleID <- rownames(cog_abund_matr_subset)
  
  cog_abund_melt <- reshape2::melt(cog_abund_matr_subset)
  cog_abund_melt$variable <- factor(cog_abund_melt$variable)
  
  cog_abund_melt$value[cog_abund_melt$value > 100000] <- 100000
  cog_abund_melt$SampleID <- factor(cog_abund_melt$SampleID, levels=c(y))
  rm(y)
  
  colnames(cog_abund_melt) <- c("SampleID", "COG", "value")
  cog_abund_melt <- merge(cog_abund_melt, metagenome.cog_annotations[,c("COG", "Name")], by="COG", all=FALSE)
  cog_abund_melt$Full.Name <- paste0(cog_abund_melt$COG, " ", cog_abund_melt$Name)
  cog_abund_melt$Full.Name <- stringr::str_trunc(cog_abund_melt$Full.Name, 50)
  cog_abund_melt$Full.Name <- as.factor(cog_abund_melt$Full.Name)
  
  temp <- ggplot(cog_abund_melt, aes(x=SampleID, y=Full.Name)) + 
    geom_tile(aes(fill=value)) +
    scale_fill_gradient2() +
    scale_y_discrete(limits=rev(levels(cog_abund_melt$Full.Name)), position="right") +
    labs(x="\nSampleID", y="COG", fill="RPKM", title=colnames(sorting_data)[i]) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y=element_text(color="black", size = 10),
          axis.text.x=element_text(color="black", size = 10, angle=90, vjust=0.5),
          axis.title=element_text(size=11),
          plot.tag=element_text(face="bold"),
          legend.title = element_text(size=11, color="black"),
          legend.position="bottom",
          legend.justification=c(0.5,0.5),
          plot.title=element_text(size=11, color="black"))
  
  sorting_list_maps[[i]] <- temp
  rm(temp)
}
names(sorting_list_maps) <- colnames(sorting_data)

# Add dendrograms (must be done manually) 
sorting_list_maps[[1]] <- sorting_list_maps[[1]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[1]]), aes(x=x,
                                                              y=7*y+39+0.5,
                                                              xend=xend,
                                                              yend=7*yend+39+0.5))
sorting_list_maps[[2]] <- sorting_list_maps[[2]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[2]]), aes(x=x,
                                                              y=7*y+63+0.5,
                                                              xend=xend,
                                                              yend=7*yend+63+0.5))
sorting_list_maps[[3]] <- sorting_list_maps[[3]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[3]]), aes(x=x,
                                                              y=7*y+88+0.5,
                                                              xend=xend,
                                                              yend=7*yend+88+0.5))
sorting_list_maps[[4]] <- sorting_list_maps[[4]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[4]]), aes(x=x,
                                                              y=7*y+39+0.5,
                                                              xend=xend,
                                                              yend=7*yend+39+0.5))
sorting_list_maps[[5]] <- sorting_list_maps[[5]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[5]]), aes(x=x,
                                                              y=7*y+33+0.5,
                                                              xend=xend,
                                                              yend=7*yend+33+0.5))
sorting_list_maps[[6]] <- sorting_list_maps[[6]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[6]]), aes(x=x,
                                                              y=7*y+17+0.5,
                                                              xend=xend,
                                                              yend=7*yend+17+0.5))
sorting_list_maps[[7]] <- sorting_list_maps[[7]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[7]]), aes(x=x,
                                                              y=7*y+151+0.5,
                                                              xend=xend,
                                                              yend=7*yend+151+0.5))
sorting_list_maps[[8]] <- sorting_list_maps[[8]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[8]]), aes(x=x,
                                                              y=7*y+47+0.5,
                                                              xend=xend,
                                                              yend=7*yend+47+0.5))
sorting_list_maps[[9]] <- sorting_list_maps[[9]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[9]]), aes(x=x,
                                                              y=7*y+25+0.5,
                                                              xend=xend,
                                                              yend=7*yend+25+0.5))
sorting_list_maps[[10]] <- sorting_list_maps[[10]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[10]]), aes(x=x,
                                                               y=7*y+19+0.5,
                                                               xend=xend,
                                                               yend=7*yend+19+0.5))
sorting_list_maps[[11]] <- sorting_list_maps[[11]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[11]]), aes(x=x,
                                                               y=7*y+20+0.5,
                                                               xend=xend,
                                                               yend=7*yend+20+0.5))
sorting_list_maps[[12]] <- sorting_list_maps[[12]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[12]]), aes(x=x,
                                                               y=7*y+17+0.5,
                                                               xend=xend,
                                                               yend=7*yend+17+0.5))
sorting_list_maps[[13]] <- sorting_list_maps[[13]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[13]]), aes(x=x,
                                                               y=7*y+26+0.5,
                                                               xend=xend,
                                                               yend=7*yend+26+0.5))
sorting_list_maps[[14]] <- sorting_list_maps[[14]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[14]]), aes(x=x,
                                                               y=7*y+79+0.5,
                                                               xend=xend,
                                                               yend=7*yend+79+0.5))
sorting_list_maps[[15]] <- sorting_list_maps[[15]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[15]]), aes(x=x,
                                                               y=7*y+21+0.5,
                                                               xend=xend,
                                                               yend=7*yend+21+0.5))
sorting_list_maps[[16]] <- sorting_list_maps[[16]] + 
  geom_segment(data=ggdendro::segment(temp_dendros[[16]]), aes(x=x,
                                                               y=7*y+15+0.5,
                                                               xend=xend,
                                                               yend=7*yend+15+0.5))

rm(cog_abund_matr, cog_abund_matr_subset, cog_abund_melt, temp_dendros)

# b) Visualize heatmaps and dendrograms ####
sorting_list_maps$Ribosome
sorting_list_maps$Flagella
sorting_list_maps$Pilus
sorting_list_maps$CRISPR
sorting_list_maps$ABC.Transport
sorting_list_maps$Secretion
sorting_list_maps$Transposase
sorting_list_maps$FOGs
sorting_list_maps$Restriction.Modification
sorting_list_maps$Cobalamin
sorting_list_maps$Spore.Capsule
sorting_list_maps$Drugs

#### EXPLORE COG DIFFERENTIAL ABUNDANCE b/w SPONGES AND WATER - SEPARATED BY TAXONOMIC GROUPS ####
# Import contig taxonomy data 
metagenome.contig_taxonomy <- read.csv("raw_data/metagenome_contig_taxonomy.csv")
metagenome.contig_taxonomy <- dplyr::distinct(metagenome.contig_taxonomy)
head(metagenome.contig_taxonomy)
metagenome.contig_taxonomy[,c(10:16)] <- NULL
colnames(metagenome.contig_taxonomy)[1] <- "Contig"

metagenome.cog_taxonomy <- merge(metagenome.contig_taxonomy,
                                 metagenome.gene_blast[,c("COG","Contig")],
                                 by="Contig",
                                 all=FALSE)

rm(metagenome.contig_taxonomy)

# Clean taxon information for clustering
metagenome.cog_taxonomy <- subset(metagenome.cog_taxonomy, domain=="Bacteria")
metagenome.cog_taxonomy <- dplyr::distinct(metagenome.cog_taxonomy)

# Convert factors to characters
metagenome.cog_taxonomy$phylum <- as.character(metagenome.cog_taxonomy$phylum)
metagenome.cog_taxonomy$class <- as.character(metagenome.cog_taxonomy$class)
metagenome.cog_taxonomy$order <- as.character(metagenome.cog_taxonomy$order)
metagenome.cog_taxonomy$family <- as.character(metagenome.cog_taxonomy$family)
metagenome.cog_taxonomy$genus <- as.character(metagenome.cog_taxonomy$genus)
metagenome.cog_taxonomy$species <- as.character(metagenome.cog_taxonomy$species)

# Replace missing taxa with the next highest classification
for(i in c(5:8)){
  for(j in c(1:nrow(metagenome.cog_taxonomy))){
    if(is.na(metagenome.cog_taxonomy[j,i])=="TRUE"){
      metagenome.cog_taxonomy[j,i] <- metagenome.cog_taxonomy[j,i-1]
    }
  }
}

# Add COG names and counts/RPKM data to the data frame
metagenome.cog_annotations$Full.Name <- paste0(metagenome.cog_annotations$COG, " ", metagenome.cog_annotations$Name)
metagenome.cog_annotations$Full.Name <- stringr::str_trunc(metagenome.cog_annotations$Full.Name, 50)
metagenome.cog_annotations$Full.Name <- as.factor(metagenome.cog_annotations$Full.Name)

metagenome.cog_taxonomy <- merge(metagenome.cog_taxonomy,
                                 metagenome.cog_annotations[,c("COG", "Name", "Class", "Full.Name")],
                                 by="COG", all=FALSE)

# Add notations for interesting COGs
# For simplicity, use only the COGs that appear in figures, which are a subset of all the COGs
# 'sorting_plot' is originally generated in the sp_figures_2022 code file
sorting_plot <- read.csv("raw_data/sorted_cogs_for_figures.csv")
sorting_plot$Group <- factor(sorting_plot$Group)
sorting_plot$Full.Name <- factor(sorting_plot$Full.Name)

merge_data <- sorting_plot
merge_data$value <- NULL
merge_data$variable <- NULL
merge_data <- dplyr::distinct(merge_data)

metagenome.cog_taxonomy <- merge(metagenome.cog_taxonomy,
                                 merge_data[,c("COG", "Group")],
                                 by="COG", all.x=TRUE, all.y=FALSE)

metagenome.cog_taxonomy <- dplyr::distinct(metagenome.cog_taxonomy)
metagenome.cog_taxonomy <- subset(metagenome.cog_taxonomy, is.na(phylum)=="FALSE")
metagenome.cog_taxonomy <- subset(metagenome.cog_taxonomy, phylum != "Other")


# Add RPKM information
metagenome.cog_taxonomy_rpkm <- merge(metagenome.cog_taxonomy,
                                      metagenome.contig_rpkm,
                                      by="Contig", all=FALSE)
metagenome.cog_taxonomy_rpkm$Length <- NULL

metagenome.cog_taxonomy_rpkm <- subset(metagenome.cog_taxonomy_rpkm, is.na(Group)=="FALSE")

# Group COGs to the family level and subset to only differentially abundant COGs.
metagenome.cog_taxonomy_rpkm_group <- metagenome.cog_taxonomy_rpkm %>%
  dplyr::group_by(COG, Full.Name, phylum, class, order, family) %>%
  dplyr::summarise_if(is.numeric, sum)

metagenome.cog_taxonomy_rpkm_group <- merge(metagenome.cog_taxonomy_rpkm_group,
                                            merge_data[,c("COG", "Group")],
                                            by="COG", all.x=TRUE, all.y=FALSE)

# Calculate the mean RPKM of each COG/contig in sponges and water samples
metagenome.cog_taxonomy_rpkm_group$mean.sponge <- rowMeans(metagenome.cog_taxonomy_rpkm_group[,c("SKE1", "SKE4", "SKE5")])
metagenome.cog_taxonomy_rpkm_group$mean.water <- rowMeans(metagenome.cog_taxonomy_rpkm_group[,c("SKEW1", "SKEW4", "SKEW5")])

rm(merge_data)

# This data is then used for a heat map (Fig. S18) to visually explore which taxa carry which contigs

##################################################### METAGENOME ANALYSIS II: MAGS ###########################################
#### IMPORT DATA AND RUN AN EXPLORATORY PCA ####
# Import MAG data. The 'master' data file includes reference genomes for the sponge-associated MAGs.
mag.master_sample_data <- read.csv("raw_data/mag.master_sample_data.csv")

# Create a data frame without the reference genomes.
mag.sample_data <- subset(mag.master_sample_data, Type !="reference")
mag.sample_data$MAG <- factor(mag.sample_data$MAG)
mag.sample_data <- mag.sample_data[order(mag.sample_data$MAG), ]
rownames(mag.sample_data) <- mag.sample_data$MAG

# Calculate mean relative abundance of each MAG in each sample type
mag.sample_data$mean.rel.sponge <- rowMeans(mag.sample_data[,c("SKE1","SKE4","SKE5")])
mag.sample_data$mean.rel.water <- rowMeans(mag.sample_data[,c("SKEW1","SKEW4","SKEW5")])
mag.sample_data$logFC.sponge <- log(mag.sample_data$mean.rel.sponge / mag.sample_data$mean.rel.water)

# Import the genome content for each MAG (number of COGs in each MAG)
mag.cog_counts <- read.csv("raw_data/mag.cog_counts.csv")
mag.cog_counts$Group <- NULL
mag.cog_counts$Name <- NULL

mag.cog_counts <- mag.cog_counts %>% dplyr::group_by(Bin, COG) %>% dplyr::count()
mag.cog_counts <- reshape2::dcast(mag.cog_counts, COG ~ Bin)
rownames(mag.cog_counts) <- mag.cog_counts$COG
mag.cog_counts$COG <- NULL
mag.cog_counts[is.na(mag.cog_counts)=="TRUE"] <- 0

# Check for significant differences in genome size or GC content
t.test(genome.size ~ Type, mag.sample_data)
t.test(gc.content ~ Type, mag.sample_data)

#### CONSTRAINED AND UNCONSTRAINED ORDINATIONS ####
# Unconstrained (PCA)
mag.pca <- mag.cog_counts
mag.pca$MAG03 <- NULL
mag.pca <- subset(mag.pca, rownames(mag.pca) !="#N/A")

# Convert to COG relative abundances
for(i in 1:ncol(mag.pca)){
  mag.pca[,i] <- 100*mag.pca[,i]/sum(mag.pca[,i])
}

# Calculate the ordination and save axis scores.
mag.pca.calc <- vegan::rda(as.data.frame(t(mag.pca)))
mag.pca.dist <- vegdist(as.data.frame(t(mag.pca)), method="euclidean")

mag.pca.data <- merge(mag.sample_data,
                      vegan::scores(mag.pca.calc, choices=c(1:3), display="sites"),
                      by=0, all=TRUE)
mag.pca.data <- subset(mag.pca.data, MAG !="MAG03")

# Try a constrained ordination that accounts for taxonomic differences and save scores.
set.seed(123987)
mag.dbrda <- cca(t(mag.pca) ~ Phylum + Type,
                 data=mag.pca.data)

mag.dbrda.data <- merge(mag.sample_data,
                        vegan::scores(mag.dbrda, choices=c(1:3), display="sites"),
                        by=0, all=TRUE)
mag.dbrda.data <- subset(mag.dbrda.data, MAG !="MAG03")

summary(mag.dbrda)

# Check the significance of these groupings by either type or phylum
# Use PERMANOVA for the unconstrained ordination and ANOVA for the constrained ordination
adonis(mag.pca.dist ~ Phylum + Type, subset(mag.sample_data, MAG !="MAG03"))

anova(mag.dbrda, by="term")
RsquareAdj(mag.dbrda)

anova(mag.dbrda, cca(t(mag.pca) ~ Phylum,
                     data=mag.pca.data))

# Use the data to plot (see the figures code) before deleting it here.
rm(mag.pca.data, mag.pca.calc, mag.pca.dist, mag.pca, mag.dbrda, mag.dbrda.data,
   mag.ords)

#### LOGISTIC REGRESSION MODELS COMPARING SPONGE- AND WATER-ASSOCIATED COGs ####
# a) Run the initial models ####
# Prepare model data and remove COGs not present in the sample.
mag.model.data <- mag.cog_counts
mag.model.data$MAG03 <- NULL
mag.model.data <- subset(mag.model.data, rownames(mag.model.data) !="#N/A")

mag.model.data$SUM <- rowSums(mag.model.data)
mag.model.data <- subset(mag.model.data, SUM > 0)
mag.model.data$SUM <- NULL

# Calculate relative abundance of each COG in the genome.
for(i in 1:ncol(mag.model.data)){
  mag.model.data[,i] <- 100*mag.model.data[,i]/sum(mag.model.data[,i])
}

# Calculate the prevalence of each COG in each group (i.e., number of sponge-associated and water-associated MAGs in which it appears)
mag.model.data.temp <- mag.model.data
mag.model.data.temp[mag.model.data.temp !=0] <- 1
mag.model.data.temp$SUM <- rowSums(mag.model.data.temp)
mag.model.data.temp$SPONGE <- rowSums(mag.model.data.temp[,c("MAG04", "MAG05", "MAG06",
                                                             "MAG07", "MAG08", "MAG09", "MAG10",
                                                             "MAG12", "MAG18", "MAG24", "MAG25")])
mag.model.data.temp$WATER <- rowSums(mag.model.data.temp[,c("MAG01", "MAG02", "MAG11",
                                                            "MAG13", "MAG14", "MAG15", "MAG16",
                                                            "MAG17", "MAG19", "MAG20", "MAG21", "MAG22", "MAG23")])

# Only model COGs present in at least half of the MAGs in each group...
cogs.for.model <- subset(mag.model.data.temp, WATER > 4 | SPONGE > 4)
mag.model.data <- subset(mag.model.data, rownames(mag.model.data) %in% rownames(cogs.for.model))
mag.model.data <- as.data.frame(t(mag.model.data))
rm(cogs.for.model)

# Add information on the MAGs for use in the model
mag.model.data <- merge(mag.sample_data[,c("Type","Phylum")],
                        mag.model.data, 
                        by=0, all=FALSE)
rownames(mag.model.data) <- mag.model.data$Row.names
mag.model.data$Row.names <- NULL

mag.model.coefficients <- data.frame(matrix(nrow=(ncol(mag.model.data)-2),
                                            ncol=6))
colnames(mag.model.coefficients) <- c("COG", "Intercept", "PhylumBacteroidetes", "PhylumProteobacteria",
                                      "PhylumVerrucomicrobia","Coeff")

# Run the models for every MAG - predict COG association as a function of phylum and COG abundance
for(i in c(3:(ncol(mag.model.data)))){
  model <- glm(Type ~ Phylum + mag.model.data[,i], 
               family="binomial"(link="logit"),
               mag.model.data)
  vector <- c(colnames(mag.model.data)[i],
              model[["coefficients"]])
  mag.model.coefficients[i,] <- vector
}

# Clean coefficients data frame
mag.model.coefficients$COG <- factor(mag.model.coefficients$COG)
for(i in c(2:6)){
  mag.model.coefficients[,i] <- as.numeric(as.character(mag.model.coefficients[,i]))
}

# Reverse coefficients for convenience, so that positive coefficients = sponge-associated
mag.model.coefficients$Coeff <- -1*mag.model.coefficients$Coeff

# Add COG annotations to the logistic regression coefficients data frame
mag.model.coefficients <- merge(mag.model.coefficients,
                                metagenome.cog_annotations,
                                by="COG", all=FALSE)
rownames(mag.model.coefficients) <- mag.model.coefficients$COG
mag.model.coefficients <- merge(mag.model.coefficients,
                                mag.model.data.temp[,c("SPONGE","WATER")],
                                by=0, all=TRUE)
rownames(mag.model.coefficients) <- mag.model.coefficients$Row.names
mag.model.coefficients$Row.names <- NULL

# Calculate the mean abundance of each COG in sponge- and water-associated MAGs
mag.model.data[,c("Type", "Phylum")] <- NULL
mag.model.data <- as.data.frame(t(mag.model.data))
mag.model.data$mean.sponge <- rowMeans(mag.model.data[,c("MAG04", "MAG05", "MAG06",
                                                         "MAG07", "MAG08", "MAG09", "MAG10",
                                                         "MAG12", "MAG18", "MAG24", "MAG25")])
mag.model.data$mean.water <- rowMeans(mag.model.data[,c("MAG01", "MAG02", "MAG11",
                                                        "MAG13", "MAG14", "MAG15", "MAG16",
                                                        "MAG17", "MAG19", "MAG20", "MAG21", "MAG22", "MAG23")])
mag.model.coefficients <- merge(mag.model.coefficients,
                                mag.model.data[,c("mean.sponge","mean.water")],
                                by=0, all=TRUE)
rownames(mag.model.coefficients) <- mag.model.coefficients$Row.names
mag.model.coefficients$Row.names <- NULL

# Remove coefficients related to phylum
mag.model.coefficients <- subset(mag.model.coefficients, is.na(Intercept)=="FALSE")
mag.model.coefficients$Intercept <- NULL
mag.model.coefficients$PhylumBacteroidetes <- NULL
mag.model.coefficients$PhylumProteobacteria <- NULL
mag.model.coefficients$PhylumVerrucomicrobia <- NULL

mag.model.coefficients <- dplyr::distinct(mag.model.coefficients)

mag.model.coefficients$logFC.COG <- log(mag.model.coefficients$mean.sponge / mag.model.coefficients$mean.water)

ggplot(mag.model.coefficients, aes(x=logFC.COG, y=Coeff)) + 
  geom_point() +
  scale_y_continuous(limits=c(-250,250))

mag.model.coefficients <- subset(mag.model.coefficients, WATER > 0 & SPONGE > 0)

# b) Test for class overrepresentation among the top 100 COGs ####
# Count total and sponge-associated COGs.
mag.class_sums <- cog_counter(mag.model.coefficients)
ggplot(mag.model.coefficients, aes(x=Coeff)) + geom_histogram() + scale_x_continuous(limits=c(-50,50))

mag.model.coefficients <- mag.model.coefficients[order(-mag.model.coefficients$Coeff), ]
rownames(mag.model.coefficients) <- seq(1:nrow(mag.model.coefficients))
top <- mag.model.coefficients[51,"Coeff"]
bottom <- mag.model.coefficients[(nrow(mag.model.coefficients)-50), "Coeff"]

# Subset to ALL differentially abundant COGs
tempA <- subset(mag.model.coefficients, Coeff > top | Coeff < bottom)
tempA <- cog_counter(tempA)
colnames(tempA) <- c("Class", "All")

# Subset to all COGs more abundant in sponges.
tempB <- subset(mag.model.coefficients, Coeff > top)
tempB <- cog_counter(tempB)
colnames(tempB) <- c("Class", "Sponge")

# Subset to all COGs more abundant in water.
tempC <- subset(mag.model.coefficients, Coeff < bottom)
tempC <- cog_counter(tempC)
colnames(tempC) <- c("Class", "Water")

# Create a summary of these data
mag.class_sums <- merge(mag.class_sums, tempA, by="Class", all=TRUE)
mag.class_sums <- merge(mag.class_sums, tempB, by="Class", all=TRUE)
mag.class_sums <- merge(mag.class_sums, tempC, by="Class", all=TRUE)
rm(tempA, tempB, tempC)

# Run the Fisher exact test... 
mag.class_fisher_test <- data.frame(matrix(nrow=24, ncol=4))
mag.class_fisher_test[,1] <- mag.class_sums$Class

mag.class_sums[is.na(mag.class_sums)] <- 0

totalM <- nrow(mag.model.coefficients) # Total COGs used in the model

for(p in 1:nrow(mag.class_sums)){
  total_sp <- mag.class_sums[p,"Total_Count"] # Total number of genes within that COG
  
  totalDE <- sum(mag.class_sums$All)
  total_sp_DE <- mag.class_sums[p,"All"] # Total number of differentially expressed genes within that COG
  
  total_sp_notDE <- total_sp - total_sp_DE
  total_DE_notsp <- totalDE - total_sp_DE
  total_notsp_notDE <- totalM - total_sp - total_DE_notsp
  
  mag.class_fisher_test[p,2] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), 
                                                   nrow=2), alternative = "greater")$p.value
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
}

# Second, for over-representation in the pool of genes more abundant in sponges
for(p in 1:nrow(mag.class_sums)){
  total_sp <- mag.class_sums[p,"Total_Count"] # Total number of genes within that COG
  
  totalDE <- sum(mag.class_sums$Sponge)
  total_sp_DE <- mag.class_sums[p,"Sponge"] # Total number of differentially expressed genes within that COG
  
  total_sp_notDE <- total_sp - total_sp_DE
  total_DE_notsp <- totalDE - total_sp_DE
  total_notsp_notDE <- totalM - total_sp - total_DE_notsp
  
  mag.class_fisher_test[p,3] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), nrow=2), alternative = "greater")$p.value
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
}

# Third, for over-representation in the pool of genes more abundant in water
for(p in 1:nrow(mag.class_sums)){
  total_sp <- mag.class_sums[p,"Total_Count"] # Total number of genes within that COG
  
  totalDE <- sum(mag.class_sums$Water)
  total_sp_DE <- mag.class_sums[p,"Water"] # Total number of differentially expressed genes within that COG
  
  total_sp_notDE <- total_sp - total_sp_DE
  total_DE_notsp <- totalDE - total_sp_DE
  total_notsp_notDE <- totalM - total_sp - total_DE_notsp
  
  mag.class_fisher_test[p,4] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), nrow=2), alternative = "greater")$p.value
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
}

colnames(mag.class_fisher_test) <- c("Class", "All", "Sponge", "Water")

rm(totalDE, totalM)

# c) Run a new ordination based on these top 100 COGs ####
# Unconstrained ordination (PCA)
tempA <- subset(mag.model.coefficients, Coeff > top | Coeff < bottom)
tempA$COG <- factor(tempA$COG)

mag.recluster <- mag.cog_counts
for(i in c(1:ncol(mag.recluster))){
  mag.recluster[,i] <- 100*mag.recluster[,i]/sum(mag.recluster[,i])
}
mag.recluster <- subset(mag.recluster, rownames(mag.recluster) %in% tempA$COG)
mag.recluster$MAG03 <- NULL
mag.recluster <- as.data.frame(t(mag.recluster))

mag.recluster.pca <- rda(mag.recluster ~ 1)

mag.recluster.data <- merge(mag.sample_data,
                            vegan::scores(mag.recluster.pca, display="sites", choices=c(1:2)),
                            by=0, all=TRUE)

mag.recluster.data <- subset(mag.recluster.data, MAG !="MAG03")
mag.recluster.data$Phylum <- factor(mag.recluster.data$Phylum)

summary(mag.recluster.pca)

# Constrained ordination (CCA)
set.seed(123987)
mag.recluster.dbrda <- cca(mag.recluster ~ Phylum + Type,
                           data=mag.recluster.data)

mag.recluster.dbrda.data <- merge(mag.sample_data,
                                  vegan::scores(mag.recluster.dbrda, choices=c(1:3), display="sites"),
                                  by=0, all=TRUE)
mag.recluster.dbrda.data <- subset(mag.recluster.dbrda.data, MAG !="MAG03")

summary(mag.recluster.dbrda)

# Check for significant differences
adonis(vegdist(mag.recluster, method="euclidean") ~ Phylum + Type, subset(mag.sample_data, MAG !="MAG03"))

anova(mag.recluster.dbrda, by="term")
RsquareAdj(mag.recluster.dbrda)

# Use these data for plots, then remove
rm(mag.recluster, mag.recluster.data, mag.recluster.dbrda.data, mag.recluster.pca, tempA)

#### IMPORT REFERENCE GENOME DATA ####
# Import the COG profiles for each reference genome (i.e., number of each COG in each genome)
reference.cog_counts_temp <- list(
  ref04 = read.table("raw_data/reference_mag_cog_profiles/cog_MAG04_pair_GCA_903891375.1.txt", sep="\t", quote=""),
  ref06 = read.table("raw_data/reference_mag_cog_profiles/cog_MAG06_pair_GCA_903886275.1.txt", sep="\t", quote=""),
  ref07 = read.table("raw_data/reference_mag_cog_profiles/cog_MAG07_pair_GCA_903845935.1.txt", sep="\t", quote=""),
  ref09 = read.table("raw_data/reference_mag_cog_profiles/cog_MAG09_pair_GCA_903931835.1.txt", sep="\t", quote=""),
  ref10 = read.table("raw_data/reference_mag_cog_profiles/cog_MAG10_pair_GCF_002251735.1.txt", sep="\t", quote=""),
  ref12 = read.table("raw_data/reference_mag_cog_profiles/cog_MAG12_pair_GCA_011046725.1.txt", sep="\t", quote=""),
  ref18 = read.table("raw_data/reference_mag_cog_profiles/cog_MAG18_pair_GCA_013140395.1.txt", sep="\t", quote=""),
  ref24 = read.table("raw_data/reference_mag_cog_profiles/cog_MAG24_pair_GCF_009939225.1.txt", sep="\t", quote=""),
  ref25 = read.table("raw_data/reference_mag_cog_profiles/cog_MAG25_pair_GCF_900167075.1.txt", sep="\t", quote="")
)

for(i in c(1:length(reference.cog_counts_temp))){
  colnames(reference.cog_counts_temp[[i]]) <- c("COG", "Gene", names(reference.cog_counts_temp)[i],
                                                "Score1", "EValue", "Name")
  rownames(reference.cog_counts_temp[[i]]) <- reference.cog_counts_temp[[i]]$COG
  reference.cog_counts_temp[[i]][,c("COG","Gene","Score1","EValue","Name")] <- NULL
}

reference.cog_counts <- transform(merge(reference.cog_counts_temp[[1]],
                                        reference.cog_counts_temp[[2]],
                                        by=0, all=TRUE), row.names=Row.names, Row.names=NULL)

for(i in c(3:length(reference.cog_counts_temp))){
  reference.cog_counts <-  transform(merge(reference.cog_counts,
                                           reference.cog_counts_temp[[i]],
                                           by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
}

rm(reference.cog_counts_temp)
reference.cog_counts[is.na(reference.cog_counts)=="TRUE"] <- 0

reference.cog_counts <-  transform(merge(metagenome.cog_annotations[,c("Class","Full.Name")],
                                         reference.cog_counts,
                                         by=0, all.x=FALSE, all.y=TRUE), row.names=Row.names, Row.names=NULL)


#### FIND COGS THAT APPEAR IN ONLY ONE COPY IN ALL MAGS AND REFERENCES ####
# Identify COGs present with one copy in all MAGs and reference genomes
mag.singles <- merge(mag.cog_counts, reference.cog_counts, by=0, all.x=TRUE, all.y=TRUE)
rownames(mag.singles) <- mag.singles$Row.names
mag.singles$Row.names <- NULL

mag.singles <- subset(mag.singles,
                      ref04==1 & ref06==1 & ref07==1 & ref09==1 & ref10==1 & ref12==1 & ref18==1 & ref24==1 & ref25==1 &
                        MAG04==1 & MAG06==1 & MAG07==1 & MAG09==1 & MAG10==1 & MAG12==1 & MAG18==1 & MAG24==1 & MAG25==1)

mag.singles <- subset(mag.singles, rownames(mag.singles) != "COG0361")
mag.singles <- rownames(mag.singles)

# Import a list showing all the genes in each MAG and the COGs they encode
mag.cogs <- list(
  mag04 = read.table("raw_data/mag_webmga_data/mag04/cog-raw.txt", sep="\t", quote=""),
  mag06 = read.table("raw_data/mag_webmga_data/mag06/cog-raw.txt", sep="\t", quote=""),
  mag07 = read.table("raw_data/mag_webmga_data/mag07/cog-raw.txt", sep="\t", quote=""),
  mag09 = read.table("raw_data/mag_webmga_data/mag09/cog-raw.txt", sep="\t", quote=""),
  mag10 = read.table("raw_data/mag_webmga_data/mag10/cog-raw.txt", sep="\t", quote=""),
  mag12 = read.table("raw_data/mag_webmga_data/mag12/cog-raw.txt", sep="\t", quote=""),
  mag18 = read.table("raw_data/mag_webmga_data/mag18/cog-raw.txt", sep="\t", quote=""),
  mag24 = read.table("raw_data/mag_webmga_data/mag24/cog-raw.txt", sep="\t", quote=""),
  mag25 = read.table("raw_data/mag_webmga_data/mag25/cog-raw.txt", sep="\t", quote=""),
  ref04 = read.table("raw_data/mag_webmga_data/MAG04_pair_GCA_903891375.1/cog-raw.txt", sep="\t", quote=""),
  ref06 = read.table("raw_data/mag_webmga_data/MAG06_pair_GCA_903886275.1/cog-raw.txt", sep="\t", quote=""),
  ref07 = read.table("raw_data/mag_webmga_data/MAG07_pair_GCA_903845935.1/cog-raw.txt", sep="\t", quote=""),
  ref09 = read.table("raw_data/mag_webmga_data/MAG09_pair_GCA_903931835.1/cog-raw.txt", sep="\t", quote=""),
  ref10 = read.table("raw_data/mag_webmga_data/ref10/cog-raw.txt", sep="\t", quote=""),
  ref12 = read.table("raw_data/mag_webmga_data/MAG12_pair_GCA_011046725.1/cog-raw.txt", sep="\t", quote=""),
  ref18 = read.table("raw_data/mag_webmga_data/MAG18_pair_GCA_013140395.1/cog-raw.txt", sep="\t", quote=""),
  ref24 = read.table("raw_data/mag_webmga_data/MAG24_pair_GCF_009939225.1/cog-raw.txt", sep="\t", quote=""),
  ref25 = read.table("raw_data/mag_webmga_data/MAG25_pair_GCF_900167075.1/cog-raw.txt", sep="\t", quote="")
)

# Subset COGs to just the single-copy COGs
for(i in 1:length(mag.cogs)){
  mag.cogs[[i]] <- subset(mag.cogs[[i]], V13 %in% mag.singles)
  mag.cogs[[i]] <- mag.cogs[[i]][order(mag.cogs[[i]]$V13), ]
  mag.cogs[[i]]$V1 <- factor(mag.cogs[[i]]$V1)
  mag.cogs[[i]]$V1 <- forcats::fct_inorder(mag.cogs[[i]]$V1)
}

# Import the sequences from each MAG
mag.faas <- list(
  mag04 = seqinr::read.fasta("raw_data/mag_protein_sequences/sponges.bins.28.faa",
                             seqtype="AA"),
  mag06 = seqinr::read.fasta("raw_data/mag_protein_sequences/70_sub.faa",
                             seqtype="AA"),
  mag07 = seqinr::read.fasta("raw_data/mag_protein_sequences/MaxBin2_out.006_sub.faa",
                             seqtype="AA"),
  mag09 = seqinr::read.fasta("raw_data/mag_protein_sequences/MaxBin2_out.111.faa",
                             seqtype="AA"),
  mag10 = seqinr::read.fasta("raw_data/mag_protein_sequences/MaxBin2_out.014.faa",
                             seqtype="AA"),
  mag12 = seqinr::read.fasta("raw_data/mag_protein_sequences/48_sub.faa",
                             seqtype="AA"),
  mag18 = seqinr::read.fasta("raw_data/mag_protein_sequences/50.faa",
                             seqtype="AA"),
  mag24 = seqinr::read.fasta("raw_data/mag_protein_sequences/117.faa",
                             seqtype="AA"),
  mag25 = seqinr::read.fasta("raw_data/mag_protein_sequences/131.faa",
                             seqtype="AA"),
  
  ref04 = seqinr::read.fasta("raw_data/ref_protein_sequences/pair_MAG04_GCA_903891375.1_proteins.faa",
                             seqtype="AA"),
  ref06 = seqinr::read.fasta("raw_data/ref_protein_sequences/pair_MAG06_GCA_903886275.1_proteins.faa",
                             seqtype="AA"),
  ref07 = seqinr::read.fasta("raw_data/ref_protein_sequences/pair_MAG07_GCA_903845935.1_proteins.faa",
                             seqtype="AA"),
  ref09 = seqinr::read.fasta("raw_data/ref_protein_sequences/pair_MAG09-10_GCA_903931835.1_proteins.faa",
                             seqtype="AA"),
  ref10 = seqinr::read.fasta("raw_data/ref_protein_sequences/new_MAG10_counterpart.faa",
                             seqtype="AA"),
  ref12 = seqinr::read.fasta("raw_data/ref_protein_sequences/pair_MAG12_GCA_011046725.1_protein.faa",
                             seqtype="AA"),
  ref18 = seqinr::read.fasta("raw_data/ref_protein_sequences/pair_MAG18_GCA_013140395.1_protein.faa",
                             seqtype="AA"),
  ref24 = seqinr::read.fasta("raw_data/ref_protein_sequences/ref_MAG24-25_GCF_009939225.1_protein.faa",
                             seqtype="AA"),
  ref25 = seqinr::read.fasta("raw_data/ref_protein_sequences/ref_MAG25_GCF_900167075.1_protein.faa",
                             seqtype="AA")
)

# Subset sequences to only the sequences for single-copy COGs
for(i in c(1:length(mag.faas))){
  mag.faas[[i]] <- subset(mag.faas[[i]], names(mag.faas[[i]]) %in% mag.cogs[[i]]$V1)
}

# Order the COGs in the same order in every file
mag.faas.order <- list()
for(j in c(1:length(mag.faas))){
  temp <- list()
  for(i in c(1:12)){
    temp[[i]] <- mag.faas[[j]][[levels(mag.cogs[[j]]$V1)[i]]]
  }
  names(temp) <- rep(names(mag.faas)[j], 12)
  #names(temp) <- paste0(names(mag.faas)[j]), "_", mag.singles)
  mag.faas.order[[j]] <- temp
  rm(temp)
}

names(mag.faas.order) <- names(mag.faas)

# Write the FASTA file with the concatenated sequence for all single-copy genes
seqinr::write.fasta(mag.faas.order[[1]], names=names(mag.faas.order[[1]]), file.out="mag_singles/mag04.singles.faa")
seqinr::write.fasta(mag.faas.order[[2]], names=names(mag.faas.order[[2]]), file.out="mag_singles/mag06.singles.faa")
seqinr::write.fasta(mag.faas.order[[3]], names=names(mag.faas.order[[3]]), file.out="mag_singles/mag07.singles.faa")
seqinr::write.fasta(mag.faas.order[[4]], names=names(mag.faas.order[[4]]), file.out="mag_singles/mag09.singles.faa")
seqinr::write.fasta(mag.faas.order[[5]], names=names(mag.faas.order[[5]]), file.out="mag_singles/mag10.singles.faa")
seqinr::write.fasta(mag.faas.order[[6]], names=names(mag.faas.order[[6]]), file.out="mag_singles/mag12.singles.faa")
seqinr::write.fasta(mag.faas.order[[7]], names=names(mag.faas.order[[7]]), file.out="mag_singles/mag18.singles.faa")
seqinr::write.fasta(mag.faas.order[[8]], names=names(mag.faas.order[[8]]), file.out="mag_singles/mag24.singles.faa")
seqinr::write.fasta(mag.faas.order[[9]], names=names(mag.faas.order[[9]]), file.out="mag_singles/mag25.singles.faa")

seqinr::write.fasta(mag.faas.order[[10]], names=names(mag.faas.order[[10]]), file.out="mag_singles/ref04.singles.faa")
seqinr::write.fasta(mag.faas.order[[11]], names=names(mag.faas.order[[11]]), file.out="mag_singles/ref06.singles.faa")
seqinr::write.fasta(mag.faas.order[[12]], names=names(mag.faas.order[[12]]), file.out="mag_singles/ref07.singles.faa")
seqinr::write.fasta(mag.faas.order[[13]], names=names(mag.faas.order[[13]]), file.out="mag_singles/ref09.singles.faa")
seqinr::write.fasta(mag.faas.order[[14]], names=names(mag.faas.order[[14]]), file.out="mag_singles/ref10.singles.faa")
seqinr::write.fasta(mag.faas.order[[15]], names=names(mag.faas.order[[15]]), file.out="mag_singles/ref12.singles.faa")
seqinr::write.fasta(mag.faas.order[[16]], names=names(mag.faas.order[[16]]), file.out="mag_singles/ref18.singles.faa")
seqinr::write.fasta(mag.faas.order[[17]], names=names(mag.faas.order[[17]]), file.out="mag_singles/ref24.singles.faa")
seqinr::write.fasta(mag.faas.order[[18]], names=names(mag.faas.order[[18]]), file.out="mag_singles/ref25.singles.faa")

# Run these sequenes in MAFFT to obtain an alignment

rm(mag.faas.order, mag.faas)

#### SPONGE- vs. WATER-ASSOCIATED MAG vs. REFERENCE: COG PRESENCE/ABSENCE AND ABUNDANCE ####
# Obtain a data frame showing the number of sponge-associated and water-associated MAGs that contain each COG
all.mag.cogs <- mag.model.data.temp
all.mag.cogs[,c(1:25)] <- NULL
all.mag.cogs[is.na(all.mag.cogs)=="TRUE"] <- 0

# Obtain a data frame showing the number of reference genomes that contain each COG
ref.mag.cogs <- reference.cog_counts
ref.mag.cogs[,c(3:ncol(ref.mag.cogs))][ref.mag.cogs[,c(3:ncol(ref.mag.cogs))] > 0] <- 1
ref.mag.cogs$REF <- rowSums(ref.mag.cogs[,c(3:ncol(ref.mag.cogs))])
ref.mag.cogs[,c(1:(ncol(ref.mag.cogs)-1))] <- NULL
ref.mag.cogs[is.na(ref.mag.cogs)=="TRUE"] <- 0

# Combine into a single data frame
all.mag.cogs <- merge(all.mag.cogs, ref.mag.cogs, by=0, all.x=TRUE, all.y=TRUE)
rownames(all.mag.cogs) <- all.mag.cogs$Row.names
all.mag.cogs$Row.names <- NULL
all.mag.cogs[,c(1:3)][is.na(all.mag.cogs[,c(1:3)])=="TRUE"] <- 0

all.mag.cogs <- merge(all.mag.cogs,
                      metagenome.cog_annotations[,c("Class","Full.Name")],
                      by=0, all.x=TRUE, all.y=FALSE)

rm(ref.mag.cogs)

# a) MAGs unique to each group (i.e., present in one group but not others ####
mag.Venn.diagrams <- list(
  uniq.to.sp.vs.wat = subset(all.mag.cogs, SPONGE > 0 & WATER == 0),
  uniq.to.wat.vs.sp = subset(all.mag.cogs, SPONGE == 0 & WATER > 0),
  uniq.to.sp.vs.ref = subset(all.mag.cogs, SPONGE > 0 & REF == 0)
)

# In sponges but not water
nrow(mag.Venn.diagrams$uniq.to.sp.vs.wat) # 383
nrow(subset(mag.Venn.diagrams$uniq.to.sp.vs.wat, SPONGE >= 4)) # 25 of these present in at least 4 samples

# In water but not sponges
nrow(mag.Venn.diagrams$uniq.to.wat.vs.sp) # 515
nrow(subset(mag.Venn.diagrams$uniq.to.wat.vs.sp, WATER >= 4)) # 30 of these present in at least 4 samples

# In sponges but not references
nrow(mag.Venn.diagrams$uniq.to.sp.vs.ref) # 178 - sponges have all the COGs present in the references
nrow(subset(mag.Venn.diagrams$uniq.to.sp.vs.ref, SPONGE >= 4)) # 3 of these present in at least 4 samples

# Explore the COG classes represented by these "unique" MAGs
sp.mag.cog.counted <- cog_counter(mag.Venn.diagrams$uniq.to.sp.vs.wat)
colnames(sp.mag.cog.counted) <- c("Class", "Uniq.to.Sponge.v.Water")
wat.mag.cog.counted <- cog_counter(mag.Venn.diagrams$uniq.to.wat.vs.sp)
colnames(wat.mag.cog.counted) <- c("Class", "Uniq.to.Water")
sp.mag.cog.counted.ref <- cog_counter(mag.Venn.diagrams$uniq.to.sp.vs.ref)
colnames(sp.mag.cog.counted.ref) <- c("Class", "Uniq.to.Sponge.v.Ref")

count <- cog_counter(subset(all.mag.cogs, SPONGE > 0 | WATER > 0))
count <- merge(count, sp.mag.cog.counted, by="Class", all=TRUE)
count <- merge(count, wat.mag.cog.counted, by="Class", all=TRUE)
count <- merge(count, sp.mag.cog.counted.ref, by="Class", all=TRUE)

count <- subset(count, is.na(Uniq.to.Sponge.v.Water)=="FALSE")
count[is.na(count)=="TRUE"] <- 0

rm(sp.mag.cog.counted, wat.mag.cog.counted, sp.mag.cog.counted.ref)

# Test for overrepresentation among COGs unique to sponges
totalM <- 3076 # Total number of available COGs

mag.uniq.to.sponge.fisher.test <- data.frame(matrix(nrow=23, ncol=3))

# vs the water-associated MAGs
for(p in 1:nrow(count)){
  total_sp <- count[p,"Total_Count"] # Total number of genes within that COG
  
  totalDE <- sum(count$Uniq.to.Sponge.v.Water)
  total_sp_DE <- count[p,"Uniq.to.Sponge.v.Water"] # Total number of differentially expressed genes within that COG
  
  total_sp_notDE <- total_sp - total_sp_DE
  total_DE_notsp <- totalDE - total_sp_DE
  total_notsp_notDE <- totalM - total_sp - total_DE_notsp
  
  mag.uniq.to.sponge.fisher.test[p,2] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), 
                                                            nrow=2), alternative = "greater")$p.value
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
}
mag.uniq.to.sponge.fisher.test$X1 <- count$Class

# vs. the reference genomes
for(p in 1:nrow(count)){
  total_sp <- count[p,"Total_Count"] # Total number of genes within that COG
  
  totalDE <- sum(count$Uniq.to.Sponge.v.Water)
  total_sp_DE <- count[p,"Uniq.to.Sponge.v.Ref"] # Total number of differentially expressed genes within that COG
  
  total_sp_notDE <- total_sp - total_sp_DE
  total_DE_notsp <- totalDE - total_sp_DE
  total_notsp_notDE <- totalM - total_sp - total_DE_notsp
  
  mag.uniq.to.sponge.fisher.test[p,3] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), 
                                                            nrow=2), alternative = "greater")$p.value
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
}

colnames(mag.uniq.to.sponge.fisher.test) <- c("Class", "Sponge.v.Water", "Sponge.v.Ref")

write.csv(mag.uniq.to.sponge.fisher.test, "mag.uniq.to.sponge.fisher.test.csv")
write.csv(subset(mag.Venn.diagrams$uniq.to.sp.vs.wat, SPONGE >= 4),
          "mag.uniq.to.sponge.COG.list.csv")

rm(count)

# b) Identify COGs present in all sponge-associated MAGs (i.e., core genome) ####
mag.core.cogs <- mag.model.data.temp

nrow(subset(mag.core.cogs, SPONGE==11)) # 134 core sponge genes
nrow(subset(mag.core.cogs, SPONGE > 0)) # Out of 2557 COGs present in sponges
134 / 2557 # 5.24% of sponge COGs present in the core

nrow(subset(mag.core.cogs, WATER==13)) # 76 core water cogs
nrow(subset(mag.core.cogs, WATER > 0)) # Out of 2689 COGs present in water
76 / 2689 # 2.8% of water COGs present in the core

nrow(subset(mag.core.cogs, SPONGE==11 & WATER==13)) # 27 universal genes

mag.core.cogs[,c(1:24)] <- NULL
mag.core.cogs <- merge(metagenome.cog_annotations[,c("COG","Class","Full.Name")],
                       mag.core.cogs, by=0, all.x=FALSE, all.y=TRUE)

mag.core.cogs <- subset(mag.core.cogs, SPONGE==11)

count <- cog_counter(mag.core.cogs)
colnames(count) <- c("Class", "Core_Count")

count2 <- cog_counter(subset(metagenome.cog_annotations, COG %in% rownames(mag.model.data.temp)))

count <- merge(count, count2, by="Class", all=TRUE)

count[is.na(count)=="TRUE"] <- 0

# Test for overrepresentation in the core COGs
totalM <- 3072 # Total number of available COGs

mag.core.sponge.fisher.test <- data.frame(matrix(nrow=25, ncol=2))

for(p in 1:nrow(count)){
  total_sp <- count[p,"Total_Count"] # Total number of genes within that COG
  
  totalDE <- 134
  total_sp_DE <- count[p,"Core_Count"] # Total number of differentially expressed genes within that COG
  
  total_sp_notDE <- total_sp - total_sp_DE
  total_DE_notsp <- totalDE - total_sp_DE
  total_notsp_notDE <- totalM - total_sp - total_DE_notsp
  
  mag.core.sponge.fisher.test[p,2] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), 
                                                         nrow=2), alternative = "greater")$p.value
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
}
mag.core.sponge.fisher.test$X1 <- count$Class

write.csv(mag.core.sponge.fisher.test, "mag.core.sponge.fisher.test.csv")
write.csv(mag.core.cogs, "mag.core.cogs.csv")

rm(count, count2, totalDE, totalM)

#### MAG: PAIRED T.TEST BETWEEN SPONGE-ASSOCIATED AND REFERENCE MAGS ####
# Get MAGs and reference genomes into the same data frame
orig <- mag.cog_counts[,c("MAG04", "MAG06", "MAG07", "MAG09", "MAG10","MAG12", "MAG18", "MAG24", "MAG25")]
ref <- reference.cog_counts
ref$Class <- NULL
ref$Full.Name <- NULL

finaldf <- transform(merge(orig, ref, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
finaldf[is.na(finaldf)=="TRUE"] <- 0

for(i in c(1:ncol(finaldf))){
  finaldf[,i] <- 100 * finaldf[,i] / sum(finaldf[,i])
}

finaldf <- as.data.frame(t(finaldf))
finaldf$Group <- c(rep("sponge", 9),
                   rep("reference", 9))

pairedtest <- data.frame(matrix(nrow=(ncol(finaldf)-1),
                                ncol=5))

for(i in c(1:(ncol(finaldf)-1))){
  test <- t.test(finaldf[,i] ~ Group, finaldf, paired=TRUE)
  pairedtest[i,1] <- colnames(finaldf)[i]
  pairedtest[i,2] <- test$statistic
  pairedtest[i,3] <- test$p.value
  pairedtest[i,4] <- mean(subset(finaldf, Group=="sponge")[,i])
  pairedtest[i,5] <- mean(subset(finaldf, Group=="reference")[,i])
}

pairedtest$adjust <- p.adjust(pairedtest$X3, method="BH")
colnames(pairedtest) <- c("COG", "t", "p", "sponge", "reference", "p.adjust")

pairedtest$spvsref <- pairedtest$sponge / pairedtest$reference

pairedtest <- merge(pairedtest, metagenome.cog_annotations[,c("COG", "Class", "Full.Name")],
                    by="COG", all.x=TRUE, all.y=FALSE)

all.mag.cogs$COG <- all.mag.cogs$Row.names

pairedtest <- merge(pairedtest, all.mag.cogs[,c("COG","WATER","SPONGE","REF")],
                    by="COG", all.x=TRUE, all.y=FALSE)
pairedtest <- subset(pairedtest, COG !="#N/A")

mag.paired.ttest <- pairedtest
rm(pairedtest, test, finaldf, orig, ref)

# Merge paired t-test results with COG regression model results
head(mag.paired.ttest)
head(mag.model.coefficients)

p1 <- mag.paired.ttest[,c("COG","p","spvsref","WATER","SPONGE","REF")]
p2 <- mag.model.coefficients[,c("COG","Coeff","mean.sponge","mean.water")]

mag.cog.test.summary <- merge(p1, p2, by="COG", all.x=TRUE, all.y=TRUE)
mag.cog.test.summary <- merge(mag.cog.test.summary, metagenome.cog_annotations[,c("COG","Full.Name")],
            by="COG", all.x=TRUE, all.y=FALSE)
rm(p1, p2)

