##################################################### 16S rRNA AMPLICON ANALYSIS ###########################################
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

sample_data <- read.csv("~/Documents/sponges_final/sponge_metadata_2020-10-8.csv")
rownames(sample_data) <- sample_data$SampleID
sample_data$Type <- factor(sample_data$Type, levels=c("sponge", "water", "biofilm"))
sample_data$Source <- factor(sample_data$Source, levels=c("Sooke", "Nanaimo", "Cowichan"))

#### IDENTIFY AND REMOVE CONTAMINANTS ####
# Store taxa names as a 'refseq' object and rename taxa with temporary names.
pseq.filter <- pseq.original

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

# Remove negative controls from phyloseq object.
pseq.filter <- subset_samples(pseq.filter, Type !="control") # Remove control samples from further analysis.
pseq.filter <- prune_taxa(taxa_sums(pseq.filter) > 0, pseq.filter)

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

# Calculate abundance-weighted OTU overlap for each pair of technical replicates.
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

# Check for differences in alpha diversity.
temp <- cbind(sample_data(ctrl.replicates), estimate_richness(ctrl.replicates))

car::leveneTest(Observed ~ Type, temp) # p = 0.5311
t.test(Observed ~ Type, temp, var.equal=TRUE) # t=0.82688, d=4, p=0.4548

car::leveneTest(Shannon ~ Type, temp) # p = 0.9205
t.test(Shannon ~ Type, temp, var.equal=TRUE) # t=-1.0837, df=4, p=0.3395

rm(temp)

# Check for differences in beta diversity using 1) bdiv.ordinations.all and 2) between- and among-sample distances
ctrl.replicate.data <- list()

for(i in c(1:2)){
  # Do this both using the standard Bray-Curtis approach and the CLR transform approach.
  if(i==1){
    temp.otu <- as.data.frame(otu_table(rarefy_even_depth(ctrl.replicates,
                                                          sample.size=min(sample_sums(ctrl.replicates)),
                                                          rngseed=700,
                                                          replace=FALSE,
                                                          trimOTUs=TRUE)))
    temp.dist <- vegdist(temp.otu, method="bray")
    temp.pca <- capscale(temp.dist ~ 1)
  }
  if(i==2){
    x1 <- transform_sample_counts(pseq.filter, function(x) 100*x/sum(x))
    x2 <- as.data.frame(colMeans(otu_table(x1)))
    x1 <- filter_taxa(x1, function(x) mean(x) > 0.001, TRUE)
    temp.otu <- subset_taxa(ctrl.replicates, taxa_names(ctrl.replicates) %in% taxa_names(x1))
    rm(x1, x2) 
    temp.otu <- as.data.frame(otu_table(microbiome::transform(ctrl.replicates, "clr")))
    temp.dist <- vegdist(temp.otu, method="euclidean")
    temp.pca <- rda(temp.otu)
  }
  
  temp.scores <- cbind(sample_data(ctrl.replicates),
                       vegan::scores(temp.pca, choices=c(1:2), display="sites"))
  colnames(temp.scores)[c(4:5)] <- c("PC1", "PC2")
  
  plot <- ggplot(temp.scores, aes(x=PC1, y=PC2)) +
    geom_point(aes(color=Type)) +
    geom_label(aes(label=SampleID))
  
  # Compare between-pair and among-pair distances
  test <- adonis(temp.dist ~ Type,
                 as.data.frame(cbind(sample_data(ctrl.replicates))), permutations=999)
  
  temp.dist <- as.matrix(temp.dist)
  vector1 <- c(temp.dist["SKE2","K0SKE2"], temp.dist["SKE3","K0SKE3"], temp.dist["SKE5","K0SKE5"]) # Between pairs.
  
  temp.dist[,"K0SKE2"] <- rep(NA)
  temp.dist[,"K0SKE3"] <- rep(NA)
  temp.dist[,"K0SKE5"] <- rep(NA)
  temp.dist[temp.dist==0] <- NA
  
  vector2 <- rowMeans(temp.dist, na.rm=TRUE)[c(1:3)] # Average distances among non-pairs.
  
  temp <- data.frame(Sample=c("SKE2", "SKE3", "SKE5"),
                     Between=vector1,
                     Among=vector2)
  
  ctrl.replicate.data[[i]] <- list()
  ctrl.replicate.data[[i]][[1]] <- plot
  ctrl.replicate.data[[i]][[2]] <- test
  ctrl.replicate.data[[i]][[3]] <- temp
  
  names(ctrl.replicate.data[[i]]) <- c("plot", "adonis", "distance")
  
  rm(temp.otu, temp.dist, temp.pca, temp.scores,
     plot, test, vector1, vector2, temp)
}

# Check for differences in general taxon abundances.
temp.replicate <- tax_glom(ctrl.replicates, taxrank="Genus")

replicate1 <- subset_samples(temp.replicate, sample_names(temp.replicate) %in% c("SKE2", "K0SKE2"))
replicate2 <- subset_samples(temp.replicate, sample_names(temp.replicate) %in% c("SKE3", "K0SKE3"))
replicate3 <- subset_samples(temp.replicate, sample_names(temp.replicate) %in% c("SKE5", "K0SKE5"))

replicate1 <- prune_taxa(taxa_sums(replicate1) > 0, replicate1)
replicate2 <- prune_taxa(taxa_sums(replicate2) > 0, replicate2)
replicate3 <- prune_taxa(taxa_sums(replicate3) > 0, replicate3)

replicate1 <- as.data.frame(t(otu_table(transform_sample_counts(replicate1, function(x) 100*x/sum(x)))))
replicate2 <- as.data.frame(t(otu_table(transform_sample_counts(replicate2, function(x) 100*x/sum(x)))))
replicate3 <- as.data.frame(t(otu_table(transform_sample_counts(replicate3, function(x) 100*x/sum(x)))))

grid.arrange(
  ggplot(replicate1, aes(x=SKE2, y=K0SKE2)) + geom_point() + scale_y_log10() + scale_x_log10(),
  ggplot(replicate2, aes(x=SKE3, y=K0SKE3)) + geom_point() + scale_y_log10() + scale_x_log10(),
  ggplot(replicate3, aes(x=SKE5, y=K0SKE5)) + geom_point() + scale_y_log10() + scale_x_log10(),
  ncol=3)

cor.test(replicate1$SKE2, replicate1$K0SKE2, method = "spearman")
cor.test(replicate2$SKE3, replicate2$K0SKE3, method = "spearman")
cor.test(replicate3$SKE5, replicate3$K0SKE5, method = "spearman")

rm(replicate1, replicate2, replicate3, temp.replicate)

# Remove replicates from analysis and move on.
pseq.filter <- subset_samples(pseq.filter, Type !="replicate")
pseq.filter <- prune_taxa(taxa_sums(pseq.filter) > 0, pseq.filter)

#### DATA PRE-PROCESSING (RAREFACTION, TAXON FILTERING, ETC.) ####
sample_data <- subset(sample_data, Type !="replicate" & Type !="control")

# Examine read count distribution to determine appropriate cutoff threshold.
temp <- as.data.frame(sample_sums(pseq.filter))
ggplot(temp, aes(x=temp[,1])) + geom_histogram() # No need to remove any samples.
rm(temp)

# Calculate read profiles for each sample
df <- as.data.frame(sample_sums(pseq.original))
df <- subset(df, !(rownames(df) %in% c("NCFILT1", "K1", "K2", "K0SKE2", "K0SKE3", "K0SKE5")))
min(df[,1])
max(df[,1])
mean(df[,1])
sd(df[,1])
rm(df)

# Remove ASVs that have a relative abundance < 0.001%
#temp <- transform_sample_counts(pseq.filter, function(x) 100*x/sum(x))
#x2 <- as.data.frame(colMeans(otu_table(temp)))
#temp <- filter_taxa(temp, function(x) mean(x) > 0.001, TRUE)
#pseq.filter <- subset_taxa(pseq.filter, taxa_names(pseq.filter) %in% taxa_names(temp))
#rm(temp, x2)
#taxa_names(pseq.filter) <- paste0("ASV", seq(ntaxa(pseq.filter)))

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

# Compile the data into a single list.
temp_data <- sample_data
temp_data$var <- paste0(temp_data$Type, temp_data$Source)

x <- as.data.frame(as.matrix(TukeyHSD(aov(Observed_Extrap ~ var, temp_data))[[1]]))


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
   adiv.signif.Sooke, adiv.signif.Nanaimo, adiv.signif.Cowichan, adiv.signif.super)

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
ggplot(subset(temp.melt, variable %in% c("Observed", "Observed_Extrap", "Chao1", "Shannon", "Shannon_Extrap", "PD", "PD_Extrap")),
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

# Calculate NMDS bdiv.ordinations.all for the Bray-Curtis, Jaccard, weighted and unweighted UniFrac distances.
for(i in c(1,2,4,5)){
  # Run the NMDS model for each of the non-Aitchison distances.
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
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq)
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
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq)
}

names(diff.abund.sp.location) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# (c) Pairwise: sponge-water ####
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
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq)
}

names(diff.abund.SPWAT) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# (d) Pairwise: sponge-biofilm #####
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
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq)
}

names(diff.abund.SPBIO) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# (e) Pairwise: water-biofilm ####
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
  
  rm(poo1, poo2, poo3, poo, vector.means, final, p, temp.transform, temp.transform.otu, temp.physeq)
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

#### ASV OVERLAP ANALYSIS #########
# Create separate phyloseq objects for (1) sponge, (2) water, (3) biofilm, (4) SKE sponge, (5) NMO sponge, (6) COW sponge
temp.physeq <- pseq.clr
#temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
#temp.physeq <- filter_taxa(temp.physeq, function(x) mean(x) > 0.01, TRUE)
#temp.physeq <- subset_taxa(temp.physeq, taxa_names(temp.physeq) %in% taxa_names(temp))

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

# Calculate PERMANOVA results (Aitchison distance) for each comparison).
temp.list <- list(temp.physeq.sponge, temp.physeq.water, temp.physeq.biofilm,
                  temp.physeq.SKE, temp.physeq.NMO, temp.physeq.COW)
temp.dist <- list()
temp.perm <- data.frame(matrix(ncol=7))
  
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
  temp.perm[i,2] <- temp[["aov.tab"]]$F.Model[1]
  temp.perm[i,3] <- temp[["aov.tab"]]$Df[1]
  temp.perm[i,4] <- temp[["aov.tab"]]$R2[1]
  temp.perm[i,5] <- temp[["aov.tab"]]$`Pr(>F)`[1]
  temp.perm[i,6] <- tempd$tab[1,4]
  temp.perm[i,7] <- tempd$tab[1,6]
  rm(temp, tempd)
}
temp.perm[,1] <- c("sponge", "water", "biofilm", "SKE", "NMO", "COW")
colnames(temp.perm) <- c("ordination", "F", "Df", "R2", "p", "Disp_F", "Disp_p")

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

river.specific <- rbind(r1, r2, r3)

river.specific <- merge(river.specific,
                        as.data.frame(cbind(tax_table(pseq.clr.sub))),
                        by=0, all=FALSE)

river.specific$water <- (river.specific$Row.names %in% water)
river.specific$biofilm <- (river.specific$Row.names %in% biofilm)
  
rm(sponge, water, biofilm, SKE_sponge, NMO_sponge, COW_sponge,
   temp.physeq.sponge, temp.physeq.water, temp.physeq.biofilm,
   temp.physeq.sponge.SKE, temp.physeq.sponge.NMO, temp.physeq.sponge.COW,
   temp.physeq.COW, temp.physeq.NMO, temp.physeq.SKE,
   temp.physeq, temp.otu)

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

##################################################### METAGENOME ANALYSIS I: COGS ###########################################
#### EXPLORE TAXONOMIC INFORMATION IN RELATION TO AMPLICONS ####
# Import taxonomic data. #
metagenome.otu <- read.csv("~/Downloads/sponges_metagenomes/otu_counts_table.csv", row.names="OTU_ID")
metagenome.taxa <- read.csv("~/Downloads/sponges_metagenomes/taxonomy_table.csv", row.names="OTU_ID")
metagenome.data <- read.csv("~/Downloads/sponges_metagenomes/sample_metadata.csv")
rownames(metagenome.data) <- metagenome.data$Sample

metagenome.otu <- otu_table(metagenome.otu, taxa_are_rows=TRUE)
metagenome.taxa <- tax_table(as.matrix(metagenome.taxa))

metagenome.pseq <- phyloseq(otu_table(metagenome.otu, taxa_are_rows=TRUE),
                            sample_data(metagenome.data),
                            tax_table(metagenome.taxa))

# Check kingdom-level abundances to see how much archaea, viruses, and viroids contribute to the sample
otu_table(transform_sample_counts(tax_glom(metagenome.pseq, taxrank="Kingdom"), 
                                  function(x) 100*x/sum(x)))
tax_table(transform_sample_counts(tax_glom(metagenome.pseq, taxrank="Kingdom"), 
                                  function(x) 100*x/sum(x)))

# Archaea, viruses, and viroids are a small fraction of the sample, so remove theme from subsequent analyses
metagenome.pseq <- subset_taxa(metagenome.pseq, Kingdom == "Bacteria")
rm(metagenome.otu, metagenome.taxa)

# Compare taxon relative abundances between metagenome samples and 16S rRNA amplicon samples (genus level)
pseq_glom <- tax_glom(metagenome.pseq, taxrank="Genus", NArm = FALSE)


#### IMPORT AND PROCESS ABUNDANCE DATA (RAW COUNTS, RPKM) ####
# a) Contig abundances extracted from BAM file ####
read_counts_SKE1 <- read.table("~/Downloads/sponges_metagenomes/counts/EmuEnvSum1_S7.txt", sep="\t")
read_counts_SKE4 <- read.table("~/Downloads/sponges_metagenomes/counts/EmuEnvSum4_S8.txt", sep="\t")
read_counts_SKE5 <- read.table("~/Downloads/sponges_metagenomes/counts/EmuEnvSum5_S9.txt", sep="\t")

read_counts_SKEW1 <- read.table("~/Downloads/sponges_metagenomes/counts/WatSum1_S10.txt", sep="\t")
read_counts_SKEW4 <- read.table("~/Downloads/sponges_metagenomes/counts/WatSum4_S11.txt", sep="\t")
read_counts_SKEW5 <- read.table("~/Downloads/sponges_metagenomes/counts/WatSum5_S12.txt", sep="\t")

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
metagenome.gene_blast <- read.csv("~/Downloads/sponges_metagenomes/cog-raw.txt", sep="\t", header=FALSE)
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

# c) COG abundances by summing all genes that map to the same COG ####
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

cog_to_class <- read.csv("~/Downloads/sponges_metagenomes/cog_contigs_proteins_function_categories.txt", sep="\t")
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

# d) COG class abundances (RPKM only) ####
temp <- merge(metagenome.cog_rpkm,
              metagenome.cog_annotations,
              by="COG", all=TRUE)
temp <- temp %>% dplyr::group_by(Class) %>% dplyr::summarise_if(is.numeric, sum)
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

#### EXPLORE DIFFERENCES AMONG SAMPLES BASED ON COGs AND COG CLASSES ####
# a) Ordination based on COG RPKM abundances ####
metagenome.sample_data <- read.csv("~/Downloads/sponges_metagenomes/sample_metadata.csv")

cog_comm_matrix <- as.data.frame(t(metagenome.cog_rpkm[,c("SKE1","SKE4","SKE5", "SKEW1", "SKEW4", "SKEW5")]))
cog_ord_PCA <- rda(cog_comm_matrix)
summary(cog_ord_PCA) # 97.37% on Axis 1, 1.84% on Axis 2

metagenome.sample_data <- cbind(metagenome.sample_data,
                                scores(cog_ord_PCA, display="sites", choices = c(1:3)))
colnames(metagenome.sample_data)[c(4:6)] <- c("COG_PC1", "COG_PC2", "COG_PC3")

ggplot(data = metagenome.sample_data, aes(COG_PC1, COG_PC2)) + 
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
  labs(x="PC1 (97.3%)", y="PC2 (1.84%)", title="PCA - COG abundances (RPKM)")

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

# b) Ordination based on COG class RPKM abundances ####
cog_comm_matrix <- as.data.frame(t(metagenome.class_rpkm[,c("SKE1","SKE4","SKE5", "SKEW1", "SKEW4", "SKEW5")]))
cog_ord_PCA <- rda(cog_comm_matrix)
summary(cog_ord_PCA) # 99.65% on Axis 1, 0.29% on Axis 2

metagenome.sample_data <- cbind(metagenome.sample_data,
                                scores(cog_ord_PCA, display="sites", choices = c(1:3)))
colnames(metagenome.sample_data)[c(10:12)] <- c("Class_PC1", "Class_PC2", "Class_PC3")

ggplot(data = metagenome.sample_data, aes(Class_PC1, Class_PC2)) + 
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
  labs(x="PC1 (99.7%)", y="PC2 (0.29%)", title="PCA - COG class abundances (RPKM)")

# Log-transform RPKM values before ordination.
cog_comm_matrix <- as.data.frame(t(metagenome.class_rpkm[,c("SKE1","SKE4","SKE5", "SKEW1", "SKEW4", "SKEW5")]))
cog_comm_matrix <- log(cog_comm_matrix+0.01)
cog_ord_PCA <- rda(cog_comm_matrix)
summary(cog_ord_PCA) # 96.87% on Axis 1, 2.74% on Axis 2

metagenome.sample_data <- cbind(metagenome.sample_data,
                                scores(cog_ord_PCA, display="sites", choices = c(1:3)))
colnames(metagenome.sample_data)[c(13:15)] <- c("Class_logPC1", "Class_logPC2", "Class_logPC3")

ggplot(data = metagenome.sample_data, aes(Class_logPC1, Class_logPC2)) + 
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
  labs(x="PC1 (96.9%)", y="PC2 (2.74%)", title="PCA - COG class abundances (log RPKM)")

rm(cog_comm_matrix, cog_ord_PCA)

# c) Heat map and dendrogram for COG class-level RPKM abundances ####
cog_totals_melt <- melt(metagenome.class_rpkm)
cog_totals_melt <- subset(cog_totals_melt, Class != "R" & Class !="S")
cog_totals_melt$Class <- factor(cog_totals_melt$Class)

cog_class_dendro <- as.dendrogram(hclust(vegdist(as.data.frame(t(metagenome.class_rpkm[,c(2:7)])), method="bray")))
plot(cog_class_dendro) # Order: SKEW4, SKEW1, SKEW5, SKE5, SKE1, SKE4
cog_class_dendro <- ggdendro::dendro_data(cog_class_dendro)

cog_totals_melt$value[cog_totals_melt$value > 100000] <- 100000
cog_totals_melt$variable <- factor(cog_totals_melt$variable, levels=c("SKEW4", "SKEW1", "SKEW5", "SKE5", "SKE1", "SKE4"))

ggplot(cog_totals_melt, aes(x=variable, y=Class)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient2() +
  scale_y_discrete(limits=rev(levels(cog_totals_melt$Class))) +
  geom_segment(data=ggdendro::segment(cog_class_dendro), aes(x=x, y=7*y+23.5, xend=xend, yend=7*yend+23.5)) +
  labs(x="\nSampleID", y="COG Class", fill="RPKM", title="Heatmap dendrogram - raw RPKM") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(color="black", size = 10),
        axis.title=element_text(size=11),
        plot.tag=element_text(face="bold"),
        legend.title = element_text(size=11, color="black"),
        legend.position="right",
        legend.justification=c(0.5,0.5),
        plot.title=element_text(size=11, color="black"))

rm(cog_class_dendro, cog_totals_melt)

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

# b) Compare sponge and water RPKM with differential abundance in mind ####
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
# a) Define a function to count the number of genes in each COG (including for genes assigned to multiple classes) ####
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

# b) Summarize the number of COGs in each class, for each group (total, diff abund, more abund in sponge, more abund in water) ####
# Calculate the total number of COGs in each group.
metagenome.cog_counts <- merge(metagenome.cog_counts,
                               metagenome.cog_annotations,
                               by="COG", all=TRUE)
metagenome.class_sums <- cog_counter(metagenome.cog_counts)
metagenome.cog_counts[,c("Gene","Name","Class")] <- NULL

# Subset to ALL differentially abundant COGs (p < 0.05)
tempA <- subset(metagenome.cog_diff_abund, abs(logFC) > 1 & p.adj < 0.01)
tempA <- cog_counter(tempA)
colnames(tempA) <- c("Class", "All")

# Subset to all COGs more abundant in sponges.
tempB <- subset(metagenome.cog_diff_abund, logFC < -1 & p.adj < 0.01)
tempB <- cog_counter(tempB)
colnames(tempB) <- c("Class", "Sponge")

# Subset to all COGs more abundant in water.
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
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
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
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
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
  rm(total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
}

colnames(metagenome.class_fisher_test) <- c("Class", "All", "Sponge", "Water")

metagenome.class_fisher_test$All.adj <- p.adjust(metagenome.class_fisher_test$All, method="BH")
metagenome.class_fisher_test$Sponge.adj <- p.adjust(metagenome.class_fisher_test$Sponge, method="BH")
metagenome.class_fisher_test$Water.adj <- p.adjust(metagenome.class_fisher_test$Water, method="BH")

#### EXPLORE DIFFERENTIAL ABUNDANCE BASED ON FUNCTIONALLY INTERESTING GROUPS ####
# a) Create heatmap/dendrogram combinations for easy visual comparison ####
sorting_data <- read.csv("~/Documents/sponges_final/prev_paper_data.csv")
sorting_list <- list()
for(i in 1:ncol(sorting_data)){
  sorting_list[[i]] <- subset(metagenome.cog_diff_abund, COG %in% sorting_data[,i])
}
names(sorting_list) <- colnames(sorting_data)

cog_abund_matr <- as.data.frame(metagenome.cog_rpkm[,c("COG", "SKE1", "SKE4", "SKE5", "SKEW1", "SKEW4", "SKEW5")])
rownames(cog_abund_matr) <- cog_abund_matr$COG
cog_abund_matr$COG <- NULL

sorting_list_maps <- list()
temp_dendros <- list()
sorting_melt_panels <- list()

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
  
  sorting_melt_panels[[i]] <- cog_abund_melt
  
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

rm(cog_abund_matr, cog_abund_matr_subset, cog_abund_melt)

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

##################################################### METAGENOME ANALYSIS II: MAGS ###########################################
#### IMPORT DATA AND RUN AN EXPLORATORY PCA ####
# Import MAG data and count number of hits per MAG per COG.
mag.sample_data <- read.csv("~/Downloads/sponges_metagenomes/mag_metadata.csv")
mag.sample_data <- mag.sample_data[order(mag.sample_data$MAG), ]
rownames(mag.sample_data) <- mag.sample_data$MAG

mag.cog_counts <- read.csv("~/Downloads/sponges_metagenomes/2020_09_30_bin_proteins_COGs.csv")
mag.cog_counts$Group <- NULL
mag.cog_counts$Name <- NULL

mag.cog_counts <- mag.cog_counts %>% dplyr::group_by(Bin, COG) %>% dplyr::count()
mag.cog_counts <- reshape2::dcast(mag.cog_counts, COG ~ Bin)
rownames(mag.cog_counts) <- mag.cog_counts$COG
mag.cog_counts$COG <- NULL
mag.cog_counts[is.na(mag.cog_counts)=="TRUE"] <- 0

# Run a PCA to see if MAGs cluster separately.
mag.pca <- mag.cog_counts
for(i in 1:ncol(mag.pca)){
  mag.pca[,i] <- 100*mag.pca[,i]/sum(mag.pca[,i])
}

mag.pca <- vegan::rda(as.data.frame(t(mag.pca)))

mag.sample_data <- cbind(mag.sample_data,
                         vegan::scores(mag.pca, choices=c(1:3), display="sites"))

ggplot(mag.sample_data, aes(x=PC1, y=PC2, color=Type)) + 
  geom_point(size=3)  #+
  #geom_text(aes(label=MAG), color="black")


#### COMPARE MAG GENETIC FUNCTIONS WITH SPONGE-ENRICHED COGS ####
mag.levels <- c("MAG05", "MAG08", "MAG07", "MAG25", "MAG24", "MAG03",
                "MAG06", "MAG04", "MAG10", "MAG12", "MAG09", "MAG18", # sponge-enriched
                "MAG02", "MAG15", "MAG19", "MAG11", "MAG13", "MAG17",
                "MAG23", "MAG16", "MAG22", "MAG01", "MAG14", "MAG21", "MAG20") # water-enriched

mag.interesting_cog_lists <- list()
for(i in 1:ncol(sorting_data)){
  temp <- subset(mag.cog_counts, rownames(mag.cog_counts) %in% sorting_data[,i])
  temp$COG <- rownames(temp)
  temp <- melt(temp)
  colnames(temp) <- c("COG", "MAG", "Value")
  temp <- merge(temp, metagenome.cog_annotations, by="COG", all=FALSE)
  temp$MAG <- factor(temp$MAG, levels=mag.levels)
  
  mag.interesting_cog_lists[[i]] <- temp
  rm(temp)
}

mag.intesting_cog_plots <- list()
for(i in 1:length(mag.interesting_cog_lists)){
  plot1 <- ggplot(mag.interesting_cog_lists[[i]], aes(x=MAG, y=COG)) +
    geom_tile(aes(fill=Value)) + 
    scale_fill_gradient(low="white", high="red") +
    geom_vline(xintercept=13.5) +
    scale_y_discrete(limits=rev(levels(mag.interesting_cog_lists[[i]]$COG))) +
    scale_x_discrete(limits=rev(levels(mag.interesting_cog_lists[[i]]$MAG))) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, color="black"),
          axis.text.y=element_text(color="black"),
          legend.position="bottom",
          panel.grid=element_blank())
  
  plot2 <- ggplot(sorting_melt_panels[[i]], aes(x=SampleID, y=Full.Name)) + 
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
  
  plots <- list(plot1, plot2)
  grobs <- list()
  heights <- list()
  for (i in 1:length(plots)){
    grobs[[i]] <- ggplotGrob(plots[[i]])
    heights[[i]] <- grobs[[i]]$heights[1:12]
  }
  maxheight <- do.call(grid::unit.pmax, heights)
  for (i in c(1,2)){
    grobs[[i]]$heights[1:12] <- as.list(maxheight)
  }
  rm(maxheight, heights, plots)
  
  mag.intesting_cog_plots[[i]] <- arrangeGrob(grobs[[1]],
                                              grobs[[2]],
                                              ncol=2)
  
  rm(plot1, plot2, grobs)
}