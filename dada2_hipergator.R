library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(vegan)
library(indicspecies)
library(ALDEx2)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#### 1st sequencing run
path <- "/blue/juliemeyer/acauvin/dendrogyra_run2/dendrogyra_run2_all/cutadapt"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Perform filtering and trimming
filt_path <- file.path(path, "filtered") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# Learn the Error Rates, it TAKES TIME!
exists<- file.exists(filtFs) & file.exists(filtRs)
errF <- learnErrors(filtFs[exists], multithread=TRUE)
errR <- learnErrors(filtRs[exists], multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# Dereplicate the filtered fastq files
exists<- file.exists(filtFs) & file.exists(filtRs)
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Inspecting the dada-class object returned by dada:
dadaFs[[1]]

# Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Inspect distribution of sequence lengths - Track reads through the pipeline
table(nchar(getSequences(seqtab)))
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track) <- sample.names
head(track)
write.table(track, "dada_read_stats1.txt",sep="\t",col.names=NA)
saveRDS(seqtab, "/blue/juliemeyer/acauvin/dendrogyra_run2/dendrogyra_run2_all/seqtab.rds")

taxa <- assignTaxonomy(seqtab.nochim, "/blue/juliemeyer/share/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# FIX the NAs in the taxa table
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]
write.table(taxon,"silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "silva_otu_table.txt",sep="\t",col.names=NA)

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("silva_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata.txt",sep="\t",header=T, row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
sample_names(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps

# remove chloroplasts and mitochondria and Eukaryota
get_taxa_unique(ps, "Family") 
get_taxa_unique(ps, "Order") 
get_taxa_unique(ps, "Kingdom") 
ps <- subset_taxa(ps, Family !="Mitochondria")
ps <- subset_taxa(ps, Order !="Chloroplast")
ps <- subset_taxa(ps, Kingdom !="Eukaryota")
ps <- subset_taxa(ps, Kingdom !="NA")
get_taxa_unique(ps, "Family") 
get_taxa_unique(ps, "Order") 
get_taxa_unique(ps, "Kingdom")
ps

# Now export cleaned otu and taxa tables from phyloseq for future reference
otu = as(otu_table(ps), "matrix")
taxon = as(tax_table(ps), "matrix")
metadata = as(sample_data(ps), "matrix")
write.table(otu,"silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)

# export ASV table as relative abundance
ps_ra<-transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
otu_ra = as(otu_table(ps_ra), "matrix")
write.table(otu_ra,"silva_nochloronomito_otu_table_RA.txt",sep="\t",col.names=NA)
