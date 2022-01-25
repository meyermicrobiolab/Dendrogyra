#Load Libraries
library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(vegan)
library(knitr)
library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(cowplot)
library(randomcoloR)
library(dplyr)
library(reshape2)
library(tibble)
library(exactRankTests)
library(nlme)
library(data.table)
library(Rmisc)
library(indicspecies)
library(viridis)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#Combine RDS files
run1 <- readRDS("/blue/juliemeyer/acauvin/dendrogyra_run1/dendrogyra_run1_all/seqtab.rds") 
run2 <- readRDS("/blue/juliemeyer/acauvin/dendrogyra_run2/dendrogyra_run2_all/seqtab.rds") 
run3 <- readRDS("/blue/juliemeyer/acauvin/dendrogyra_run3/dendrogyra_run3_all/seqtab.rds") 

st.all <- mergeSequenceTables(run1, run2, run3, repeats="sum") # You may get the message "Duplicated sample names detected in the sequence table row names." to let you know that there are duplicate names across samples - it is not an error, just a message. If you have run a sample on more than one sequencing run, the ASV counts will be added together.
#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(st.all)
# Combine read stats from 4 runs and add chimera summary - This does not combine duplicate rows from repeated samples.
stat1 <- read.table("/blue/juliemeyer/acauvin/dendrogyra_run1/dada_read_stats1.txt",sep="\t",header=TRUE, row.names=1)
stat2 <- read.table("/blue/juliemeyer/acauvin/dendrogyra_run2/dada_read_stats1.txt",sep="\t",header=TRUE, row.names=1)
stat3 <- read.table("/blue/juliemeyer/acauvin/dendrogyra_run3/dendrogyra_run3_all/dada_read_stats1.txt",sep="\t",header=TRUE, row.names=1)

stats.all<-bind_rows(stat1, stat2,stat3) 
write.table(stats.all, "dada_read_stats_all.txt",sep="\t",col.names=NA)
# Track reads through the pipeline
# As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline
rowSums(seqtab.nochim)
# need to write this out to add to dada read stats
# SAVE the non-chimeric sequence variant table SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS
saveRDS(seqtab.nochim, file="~/blue/juliemeyer/acauvin/dendrogyra_run_all/dendrogyra_seqtab_all.rds")



###Assign Taxonomy
## Assign taxonomy in DADA2


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
samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
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
ps #20179 taxa and 197 samples
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