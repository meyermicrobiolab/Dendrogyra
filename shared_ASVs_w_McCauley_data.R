library(metagMisc)
library(nationalparkcolors)
library(phyloseq)
library(microbiome)
library(CoDaSeq)
library(zCompositions)
library(cowplot)
library(googleway)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(vegan)
library(tidyverse)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(corncob)
library(DESeq2)
library(Biostrings)
library(seqateurs)
library(DECIPHER)
library(phangorn)
library(seqinr)
library(ape)
library(ggtree)
library(purrr)
library(digest)


# use an ASV table that has been filtered to remove low-abundance ASVs, no blanks
otu <- read.table("silva_nochloronomito_otu_table_ps5_filt.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5_filt.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_noblanks_filt.txt",sep="\t",header=T,row.names=1)
samples$Region[samples$Region == 'Curacao'] <- 'Curaçao'
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps5 <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                sample_data(samples), 
                tax_table(taxon))

#############################################################
#Endozoicomonas Phylogenetic Analyses
############################################################
#Read in McCauley Phyloseq Object
v4<- readRDS("V4 Library.rdata")
v4

#Have to subset the full V4 library by those that only use the same primers
v4<- subset_samples(v4, Primer_Region=="515F-806R")
v4

#Subset by Atlantic Ocean
unique(v4@sam_data$Original_Waterbody)
v4<- subset_samples(v4, Original_Waterbody=="Southern Atlantic Ocean" | Original_Waterbody=="Caribbean Sea" | 
                      Original_Waterbody== "Northern Atlantic Ocean")
v4

v4<- subset_samples(v4, Collection=="Field" )
v4

v4<- subset_samples(v4, Condition=="Healthy")
v4

#Prune taxa no longer in samples
v4 <- prune_taxa(taxa_sums(v4) > 0, v4) 
v4 

mccauley_endo<- subset_taxa(v4, Genus == "Endozoicomonas")
mccauley_endo #51 ASVs

colnames(mccauley_endo@otu_table) <- paste('McCauley', colnames(mccauley_endo@otu_table), sep = '_')
row.names(mccauley_endo@tax_table)<- colnames(mccauley_endo@otu_table)
mccauley_endo@refseq@ranges@NAMES<- taxa_names(mccauley_endo)

dendro_endo<- subset_taxa(ps5, Genus == "Endozoicomonas")
dendro_endo #7 ASVs
dendro_seq_set<-DNAStringSet(row.names(dendro_endo@tax_table))
dendro_seq_set
colnames(dendro_endo@otu_table) <- paste('Dendro_ASV',1:7, sep="")
row.names(dendro_endo@tax_table)<- colnames(dendro_endo@otu_table)
names(dendro_seq_set)<- taxa_names(dendro_endo)
dendro_seq_set
view(dendro_endo@tax_table)
dendro_endo<- merge_phyloseq(dendro_endo, dendro_seq_set)

endozoicomonas_all_otu_table<- merge_phyloseq_pair(mccauley_endo@otu_table, dendro_endo@otu_table)
endozoicomonas_all_tax_table<- merge_phyloseq_pair(mccauley_endo@tax_table, dendro_endo@tax_table)
endozoicomonas_all_sam_data<- merge_phyloseq_pair(mccauley_endo@sam_data, dendro_endo@sam_data)
endozoicomonas_all_refseq<- merge_phyloseq_pair(mccauley_endo@refseq, dendro_endo@refseq)

endozoicomonas_all<- merge_phyloseq(endozoicomonas_all_otu_table, endozoicomonas_all_tax_table, 
                                    endozoicomonas_all_sam_data, endozoicomonas_all_refseq)

endo_all_fas<- ps_to_fasta(endozoicomonas_all, out.file = "endozoicomonas_all_fasta.fa", seqnames = "unique")

endo_all_fas<- "~/Documents/dendrogyra_all_updated/filtered_by_500/endozoicomonas_all_fasta.fa"
endo_all_seqs <- readDNAStringSet(endo_all_fas)
names(endo_all_seqs)[1:51] <- paste(taxa_names(mccauley_endo))
names(endo_all_seqs)[52:58] <- paste(taxa_names(dendro_endo))
names(endo_all_seqs)

seqs <- OrientNucleotides(c(endo_all_seqs))
aligned <- AlignSeqs(seqs)
BrowseSeqs(aligned, highlight=0)

# Identify trimming positions
start_trim <- 19  # Trim from position 19 at the start
end_trim <- 227   # Trim up to position 100 at the end
trimmed_seqs <- subseq(aligned, start=start_trim, end=end_trim)
BrowseSeqs(trimmed_seqs, highlight=0)

unique_seqs <- unique(trimmed_seqs)
identical_seq_groups <- vector("list", length(unique_seqs))
for (i in seq_along(unique_seqs)) {
  identical_indices <- which(trimmed_seqs == unique_seqs[i])
  identical_seq_groups[[i]] <- names(trimmed_seqs)[identical_indices]
}

names(identical_seq_groups) <- sapply(identical_seq_groups, function(x) x[1])

identical_seq_groups

#Get McCauley samples shared w Dendrogyra
asvs_of_interest<- c("McCauley_ASV73", "McCauley_ASV630", "McCauley_ASV1157", "McCauley_ASV20697",
         "McCauley_ASV28684")

otu_table_endozoicomonas<- otu_table(endozoicomonas_all)
otu_table_endozoicomonas<- as.data.frame(otu_table_endozoicomonas)

df_subset <-subset(otu_table_endozoicomonas, select = asvs_of_interest)
# Calculate row sums for numeric columns only
rows_sum_gt_zero <- rowSums(df_subset) > 0
# Subset the dataframe
df_subset <- df_subset[rows_sum_gt_zero, ]
endo_samples<- row.names(df_subset)

# Subset the phyloseq object
ps_subset <- prune_samples(endo_samples, endozoicomonas_all)

species_w_shared_endo<- unique(ps_subset@sam_data$Scientific_Name)
print(species_w_shared_endo)

#Make a venn diagram

library(VennDiagram)

pdf("Endozoicomonas_ASVs_Shared_venndiagram.pdf", height=8, width=11)
grid.newpage() 
draw.pairwise.venn(area1=45, area2=6,cross.area=3,fill=c("#432371","#faae7b"))
dev.off()

####################################################
#Candidatus Amoebophilus ASVs
####################################################

# use an ASV table that has been filtered to remove low-abundance ASVs, no blanks
otu <- read.table("silva_nochloronomito_otu_table_ps5_filt.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5_filt.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_noblanks_filt.txt",sep="\t",header=T,row.names=1)
samples$Region[samples$Region == 'Curacao'] <- 'Curaçao'
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps5 <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                sample_data(samples), 
                tax_table(taxon))

#############################################################
#Candidatus Phylogenetic Analyses
############################################################
#Read in McCauley Phyloseq Object
v4<- readRDS("V4 Library.rdata")
v4

#Have to subset the full V4 library by those that only use the same primers
v4<- subset_samples(v4, Primer_Region=="515F-806R")
v4

#Subset by Atlantic Ocean
unique(v4@sam_data$Original_Waterbody)
v4<- subset_samples(v4, Original_Waterbody=="Southern Atlantic Ocean" | Original_Waterbody=="Caribbean Sea" | 
                      Original_Waterbody== "Northern Atlantic Ocean")
v4

v4<- subset_samples(v4, Collection=="Field" )
v4

v4<- subset_samples(v4, Condition=="Healthy")
v4

#Prune taxa no longer in samples
v4 <- prune_taxa(taxa_sums(v4) > 0, v4) 
v4 

mccauley_candi<- subset_taxa(v4, Genus == "Candidatus Amoebophilus")
mccauley_candi #45 ASVs
colnames(mccauley_candi@otu_table) <- paste('McCauley', colnames(mccauley_candi@otu_table), sep = '_')
row.names(mccauley_candi@tax_table)<- colnames(mccauley_candi@otu_table)
mccauley_candi@refseq@ranges@NAMES<- taxa_names(mccauley_candi)

dendro_candi<- subset_taxa(ps5, Genus == "Candidatus Amoebophilus")
dendro_candi #5 ASVs
dendro_seq_set<-DNAStringSet(row.names(dendro_candi@tax_table))
dendro_seq_set
colnames(dendro_candi@otu_table) <- paste('Dendro_ASV',1:5, sep="")
row.names(dendro_candi@tax_table)<- colnames(dendro_candi@otu_table)
names(dendro_seq_set)<- taxa_names(dendro_candi)
dendro_seq_set
view(dendro_candi@tax_table)
dendro_candi<- merge_phyloseq(dendro_candi, dendro_seq_set)

candidatus_all_otu_table<- merge_phyloseq_pair(mccauley_candi@otu_table, dendro_candi@otu_table)
candidatus_all_tax_table<- merge_phyloseq_pair(mccauley_candi@tax_table, dendro_candi@tax_table)
candidatus_all_sam_data<- merge_phyloseq_pair(mccauley_candi@sam_data, dendro_candi@sam_data)
candidatus_all_refseq<- merge_phyloseq_pair(mccauley_candi@refseq, dendro_candi@refseq)

candidatus_all<- merge_phyloseq(candidatus_all_otu_table, candidatus_all_tax_table, 
                                    candidatus_all_sam_data, candidatus_all_refseq)

candi_all_fas<- ps_to_fasta(candidatus_all, out.file = "candidatus_all_fasta.fa", seqnames = "unique")

candi_all_fas<- "~/Documents/dendrogyra_all_updated/filtered_by_500/candidatus_all_fasta.fa"
candi_all_seqs <- readDNAStringSet(candi_all_fas)
names(candi_all_seqs)[1:45] <- paste(taxa_names(mccauley_candi))
names(candi_all_seqs)[46:50] <- paste(taxa_names(dendro_candi))
names(candi_all_seqs)

seqs <- OrientNucleotides(c(candi_all_seqs))
aligned <- AlignSeqs(seqs)
BrowseSeqs(aligned, highlight=0)

# Identify trimming positions
start_trim <- 19  # Trim from position 19 at the start
end_trim <- 227   # Trim up to position 100 at the end
trimmed_seqs <- subseq(aligned, start=start_trim, end=end_trim)
BrowseSeqs(trimmed_seqs, highlight=0)

unique_seqs <- unique(trimmed_seqs)
identical_seq_groups <- vector("list", length(unique_seqs))
for (i in seq_along(unique_seqs)) {
  identical_indices <- which(trimmed_seqs == unique_seqs[i])
  identical_seq_groups[[i]] <- names(trimmed_seqs)[identical_indices]
}

names(identical_seq_groups) <- sapply(identical_seq_groups, function(x) x[1])

identical_seq_groups

#Get McCauley samples shared w Dendrogyra
asvs_of_interest<- c("McCauley_ASV124", "McCauley_ASV268", "McCauley_ASV258", "McCauley_ASV4975")

otu_table_candidatus<- otu_table(candidatus_all)
otu_table_candidatus<- as.data.frame(otu_table_candidatus)

df_subset <-subset(otu_table_candidatus, select = asvs_of_interest)
# Calculate row sums for numeric columns only
rows_sum_gt_zero <- rowSums(df_subset) > 0
# Subset the dataframe
df_subset <- df_subset[rows_sum_gt_zero, ]
candi_samples<- row.names(df_subset)

# Subset the phyloseq object
ps_subset <- prune_samples(candi_samples, candidatus_all)

species_w_shared_candi<- unique(ps_subset@sam_data$Scientific_Name)
print(species_w_shared_candi)

pdf("Candidatus_ASVs_Shared_venndiagram.pdf", height=8, width=11)
grid.newpage() 
draw.pairwise.venn(area1=43, area2=5,cross.area=3,fill=c("#432371","#faae7b"))
dev.off()

library(clipr)
write_clip(mccauley_endo@refseq$McCauley_ASV73)
write_clip(mccauley_endo@refseq$McCauley_ASV1157)
write_clip(mccauley_endo@refseq$McCauley_ASV28684)

write_clip(mccauley_candi@refseq$McCauley_ASV124)
write_clip(mccauley_candi@refseq$McCauley_ASV258)
write_clip(mccauley_candi@refseq$McCauley_ASV4975)







