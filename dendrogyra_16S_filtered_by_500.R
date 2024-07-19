#Load libraries
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
library(bioconductor)

names(park_palettes)

pal<- park_palette("GrandTeton")
GrandTeton = c("#F0EEE2", "#5B6C88", "#48594E", "#A8D0CF", "#BABBB1")


#Load the data and create the phyloseq object
# Create phyloseq object from otu and taxonomy tables without chloroplast and mitochondria
otu <- read.table("silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)
samples$Region[samples$Region == 'Curacao'] <- 'Curaçao'
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps


#Remove problematic samples from ps
ps<- remove_samples(c("17818", "17797"), ps)
summarize_phyloseq(ps)


#Remove samples that have less than 500 reads
ps <- prune_samples(sample_sums(ps) >= 500, ps)
ps #7559 taxa across 132 samples
#Prune taxa no longer in samples
ps <- prune_taxa(taxa_sums(ps) > 0, ps) 
ps #7451 taxa in 132 samples

summarize_phyloseq(ps)
ps_ra<- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))


#Get number of unique taxa
get_taxa_unique(ps, "Phylum") #54
get_taxa_unique(ps, "Order") #307
get_taxa_unique(ps, "Class") #130
get_taxa_unique(ps, "Family") #544
length(get_taxa_unique(ps, "Genus")) #1054
ps5<-filter_taxa(ps, function(x) mean(x) >5, TRUE)
ntaxa(ps5) #183
ps10<-filter_taxa(ps, function(x) mean(x) >10, TRUE)
ntaxa(ps10) #102
get_taxa_unique(ps, "Genus") #1054
get_taxa_unique(ps5, "Genus") #97
get_taxa_unique(ps10, "Genus") #65

# Now export filtered otu and taxa tables from phyloseq for future reference after removing sample blanks

#Write out ps5 to tables
otu_ps5 = as(otu_table(ps5), "matrix")
taxon_ps5 = as(tax_table(ps5), "matrix")
metadata = as(sample_data(ps5), "matrix")
write.table(otu_ps5,"silva_nochloronomito_otu_table_ps5_filt.txt",sep="\t",col.names=NA)
write.table(taxon_ps5,"silva_nochloronomito_taxa_table_ps5_filt.txt",sep="\t",col.names=NA)
write.table(metadata,"metadata_noblanks_filt.txt",sep="\t",col.names=NA)
ps5_ra<-transform_sample_counts(ps5, function(OTU) OTU/sum(OTU))
otu_ps5_ra = as(otu_table(ps5_ra), "matrix")
write.table(otu_ps5_ra,"silva_nochloronomito_otu_table_ps5_RA_filt.txt",sep="\t",col.names=NA)


ps10 #102 taxa and 132 samples
otu_ps10 = as(otu_table(ps10), "matrix")
taxon_ps10 = as(tax_table(ps10), "matrix")
write.table(otu_ps10,"silva_nochloronomito_otu_table_ps10_filt.txt",sep="\t",col.names=NA)
write.table(taxon_ps10,"silva_nochloronomito_taxa_table_ps10_filt.txt",sep="\t",col.names=NA)


###########################################################
###########################################################
#Beta Diversity Calculations
###########################################################
###########################################################

## Perform center-log-ratio transformation on ASVs and calculate Aitchison Distance and principal components

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


# First, replace 0 values with an estimate (because normalization is taking log, can't have 0)
# Also transposing here, need samples as rows
d.czm <- cmultRepl(t(otu), method="CZM", label=0, z.warning = 0.999)
# Perform the center-log-ratio (CLR) transformation 
d.clr <- codaSeq.clr(d.czm)
# transpose matrix of CLR transformed data for ordination and dendrogram
E.clr <- t(d.clr)
# plot compositional PCA biplot (perform a singular value decomposition)
d.pcx <- prcomp(E.clr)
# calculate percent variance explained for the axis labels
pc1 <- round(d.pcx$sdev[1]^2/sum(d.pcx$sdev^2),2)
pc2 <- round(d.pcx$sdev[2]^2/sum(d.pcx$sdev^2),2)
xlab <- paste("PC1: ", pc1, sep="")
ylab <- paste("PC2: ", pc2, sep="")
#biplot(d.pcx, cex=c(0.6,0.4), var.axes=F,scale=1, xlab=xlab, ylab=ylab)
summary(d.pcx)
str(d.pcx)
screeplot(d.pcx)
# replot PCA with ggplot2 (showing samples only)
df_out <- as.data.frame(d.pcx$x)
theme_set(theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))

cols<-c("Belize"="#003459", "Curaçao"="#75dddd")

pdf("PCA_Region_filt_500.pdf",bg = "white",width=8.5)
p<-ggplot(df_out,aes(x=PC1,y=PC2,fill=samples$Region, shape=samples$Region))
p<-p+geom_point( size=3.5)+
  theme(axis.title = element_text(size=14))+
  theme(axis.text=element_text(size=12))+
  theme(legend.title = element_text(size=14))+
  theme(legend.text = element_text(size=12))+
  theme(legend.position = "top")+
  theme(legend.title = element_blank()) +
  theme(legend.key.width = unit(1.0,  unit = "cm")) +
  guides(fill = guide_legend(byrow = TRUE)) 
p<- p+ labs(x=xlab, y=ylab, fill="Region") + coord_fixed()+
  labs(fill = "Region", shape="Region")+
  scale_fill_manual(values=cols, labels = c('Belize', 'Curaçao'))+
  scale_shape_manual(values=c(21,22), labels = c('Belize', 'Curaçao'))+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))
p
dev.off()

####### Use phyloseq/vegan to perform ANOSIM/PERMANOVA

# set metadata as factors for anosim
region<- as.character(samples$Region)

# anosim between groups by Region using Aitchison distance
dist.clr <- dist(E.clr)
ano.region <- anosim(dist.clr, region, permutations=999)
ano.region
pdf("DCYL_Healthy_Region_filt.pdf",width=8.5)
plot(ano.region)
dev.off()
print(ano.region)

###################PERMANOVA#################
dist.clr <- dist(E.clr)
perm<-adonis2(dist.clr~region,as(sample_data(ps5),"data.frame"))
print(perm)

write.table(perm, "permanova_results.txt")

### Beta Diversity Dispersions - Region
#calculate multivariate dispersions based on Region
mod <-betadisper(dist.clr, region)
#one way anova
anova(mod)
#boxplots
pdf("beta_dispersion_Region_filt_500.pdf",width=8.5)
plot(mod)
dev.off()

pdf("beta_dispersion_Region_filt_500.pdf",width=8.5)
boxplot(mod)
dev.off()


## Compute mean distance to centroid per group
#this just prints values on the console
tapply(mod$distances, region, mean)
## Same, but variance instead
tapply(mod$distances, region, var)

#Get the distances to centroid from the model
mod$distances
dis <- mod$distances
#melt
dis.melt <- melt(dis)
#move rownames to columns so we can merge the dispersion values and metadata
dis.melt$Sample <- rownames(dis.melt)
samples$Sample <- rownames(samples)
#merge metadata and dispersion 
dis.region <- merge(samples, dis.melt)
#rename column
colnames(dis.region)[6] <- "distance"

#run linear model to test significance of region
distlm <-lm(distance~Region, data=dis.region)
summary(distlm)
anova(distlm)

t.test(distance~Region, data=dis.region)
wilcox.test(distance~Region, data = dis.region)

#plot average dispersion by group, with all points shown                               

pdf("Region_DistanceToCentroid_filt_500.pdf",width=11,height=8.5)
p2<-ggplot(dis.region, aes(x=Region, y=distance)) +
  geom_boxplot(mapping =aes(x=Region,y=distance, fill=Region, color="black"), alpha = 0.4,  outlier.shape=NA) +
  geom_point(position=position_jitter(width=.1, height=0),aes(colour="black",fill=Region),size=4, pch=21, colour = "black") +
  scale_fill_manual(values = cols, labels = c('Belize', 'Curaçao')) +
  scale_color_manual(values = c( "black")) +
  ylab("Distance to Centroid")+
  xlab("Region")+
  guides(color="none") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(colour="black"))+
  theme(axis.text.y = element_text(colour="black"))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  theme(strip.text.x=element_text(size=16, face="bold"))+
  theme(legend.position = "none")
p2
dev.off()


###############################################
###############################################
#Taxonomic Differences
###############################################
###############################################

#### Making bar charts using ps5 (low abundance filtered) with n samples (n taxa) 

ps5 #183 taxa in 132 samples
ps5_ra
#figure out how many colors you need
get_taxa_unique(ps5_ra, "Order") #62
get_taxa_unique(ps5_ra, "Class") #35
get_taxa_unique(ps5_ra, "Family") #76


#Get top 20 taxa
ps20<-prune_taxa(names(sort(taxa_sums(ps5_ra),TRUE)[1:20]), ps5_ra)
ps20

#Coerce Top 20 taxa as data frame
# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps20), "matrix")
# transpose if necessary
if(taxa_are_rows(ps20)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
colnames(OTUdf)
top20<- colnames(OTUdf)
write(top20, "top20ASVs.txt")

#Extract Top 20 ASV table as a dataframe
ASV1<- as(tax_table(ps20), "matrix")
ASVdf<- as.data.frame(ASV1)
write.csv(ASVdf, "top20Taxa.csv", row.names=TRUE)

#Write out ps20
otu.ps20 = as(otu_table(ps20), "matrix")
taxon.ps20 = as(tax_table(ps20), "matrix")
write.table(otu.ps20,"ps20_ASVs.txt",sep="\t",col.names=NA)
write.table(taxon.ps20,"ps20_taxon.txt",sep="\t",col.names=NA)

########Creating the bubble plot

#### IMPORTANT - remove any dashes in sample names on the ASV doc before starting. 
#IMPORTANT - Also don't forget to switch your data in Excel transpose data on ASV doc so sample names are at the top.
#You will have to do this manually in Excel.

otu_tab_file <- read.table("silva_nochloronomito_otu_table_ps5_filt.txt",sep="\t",header=TRUE) #you will have to open your ps20_ASVs.txt table and transpose the data before doing this step
tax_file <- read.table("silva_nochloronomito_taxa_table_ps5_filt.txt", sep="\t",header=TRUE)
metadat<- read.table("metadata_noblanks_filt.txt", sep="\t", header=TRUE)
colnames(metadat)[1]<- "Sample.ID"

rownames(otu_tab_file) <- otu_tab_file[,1]
otu_tab_file<- otu_tab_file[,-1]
otu_tab_file<- t(otu_tab_file)
colnames(metadat)[1]<- "Sample.ID"

# List of all replicates, it's important that they are in the same order and match the names in the OTU table.
replicates_list <- metadat$Sample.ID
replicates_list

replicates_groups <- metadat$Region

# Do these two list have the same length?
length(replicates_list) == length(replicates_groups) 

# Taxonomic level that will be aggregated and displaeyed in the bubble plot.
# It should match a taxonomic level from the taxonomic annotation file.
# If you want to display OTUs/ASVs, simply write 'OTU'.
tax_aggr <- "Genus"

# Number of taxonomic bins to be displayed in the bubble plot.
tax_number <- 15

# Taxonomic level that will be used for colouring bubbles (i.e. categorical variable)
# It should match a taxonomic level from the taxonomic annotation file
# (and it should be a higher taxonomic level that the one from 'tax_aggr').
tax_col <- "Phylum"

# Filename for saving the bubble plot (svg or pdf)
file_name <- "bubble_plot2.svg"

# Dimension (in inch) of this file (first number is the width, the second number is the height)
plot_dim <- c(6,6)

otu_tab <- otu_tab_file
names(rownames(otu_tab)) <- "OTU"
colnames(otu_tab)[1:132]<- metadat$Sample.ID

otu_tab[1:5,1:5]

tax_tab <- tax_file
names(tax_tab)[1] <- "OTU"

tax_tab[1:5,1:4]

###Parse the taxonomic file
#This code removes useless annotations from the taxonomic annotation file, such as 'uncultured bacteria', or 'unkwown family', and fills the blanks by copying the previous taxonomic levels (while adding an 'unassigned' mention). Do this, even though you've already removed NAs from the data. The bubble plot won't work otherwise.

# Delete all cells containing 'uncultured' or 'unkNown'
for (col in 2:ncol(tax_tab)) {
  for (row in 1:nrow(tax_tab)) {
    if (grepl("uncultured",tax_tab[row,col],ignore.case = TRUE)) {
      tax_tab[row,col] <- ""
    }
    if (grepl("unknown",tax_tab[row,col],ignore.case = TRUE)) {
      tax_tab[row,col] <- ""
    }
  }
}

# Replace empty cells by 'NA'
tax_tab2 <- as.data.frame(apply(tax_tab, 2, function(x) gsub("^$|^ $", NA, x)))

# Remove columns containing only 'NA'
col_to_remove <- c()

for (col in 2:ncol(tax_tab2)) {
  x <- sum(is.na(tax_tab2[,col]))/nrow(tax_tab2)
  if (x == 1) {
    col_to_remove <- c(col_to_remove, col)
  }
}

if (length(col_to_remove) > 0) {
  tax_tab3 <- tax_tab2[,-col_to_remove]
} else {
  tax_tab3 <- tax_tab2
}

# Set taxonomic annotations as character variables

for (col in 2:ncol(tax_tab3)) {
  tax_tab3[,col] <- as.character(tax_tab3[,col])
}

# Fill all NAs

for (col in 2:ncol(tax_tab3)) {
  for (row in 1:nrow(tax_tab3)) {
    if (is.na(tax_tab3[row,col])) {
      if (!grepl("OTU", tax_tab3[row,col-1]) & !grepl("unassigned", tax_tab3[row,col-1])) {
        tax_tab3[row,col] <- paste0("unassigned ", tax_tab3[row,col-1])
      } else {
        tax_tab3[row,col] <- tax_tab3[row,col-1]
      }
    }
  }
}

###Compute the relative abundance of OTUs for each sample
#You can run this command even if your OTU table already contains relative abundances instead of reads counts. It won't hurt.

otu_counts <- colSums(otu_tab[,-1])
otu_tab2 <- otu_tab
otu_tab2 <- sweep(otu_tab[,-1], 2, otu_counts, `/`)
otu_tab2[is.na(otu_tab2)] <- 0

# Check that the sum of relative abundances for each sample is 1.
colSums(otu_tab2[,-1])


###Merge the OTU and taxonomic tables together
m <- merge(otu_tab2, tax_tab3)

# Has the merged table the expected dimension?
dim(m)


###Aggregate the table to taxonomic level defined in the variable 'tax_aggr'
# First, we should save in a column the taxonomic information needed for computing the bubble plot
taxonomy <- c()
for (row in 1:nrow(m)) {
  taxonomy <- c(taxonomy, paste0(m[row,names(m)==tax_col], ";", m[row,names(m)==tax_aggr]))
}

# Subset from the merged table the selected samples only
# Subset from the merged table the selected samples only
m2 <- m[,names(m) %in% replicates_list]

# Aggregate 'm2' based on the selected taxonomic level
m3 <- aggregate(m2, by=list(taxonomy), FUN=)
dim(m3)
m3[1:5,1:4]



###Sort the table by decreasing size of taxonomic groups
#Here we only keep the top n taxonomic groups, as defined in the tax_number variable. All the others taxonomic groups will be pooled together in a new bin labelled 'Other'.

# Sort the taxonomic table

if (tax_number > nrow(m3)) {
  tax_number <- nrow(m3)
}

m3$average <- rowMeans(m3[,-1])
m3.sorted <- m3[order(-m3$average),]

# Aggregate the smaller taxonomic bins together
m3.sorted$selection <- rep("discarded", nrow(m3.sorted))
m3.sorted$selection[1:tax_number] <- "retained"
m3.sorted$Group.1[m3.sorted$selection == "discarded"] <- "Other;Other"
m3.sorted$average <- NULL
m3.sorted$selection <- NULL
m4 <- aggregate(m3.sorted[,-1], by=list(taxonomy=m3.sorted$Group.1), FUN=sum)

# What is the relative abundances of the taxonomic bins that were pooled together in the 'Other' bin?
m4[m4$taxonomy == "Other;Other", -1]


mean(as.numeric(m4[m4$taxonomy == "Other;Other", -1]))

# If you find these numbers too big, you can simply increase the value of the 'tax_number' variable,
# Or alternatively choose a higher taxonomic level to display


###Transpose 'm4'
n <- m4$taxonomy
m4.t <- as.data.frame(t(m4[,-1]))
colnames(m4.t) <- n
m4.t$sample <- rownames(m4.t)
rownames(m4.t) <- NULL


###Calculate the mean and the standard deviation for each sample
m4.t$replicate <- rep(NA, nrow(m4.t))
for (line in 1:(nrow(m4.t))){
  m4.t$replicate[line] <- replicates_list[m4.t$sample[line] == replicates_list]
}

# Compute the mean
m4.t.mean <- aggregate(m4.t[,1:(ncol(m4.t)-2)],
                       by = list(m4.t$replicate),
                       FUN = "mean")
names(m4.t.mean)[1] <- "sample"

dim(m4.t.mean)                              
## [1]  132 17

# Compute the standard deviation                             
m4.t.sd <- aggregate(m4.t[,1:(ncol(m4.t)-2)],
                     by = list(m4.t$replicate),
                     FUN = "sd")
names(m4.t.sd)[1] <- "sample"

dim(m4.t.sd)


###Melt and merge the two dataframes
# Melt the dataframes
molten.mean <- melt(m4.t.mean, id.vars = "sample")
molten.mean$id <- paste0(molten.mean$sample, "-", molten.mean$variable)

molten.sd <- melt(m4.t.sd, id.vars = "sample")
molten.sd$id <- paste0(molten.sd$sample, "-", molten.sd$variable)

# Merge the dataframes
molten <- merge(molten.mean, molten.sd, by.x = "id", by.y = "id")
#You'll get an error here, so far it's still fine and the plot still works



###Final rearragement of the dataframe
#Few last steps before we can plot the data!

molten$id <- NULL
molten$sample.y <- NULL
molten$variable.y <- NULL
names(molten) <- c("sample", "taxonomy", "mean", "sd")

molten$tax_col <- str_split_fixed(molten$taxonomy, ";", 2)[,1]
molten$tax_bin <- str_split_fixed(molten$taxonomy, ";", 2)[,2]

# Reorder the taxonomic annotation for the plot
molten <- molten[order(molten$tax_col),]
tax_levels <- as.character(molten$tax_bin[!duplicated(molten$tax_bin)])
tax_levels <- tax_levels[tax_levels != "Other"]
tax_levels <- c(tax_levels, "Other")
molten$tax_bin <- factor(molten$tax_bin, levels = rev(tax_levels))

# Reorder the samples for the plot
molten$sample <- factor(molten$sample, levels = replicates_list)

# Remove null values
molten2 <- molten[molten$mean > 0,]
colnames(molten2)[1]<- "Sample.ID"

molten3 <- merge(molten2, metadat, by = "Sample.ID")

cols2<- c("Bacteria"="#001D29", "Bacteroidota"="#005F73", 
          "Firmicutes" = "#94D2BD", "Other" = "#9B2226",
          "Proteobacteria" = "#EE9B00")

region.labs <- as_labeller(c("Belize", "Curaçao"))
###Plot 

pdf("Bubble_Genus_DCYLHealthy_15taxa_Region_filt.pdf",width=25, height=10)
bubble_plot <- ggplot(molten3,aes(Sample.ID,tax_bin)) +
  geom_point(aes(size=mean, fill=tax_col),shape=21,color="black") +
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_text(size = 24, face = "bold")) +
  ylab("Lowest Taxonomic Assignment") +
  xlab("Samples") +
  scale_fill_manual(values = cols2, name = "Phylum", labels = 
                      c("Bacteria", "Bacteroidota", "Mycoplasmatota", "Other", "Proteobacteria")) +
  scale_size(range = c(2, 20),name = "Relative abundance")+
  facet_grid(.~Region, scales = "free", labeller = labeller (region.labs))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(face = "bold"))+
  theme(strip.text = element_text(size = 24))+
  theme(axis.title.y = element_text(size = 24, face = "bold"))+
  theme(axis.title.x = element_text(size = 24, face = "bold"))+
  theme(axis.text.y = element_text(size=20, color = "black"))+
  theme(legend.title = element_text(size = 20, face = "bold"))+
  theme(legend.text = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(size=10)))
bubble_plot

dev.off()


#####################################
############Analysis of Community to detect differentially abundant taxa
#####################################

library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
source("ancom_v2.1.R")
otu <- read.table("silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps

ps<- remove_samples(c("17818", "17797"), ps)

#Remove samples that have less than 1000 reads
ps <- prune_samples(sample_sums(ps) >= 500, ps)
ps #7559 taxa across 104 samples
#Prune taxa no longer in samples
ps <- prune_taxa(taxa_sums(ps) > 0, ps) 
ps #7451 taxa in 132 samples

get_taxa_unique(ps, "Family") #544

otu<-t(otu)
metadata<-samples
metadata$Sample.ID <- row.names(metadata)
rownames(metadata) <- NULL

feature_table = otu; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
my_prepro = feature_table_pre_process(feature_table, metadata, sample_var, group_var, 
                                      out_cut, zero_cut, lib_cut, neg_lb)
feature_table = my_prepro$feature_table # my_preprocessed feature table
metadata = my_prepro$meta_data # my_preprocessed metadata
struc_zero = my_prepro$structure_zeros # Structural zero info

# Step 2: ANCOM by Region
main_var = "Region"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
t_start = Sys.time()
res = ANCOM(feature_table, metadata, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s
write.table(res$out, "ANCOM_Region.txt",sep="\t",col.names=NA)

n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))

# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")
# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.9"], label = "W[0.9]")
fig = res$fig +  
  geom_hline(yintercept = cut_off["detected_0.9"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig  

#The x-axis is the effect size, looking for ASVs that are high-W stat and high/low CLR mean difference

########################################
#Try the other ANCOM Method where you read the function in
########################################

ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
}



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected

  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}
  
#from phyloseq -use all of the data, not relative abundance, and  not filtered beyond mitochondria and choloroplasts 
#Make the otutable
#group the otus based on family
#this is probaby a little overkill, but it's how I've done it in the past - all you need is the otutable really, so probably that could have been done in one step, but i like this for future plot-making
dat <- tax_glom(ps, taxrank = "Genus") #at the Genus level

#melt the data, so it's like a dataframe
datm <- psmelt(dat)

#Cast the new datatable with columns that are of interest
datc <- reshape2::dcast(datm, Sample + Region + Reef +Season ~ Genus, value.var = 'Abundance', fun.aggregate = sum)

dim(datc) #dimensions of the table

otud <- datc[,c(1,5:1058)] #select the first column, and then all of the taxa columns 
colnames(otud)[1] <- "Sample.ID" #rename the first column to Sample.ID - this is to match ANCOM syntax

metadat <- sample_data(ps) #get the sample data
metadat <- as.data.frame(as.matrix(metadat)) #make into into a matrix

library(dplyr)
library(reshape2)
metadat <- tibble::rownames_to_column(metadat, "Sample.ID") #make sure the sample names column is called Sample.ID

names(otud) <- make.names(names(otud)) #get the names from the table for otud
otu_test <- otud #rename otud to otu_test, for syntax in ANCOM

metadat <- select(metadat, c("Sample.ID","Reef", "Region", "Season")) # use select to only use treatment columns of interest
map_test <- metadat #rename map_TEst
Vardat <- map_test #specify that this for Vardat - ANCOM syntax

### ANCOM test REGION - not adjusted, more than 2 levels = Kruskal Wallis
comparison_test_treat=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                                 Vardat=map_test, #calling the metadata
                                 adjusted=FALSE, #true if covariates are to be included for adjustment
                                 repeated=FALSE, #repeated measure
                                 main.var="Region", #main variable or factor
                                 adj.formula= NULL, #other factors to include
                                 repeat.var=FALSE, #repeated measure
                                 long = FALSE, #longitudinal study
                                 multcorr=2,
                                 sig=0.05, #significance level
                                 prev.cut=0.90) #OTUs with proportion of zeroes greater than prev.cut are not included in the analysis

res <- comparison_test_treat$W.taxa #taxa that significantly vary across factor level of interest
write.table(res,"ANCOM_genus_KruskallWallis_Region.txt",sep="\t",col.names=NA)
res2 <- res[which(res$detected_0.7==TRUE),] 

sig_region <- glue::glue_collapse(droplevels(factor(res2$otu.names)), sep = ", ") #this is to get a list of the genera that are different
print(sig_region)


datc_relabund <-  sweep(datc[,5:1014], 1, rowSums(datc[,5:1014]), '/')
datc_relnames <- cbind(datc[,1:4],datc_relabund)
#only select the significant families
sig.groups<- colnames(datc_relnames)
sig_dis <- select(datc_relnames, Sample, Region, Reef, Season, Proteobacteria, Endozoicomonas, Spiroplasma, Alphaproteobacteria, "Candidatus Amoebophilus", Maritimimonas, Fulvivirga, "JGI 0000069-P22")
sig_long <- melt(sig_dis, id.vars=c("Sample", "Region", "Reef", "Season"),variable.name="Genus",value.name="Proportion")

sum_sig <- Rmisc::summarySE(sig_long, measurevar = "Proportion", groupvars = c("Region", "Genus"), na.rm=TRUE)

pdf("ANCOM_Genera_Region.pdf",width=11)
gens_region <- ggplot(sum_sig, aes(x=Genus, y=Proportion))+
  geom_point(size=8,pch = 21, aes(fill=Region))+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(size=16))+
  theme(axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=20, face = "bold"))+
  theme(axis.title.y=element_text(size=20, face = "bold"))+
  theme(legend.position = "right")+
  geom_errorbar(aes(ymin=Proportion-se, ymax=Proportion+se), width=.1)+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.text = element_text(size=16))+
  theme(legend.title = element_text(size = 20, face="bold"))+
  scale_fill_manual(values = cols)+
  theme(axis.text.x = element_text(color = "black"))+
  theme(axis.text.y = element_text(color = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))+
  theme(axis.title.y = )
  scale_x_discrete(limits = rev(c("Alphaproteobacteria", "Candidatus Amoebophilus", "Maritimimonas",
                                  "Proteobacteria", "Endozoicomonas", "Spiroplasma", "Fulvivirga", "JGI 0000069-P22")))
gens_region
dev.off()


##################################
#Prevalence 
#################################

ps

ps.prev.75<- phyloseq_filter_prevalence(ps, prev.trh = 0.75, abund.trh = NULL)  
ps.prev.75
colnames(ps.prev.75@otu_table) #prints out ASVs, You can BLAST this

ps.prev.70<- phyloseq_filter_prevalence(ps, prev.trh = 0.7, abund.trh = NULL)  
ps.prev.70
colnames(ps.prev.70@otu_table) #prints out ASVs, You can BLAST this
write.table(ps.prev.70@otu_table, "ps_prev_70_otu.txt", sep="\t")
write.table(ps.prev.70@tax_table, "ps_prev_70_taxa.txt", sep="\t")


######Look at core microbiome by region
ps.belize<- subset_samples(ps, Region=="Belize")
ps.belize

ps.curacao<- subset_samples(ps, Region=="Curacao")
ps.curacao

###Find ASVs present in over 70% of belize samples
ps.belize.prev.70<- phyloseq_filter_prevalence(ps.belize, prev.trh = 0.7, abund.trh = NULL) 
ps.belize.prev.70
colnames(ps.belize.prev.70@otu_table)

########Find ASVs present in over 70% of curacao samples
ps.curacao.prev.70<- phyloseq_filter_prevalence(ps.curacao, prev.trh = 0.7, abund.trh = NULL) 
ps.curacao.prev.70
colnames(ps.curacao.prev.70@otu_table)



#############################################
#Plot alpha diversity metrics
#############################################

ps
richness<- plot_richness(ps, x="Region", measures=c("Chao1", "Shannon"), title = "Alpha Diversity")


Regionlabs <- c("Belize", "Curaçao")

pdf("Alpha_Div_Boxplot_DCYL_filt.pdf",width=13)
richness<- plot_alpha(ps, x="Region", measures=c("Chao1", "Shannon"), title = "Alpha Diversity", estimators = NULL)
richness$layers <- richness$layers[-1]
richness + geom_boxplot(aes(fill=Region), alpha = 0.4,  outlier.shape=NA)+
  geom_jitter(position=position_jitter(width=.1,height=0),aes(colour="black",fill=Region),size=3, pch=21, colour = "black")+
  scale_fill_manual(values = cols, labels = c('Belize', 'Curaçao'))+
  scale_color_manual(values = c("#003459", "black")) +
  theme(strip.background =element_rect(fill="white"))+
  theme(axis.text.x = element_text(colour="black"))+
  theme(axis.text.y = element_text(colour="black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1.1, vjust=1.1))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  theme(strip.text.x=element_text(size=16, face="bold"))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(size = 12))+
  guides(color="none")+
  scale_x_discrete(labels=Regionlabs)
dev.off()


##########Statistical analysis of shannon diversity

data_adiv<- estimate_richness(ps, measures  = c("Chao1", "Shannon", "InvSimpson"))
shapiro.test(data_adiv$Shannon)
shapiro.test(data_adiv$Chao1)
shapiro.test(data_adiv$InvSimpson)

rownames(data_adiv) <- sub(".*X", "", rownames(data_adiv))
data_adiv <- tibble::rownames_to_column(data_adiv, "Sample.ID") 
data_adiv<- merge(data_adiv, metadat, by.x = "Sample.ID")

kruskal.test(Shannon ~ Region, data=data_adiv) 
kruskal.test(Chao1 ~ Region, data=data_adiv)
kruskal.test(InvSimpson ~ Region, data=data_adiv)

ps
length(get_taxa_unique(ps, taxonomic.rank = "Genus"))


ps.curacao
write.table(ps.curacao@otu_table, "ps.curacao_otu_table.txt", sep="\t", col.names = NA)

ps.belize
write.table(ps.belize@otu_table, "ps.belize_otu_table.txt", sep="\t", col.names = NA)

ps5_ra.belize<- subset_samples(ps5_ra, Region=="Belize")
write.table(ps5_ra.belize@otu_table, "ps5_ra.belize_otu_table.txt", sep="\t", col.names = NA)
summarize_phyloseq(ps)


#################################################
#Endozoicomonas AVSs
#################################################
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



endo<- subset_taxa(ps5_ra, Genus=="Endozoicomonas")
endo

otu.endo<- as(otu_table(endo), "matrix")
taxon.endo<- as(tax_table(endo), "matrix")
meta.endo<- as(sample_data(endo), "matrix")
otu.endo<- as.data.frame(otu.endo)
otu.endo<- rownames_to_column(otu.endo, var="Sample")

#Export Endozoicomonas ASVS
write.table(otu.endo, "Endozoicomonas_ASVs.txt", sep="\t", col.names=NA)
write.csv(otu.endo, "Endozoicomonas_ASVs.csv")

#Make ASv names
names(otu.endo)[2] <- "ASV1"
names(otu.endo)[3] <- "ASV2"
names(otu.endo)[4] <- "ASV3"
names(otu.endo)[5] <- "ASV4"
names(otu.endo)[6] <- "ASV5"
names(otu.endo)[7] <- "ASV6"
names(otu.endo)[8] <- "ASV7"


meta.endo<- as.data.frame(meta.endo)
meta.endo<- rownames_to_column(meta.endo, var="Sample")
otu.endo.meta<- merge(meta.endo, otu.endo, "Sample")
otu.endo.long<- melt(otu.endo.meta, id.vars=c("Sample", "Reef", "Region", "Collection.Date", "Season"), variable.name = "ASV", value.name="Proportion" )
write.csv(otu.endo.long, "otu_endo_long.csv")

otu.endo.long<- otu.endo.long[!(otu.endo.long$Sample == "17852"),]
dim(otu.endo.long)

library(randomcoloR)
n <- 7
palette <- distinctColorPalette(n)
pdf("Endozoicomonas_bars_2.pdf", width=11)
p1<- ggplot(otu.endo.long, aes(x=Proportion, y=Reef))+
  geom_bar(aes(fill=ASV), stat="summary", position="stack", colour="black")+
  scale_fill_manual(values=palette, name="Endozoicomonas ASV")+
  facet_grid(Region~., scales="free", space="free")+
  theme_classic() +
  theme(strip.text = element_text(face="bold", size = 22))+
  theme(strip.background = element_rect( fill="white"))+
  theme(legend.title=element_blank()) +
  theme(text=element_text(size=14))+
  theme(axis.text.y = element_text(color="black", size = 16))+
  theme(axis.text.x = element_text(color="black", size=16))+
  theme(axis.title.x = element_text(size = 20, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  labs(fill='Endozoicomonas ASV') +
  xlab('Relative Abundance')

p1
dev.off()

##################################################
#Candidatus ASVs
##################################################
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

candi<- subset_taxa(ps5_ra, Genus=="Candidatus Amoebophilus")
candi

otu.candi<- as(otu_table(candi), "matrix")
taxon.candi<- as(tax_table(candi), "matrix")
meta.candi<- as(sample_data(candi), "matrix")
otu.candi<- as.data.frame(otu.candi)
otu.candi<- rownames_to_column(otu.candi, var="Sample")

#Export Endozoicomonas ASVS
write.table(otu.candi, "Candidatus_ASVs.txt", sep="\t", col.names=NA)
write.csv(otu.candi, "Candidatus_ASVs.csv")

#Make ASv names
names(otu.candi)[2] <- "ASV1"
names(otu.candi)[3] <- "ASV2"
names(otu.candi)[4] <- "ASV3"
names(otu.candi)[5] <- "ASV4"
names(otu.candi)[6] <- "ASV5"


meta.candi<- as.data.frame(meta.candi)
meta.candi<- rownames_to_column(meta.candi, var="Sample")
otu.candi.meta<- merge(meta.candi, otu.candi, "Sample")
otu.candi.long<- melt(otu.candi.meta, id.vars=c("Sample", "Reef", "Region", "Collection.Date", "Season"), variable.name = "ASV", value.name="Proportion" )
write.csv(otu.candi.long, "otu_candi_long.csv")

otu.candi.long<- otu.candi.long[!(otu.candi.long$Sample == "17852"),]
dim(otu.candi.long)

library(randomcoloR)
n <- 5
palette <- distinctColorPalette(n)
pdf("Candidatus_bars_2.pdf", width=11)
p1<- ggplot(otu.candi.long, aes(x=Proportion, y=Reef))+
  geom_bar(aes(fill=ASV), stat="summary", position="stack", colour="black")+
  scale_fill_manual(values=palette, name="Candidatus ASV")+
  facet_grid(Region~., scales="free", space="free")+
  theme_classic() +
  theme(strip.text = element_text(face="bold", size = 22))+
  theme(strip.background = element_rect( fill="white"))+
  theme(legend.title=element_blank()) +
  theme(text=element_text(size=14))+
  theme(axis.text.y = element_text(color="black", size = 16))+
  theme(axis.text.x = element_text(color="black", size=16))+
  theme(axis.title.x = element_text(size = 20, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  labs(fill='Candidatus ASV') +
  xlab('Relative Abundance')

p1
dev.off()


#################################################
#Spiroplasma AVSs
#################################################

spiro<- subset_taxa(ps5_ra, Genus=="Spiroplasma")
spiro

otu.spiro<- as(otu_table(spiro), "matrix")
taxon.spiro<- as(tax_table(spiro), "matrix")
meta.spiro<- as(sample_data(spiro), "matrix")
otu.spiro<- as.data.frame(otu.spiro)
otu.spiro<- rownames_to_column(otu.spiro, var="Sample")

#Export Endozoicomonas ASVS
write.table(otu.spiro, "Spiroplasma_ASVs.txt", sep="\t", col.names=NA)
write.csv(otu.spiro, "Spiroplasma_ASVs.csv")

#Make ASv names
names(otu.spiro)[2] <- "ASV1"
names(otu.spiro)[3] <- "ASV2"


meta.spiro<- as.data.frame(meta.spiro)
meta.spiro<- rownames_to_column(meta.spiro, var="Sample")
otu.spiro.meta<- merge(meta.spiro, otu.spiro, "Sample")
otu_spiro_long<- melt(otu.spiro.meta, id.vars=c("Sample", "Reef", "Region", "Collection.Date", "Season"), variable.name = "ASV", value.name="Proportion" )

pdf("Spiroplasma_bars.pdf", width=11)
p1<- ggplot(otu_spiro_long, aes(x=Sample, y=Proportion))+
  geom_bar(aes(fill=ASV), stat="identity", position="stack")+
  facet_grid(.~Reef, scales="free", space="free")+
  theme(strip.text = element_text(face="bold"))+
  theme(axis.text.x=element_text(angle=90, size=7))+
  theme(legend.title=element_blank()) +
  theme(text=element_text(size=14))
p1
dev.off()



###Box and whiskers
# Group the ASVs based on genus
dat <- tax_glom(ps, taxrank = "Genus")

# Melt the data, so it's like a dataframe
datm <- psmelt(dat)

# Cast the new datatable with columns that are of interest. This will need to be modified for your project. In the example below the metadata columns are coral, site, susceptibility. For Dendrogyra, you will just have country and maybe reef

datc <- data.table::dcast(datm, Sample + Region + Reef  ~ Genus, value.var = 'Abundance', fun.aggregate = sum)

# Get dimensions of the table so you have the number of columns for your dataset
dim(datc)

#Calculate relative abundance; remember to adjust number of columns for your dataset. In the example below there are 224 columns total and the first 4 columns are metadata

datc_relabund <-  sweep(datc[,4:1013], 1, rowSums(datc[,4:1013]), '/')
datc_relnames <- cbind(datc[,1:3],datc_relabund)


#select Endozoicomonas (or whatever genus you want); repeat this section for Spiroplasma
spiro <- select(datc_relnames, Sample, Region, Reef, Spiroplasma)

#change column name from Spiroplasma to Proportion
names(spiro)[4] <- 'Proportion'

pdf("Spiroplasma_byRegion.pdf",width=8.5)
p <- ggplot(spiro, aes(x=Region, y=Proportion))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(position=position_jitter(width=.1,height=0),aes(color=Region),size=3)+
  theme_bw()+
  coord_flip()+
  scale_color_manual(values=cols)+
  theme(axis.title.x=element_blank())+
  theme(text=element_text(size=14))+
  theme(strip.text.y=element_text(size=14))+
  ylab("Relative Abundance of Spiroplasma")
  p
dev.off()

otu.spiro.long.asv.1.2<- otu_spiro_long[otu_spiro_long$ASV=="ASV1" | otu_spiro_long$ASV=="ASV2", ]


pdf("Spiroplasma_ASVs_1and2_Box_Whiskers_2.pdf",width=8.5)
p <- ggplot(otu.spiro.long.asv.1.2, aes(x=ASV, y=Proportion))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(position=position_jitter(width=.1,height=0), aes(fill=Region),size=4, pch=21)+
  theme_bw()+
  scale_fill_manual(values=cols)+
  theme(axis.title.x=element_blank())+
  theme(text=element_text(size=14))+
  theme(strip.text.y=element_text(size=14))+
  ylab("Relative Abundance of Spiroplasma")+
  theme(axis.text.x=element_text(color="black", size = 16))+
  theme(axis.text.y=element_text(color="black", size = 16))+
  theme(axis.title.x = element_text(face="bold", size = 20))+
  theme(axis.title.y = element_text(face="bold", size = 20))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(size = 12))
p
dev.off()


####################################
#Core Microbiome
####################################


library(microbiome)
library(Rmisc)

ps5_ra.genus<- tax_glom(ps5_ra, taxrank = "Genus")

head(prevalence(ps5_ra, detection = 00/100, sort = TRUE))
pseq.core <- core(ps5_ra.genus, detection = 0, prevalence = 0.75)


all.core.melt <- psmelt(pseq.core)
unique(all.core.melt$OTU)

### Means and SE for Core microbiomes ####
all.core.melt_prev75_SE <- summarySE(all.core.melt, measurevar = "Abundance", groupvars = "Genus")
all.core.melt_prev75_SE

write.table(all.core.melt_prev75_SE, "core_microbiome_75.txt")

#########Core Microbiome by ASV rather than Genus
head(prevalence(ps_ra, detection = 00/100, sort = TRUE))
pseq.core <- core(ps_ra, detection = 0, prevalence = 0.75)


all.core.melt <- psmelt(pseq.core)
unique(all.core.melt$OTU)

### Means and SE for Core microbiomes ####
all.core.melt_prev75_SE <- summarySE(all.core.melt, measurevar = "Abundance", groupvars = "Genus")
all.core.melt_prev75_SE

write.table(all.core.melt_prev75_SE, "core_microbiome_75.txt")


##########################################
#Try Different bubble plot
##########################################

##find relative abundance
#Pick only the top 15 most abundant taxa
genus.sum = tapply(taxa_sums(ps5_ra), tax_table(ps5_ra)[, "Genus"], sum, na.rm=TRUE)
top15genus = names(sort(genus.sum, TRUE))[1:15]
ps_ra15 = prune_taxa((tax_table(ps5_ra)[, "Genus"] %in% top15genus), ps5_ra)

ps_ra15 #taxa 51
get_taxa_unique(ps_ra15, "Genus")

###use a nested for-loop to pull relative abundance data out of phyloseq object (by taxonomic givision)
graphgenus <- unique(tax_table(ps_ra15)[,6])
graphgenus <- as.vector(graphgenus)
genus.sums.m <- matrix(nrow=132, ncol=15) # Rows is number of samples
rownames(genus.sums.m) <- sample_names(ps_ra15) # specify row names
colnames(genus.sums.m) <- graphgenus # specify column names

for (i in sample_names(ps_ra15)) { 
  for (j in graphgenus) {
    genus.sums.m[[i,j]] <- sum(otu_table(ps_ra15)[sample_names(ps_ra15)==i,which(tax_table(ps_ra15)[,6]==j)]) # this is calculating the average percent abundance of each genus
  }
}

###Make data frames (in long format) from the matix of abundances to use with ggplot
flux.genus.wide <- t(genus.sums.m)
taxonrep <- rep(rownames(flux.genus.wide), each=132)
tax <- as.character(taxonrep)
grouprep <- colnames(flux.genus.wide)
SampleID <- rep(grouprep,times=15)
abund <- as.numeric(t(flux.genus.wide))
longdf <- cbind.data.frame(tax,abund,SampleID)


#Make sample data table
metadat_rel <- as.data.frame(as.matrix(sample_data(ps5_ra)))
metadat_rel

metadat.rel <- cbind(rownames(metadat_rel), metadat_rel)
rownames(metadat.rel) <- NULL
colnames(metadat.rel) <- c("SampleID", "Reef", "Region","Collection.Date", "Season")
head(metadat.rel)

#merge longdf with metadata table by SampleID
long.metadata_df <- merge(metadat.rel, longdf, by="SampleID")
colnames(long.metadata_df)[6]<- "Genus"

#Get Phyla for each genus
ps_ra15_taxtable<- data.frame(tax_table(ps_ra15))
long.metadata_df<- merge(long.metadata_df, ps_ra15_taxtable, by="Genus")

#create a table
write.csv(long.metadata_df, "Top15RelativeAbundanceTable.csv")

#Assign colors per phylum
length(unique(long.metadata_df$Phylum))
unique(long.metadata_df$Phylum)

cols <- c("Bacteroidota"="#2A9D8F","Bacteria"="#264653","Proteobacteria"="#E9C46A",
          "Firmicutes"="#F4A261")
library(stringr)

pdf("Bubbleplot_Top15.pdf",width=20,height=8.5)
bubplot <- ggplot(long.metadata_df, aes(x=SampleID, y=Genus, size=abund, fill=Phylum))+
  geom_point(shape=21)+
  scale_size(range = c(0, 17), name="Relative Abundance")+
  theme(axis.text.x = element_blank())+
  facet_grid(.~Region, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = cols, labels = c(str_wrap("Unclassified Bacterial Phylum", width=20), "Bacteroidota", "Firmicutes", "Proteobacteria")) +
  theme(strip.background = element_blank())+
  xlab("Sample") +
  ylab("Lowest Taxonomic Assignment") +
  scale_y_discrete(limits= c("Vibrio", "Tistlia", "Thalassobaculales", "Ruegeria",
                             "Proteobacteria", "Nisaea", "Dstr-E11", "Spiroplasma",
                             "Fulvivirga", "Flavobacteriaceae", "Cyclobacteriaceae",
                             "Candidatus Amoebophilus", "Amoebophilaceae", "Bacteria"),
                    labels= c(expression(italic("Vibrio")), expression(italic("Tistlia")), 
                                              expression(italic("Thalassobaculales")), 
                                              expression(italic("Ruegeria")), "Unclassified Proteobacteria",
                                              expression(italic("Nisaea")), "Dstr-E11", expression(italic("Spiroplasma")), 
                                              expression(italic("Fulvivirga")), "Flavobacteriaceae",
                                              "Cyclobacteriaceae", expression(italic("Candidatus") * " Amoebophilus"), 
                                              "Amoebophilaceae", "Unclassified Bacterial Phylum"))+
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_text(size = 22, face = "bold")) +
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(face = "bold"))+
  theme(strip.text = element_text(size = 22))+
  theme(axis.title.y = element_text(size = 22, face = "bold"))+
  theme(axis.title.x = element_text(size = 22, face = "bold"))+
  theme(axis.text.y = element_text(size=18, color = "black"))+
  theme(legend.title = element_text(size = 18, face = "bold"))+
  theme(legend.text = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(size=7.5, vjust=30.5)))

bubplot
dev.off()


##############################################
## DESeq2 analysis
##############################################
# check to see that you are using the full dataset (ps), not the low abundance ASV filtered dataset (ps5)
ps

# Define the order of the conditions for testing
# In this order, the positive fold change values are increased in Belize
sample_data(ps)$Region<-factor(sample_data(ps)$Region,levels=c("Curaçao","Belize"))
head(sample_data(ps)$Region, 10)

# DESEQ2 analysis on Dcyl Samples
dds = phyloseq_to_deseq2(ps, ~ Region)
dds
#filter rows with very few counts
dds <- dds[ rowSums(counts(dds)) > 5, ]
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- DESeq(dds,test="Wald", fitType="parametric")
res = results(dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
#save table of results
sig = as(sigtab, "matrix")
write.table(sig,"DESeq2_results.txt",sep="\t",col.names=NA)
sig_dif<- read.table("DESeq2_results.txt", header = TRUE)
colnames(sig_dif)[1]<- "ASV"

#Count how many phyla
phyla<-unique(sig_dif$Phylum)
print(phyla)
cols <- c("Bacteroidota"="#2A9D8F","Bacteria"="#264653","Proteobacteria"="#E9C46A",
          "Firmicutes"="#F4A261")


#Create Figure for Significantly Different Taxa
pdf("FigureDEtaxa.pdf",width=11,height=8)
p1 <- ggplot(sig_dif,aes(log2FoldChange,Genus,fill=Phylum))+
  geom_point(size=5, shape=21)+
  geom_errorbarh(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE),height=0.2)+
  scale_fill_manual(values=cols, labels = c("Unclassified Bacterial Phylum", "Bacteroidota", "Firmicutes", "Proteobacteria"))+
  theme(plot.title = element_text(face="italic"))+
  labs(y= "Lowest Taxonomic Assignment")+
  theme(axis.title=element_text(size=14))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=16))+
  theme(axis.text=element_text(size=14, color="black"))+
  theme(legend.position = "right")+
  scale_y_discrete(limits=rev,
                   labels = c(expression(italic("Spiroplasma")), "Unclassified Proteobacteria",
                              expression(italic("Fulvivirga")), "Flavobacteriaceae",
                              expression(italic("Enhydrobacter")), expression(italic(
                                "Endozoicomonas")), expression(italic("Candidatus") * " Amoebophilus"),
                              "Unclassified Bacterial Phylum"))+
  geom_vline(xintercept=0, linetype="dashed")
p1
dev.off()


###############################################
#Are there Spiroplasma ASVs in other cnidarian microbiomes?
###############################################
#Read in McCauley Phyloseq Object
v4<- readRDS("V4 Library.rdata")
v4

#Have to subset the full V4 library by those that only use the same primers
v4<- subset_samples(v4, Primer_Region=="515F-806R")
v4
#Prune taxa no longer in samples
v4 <- prune_taxa(taxa_sums(v4) > 0, v4) 
v4 

mccauley_spiro<- subset_taxa(v4, Genus == "Spiroplasma")
mccauley_spiro

#Remove samples that have less than 500 reads
mccauley_spiro <- prune_samples(sample_sums(mccauley_spiro) > 0, mccauley_spiro)
mccauley_spiro 

View(mccauley_spiro)
