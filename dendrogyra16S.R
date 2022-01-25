#DADA2 Dendrogyra Batch 1 Run 3 Processing
install.packages("devtools")
devtools::install_version("Matrix", version = "1.3.2")
BiocManager::install("Matrix")
library(Matrix)
BiocManager::install("ShortRead")

#Install and load package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")
library(dada2); packageVersion("dada2")

install.packages("R.utils")
library(R.utils)

remove.packages("Matrix")

path <- "~/Documents/Dendrogyra_Batch1/cutadapt" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq.gz", full.names = TRUE))

fnFs

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
length(sample.names)
sample.names

library(strex)

forward.names<- basename(fnFs)
forward.names

reverse.names<- basename(fnRs)
reverse.names

write.csv(sample.names, file="sample_names.csv", row.names = F)
write.csv(forward.names, file="forward_names.csv", row.names = F)
write.csv(reverse.names, file="reverse_names.csv", row.names = F)

#Plot quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

write.csv(sample.names, file="sample_names.csv", row.names=F)

#Filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#Try Looking at Error Rates
errF <- learnErrors(fnFs, multithread=TRUE)
errR<- learnErrors(fnRs, multithread = TRUE)

plotErrors(errF, nominalQ=TRUE)

#Try Sequence Algorithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
