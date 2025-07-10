#Load packages
library(dada2)

#Import sequence data
path <- "data"
list.files(path)

#Sort forward and reverse files
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

#Plot quality profiles
plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])

#Extract sample IDs
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Create a subdirectory for filtered files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#Filter and trim poor reads, PhiX and ambiguous bases
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 160), 
                     maxEE = c(2,2), truncQ = 2, maxN = 0, rm.phix = TRUE, 
                     compress = TRUE, multithread = TRUE)
View(out)

#Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

#Dereplication
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

#Name the dereplicated class objects by the sample IDs
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample inferencing
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
dadaFs[[1]]
dadaRs[[2]]

#Merge the denoised forward and reverse files
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
head(mergers[[1]])

#Generate sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                    multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v128_train_set.fa.gz", multithread = TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v128.fa.gz")

#Evaluate pipeline accuracy
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.mock), 
    "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

#Create a simple metadata dataframe
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
write.csv(samdf, "data/samdf.csv")

#Save all outputs
saveRDS(seqtab.nochim, "output/seqtab.nochim.rds")
saveRDS(taxa, "output/taxa.rds")



