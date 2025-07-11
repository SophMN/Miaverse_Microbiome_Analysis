#Load packages
library(mia)
library(miaViz)
library(patchwork)
library(scater)

#Import data
samdf <- read.csv("data/samdf.csv")
seqtab.nochim <- readRDS("output/seqtab.rds")
taxa <- readRDS("output/taxa.rds")
tse <- readRDS("output/tse.rds")

###Creating a TreeSummarisedExperiment object
#Transpose the ASV count table
samples.out <- rownames(seqtab.nochim)
counts <- t(seqtab.nochim)
counts <- as.matrix(counts)
assays <- SimpleList(counts = counts)

#Convert colData and rowData into DataFrame
samdf <- DataFrame(samdf)
taxa <- DataFrame(taxa)

#Create TSE
tse <- TreeSummarizedExperiment(assays = assays, 
                                colData = samdf, 
                                rowData = taxa)
tse
colData(tse)
dim(colData(tse))
rowData(tse)
dim(rowData(tse))
assays(tse)

#Convert TSE into a phyloseq object
MiSeqSop_ps <- convertToPhyloseq(tse)
MiSeqSop_ps
saveRDS(MiSeqSop_ps, "output/miseqsop_ps.rds") 

#Remove mock community sample
tse <- tse[, colnames(tse) != "Mock"]
tse

#Extract reference sequences
dna <- Biostrings::DNAStringSet(rownames(tse))

#Add reference sequences to the RefSeq slot
referenceSeq(tse) <- dna
tse

#Convert rownames into ASV_number format
rownames(tse) <- paste0("ASV", "", seq(nrow(tse)))
tse

#Check the number of kingdoms in the rowData
unique(rowData(tse)$Kingdom)

#Assess the distribution of different phyla and genera
table(rowData(tse)$Phylum)
table(rowData(tse)$Genus)

##Agglomerate at the genus and species levels
#Genus
altExp(tse, "Genus") <- agglomerateByRank(tse, rank = "Genus", na.rm = TRUE)
tse

#Species
altExp(tse, "Species") <- agglomerateByRank(tse, rank = "Species", na.rm = TRUE)
tse

#Melt the agglomerated species TSE into a dataframe
tse_species <- altExp(tse, "Species")
species_df <- meltSE(tse_species, add_row_data = TRUE, 
                     add_col_data = TRUE, assay.type = "counts")
View(species_df)

#Export
write.csv(species_df, "output/species_df.csv", row.names = FALSE)

##Visualise the relative abundance of the top genera
tse_genus <- altExp(tse, "Genus")
tse_genus

#Transform counts into relative abundances
tse_genus <- transformAssay(tse_genus, assay.type = "counts", 
                            method = "relabundance")
assays(tse_genus)

#Get the top 20 genera
top_genera <- getTop(tse_genus, top = 20)
tse_genus_top <- tse_genus[rowData(tse_genus)$Genus %in% top_genera, ]
tse_genus_top

#Plot barplot
#Rename column name X appropriately
colnames(colData(tse_genus_top))[colnames(colData(tse_genus_top)) == "X"] <- "SampleID"

#Plot
plot <- plotAbundance(tse_genus_top, assay.type = "relabundance", 
              group = "Genus", order.row.by = "abund", col.var = "When", 
              facet.cols = TRUE, scales = "free_x")
plot

##Calculating beta diversity
#Convert counts to relative abundance
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
assays(tse)

#Calculate Bray-Curtis dissimilarity distances
tse <- runNMDS(tse, FUN = getDissimilarity, method = "bray", 
              assay.type = "relabundance", name = "NMDS", ncomponents = 10)
tse

#Generate NMDS plot
plotReducedDim(tse, "NMDS", color_by = "When")

#Save output
saveRDS(tse, "output/tse.rds")
saveRDS(tse_genus_top, "output/tse_genus_top.rds")

#Exploratory data analysis and QC with PCA
tse <- readRDS("output/tse.rds")

#Centre-log ratio transformation
tse <- transformAssay(x = tse, assay.type = "counts", 
                      method = "clr", pseudocount = TRUE, 
                      name = "clr")
tse

#Run PCA
tse <- runPCA(tse, ncomponents = 10, assay.type = "clr")
tse

#Plot PCA
plotReducedDim(tse, "PCA", color_by = "When")




