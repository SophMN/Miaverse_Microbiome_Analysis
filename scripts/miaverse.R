#Load packages
library(mia)
library(miaViz)
library(patchwork)
library(scater)
library(knitr)
library(ggplot2)
library(microbiomeDataSets)

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

##Taxonomic information
#Accessing taxonomic information
checkTaxonomy(tse)

#Check the taxonomy ranks in rowData
taxonomyRanks(tse)

#Subset rowData to the needed columns
rowData(tse)[, taxonomyRanks(tse)]

#Check for empty values in a selected taxonomy rank
all(!taxonomyRankEmpty(tse, rank = "Kingdom"))
all(!taxonomyRankEmpty(tse, rank = "Genus"))
all(!taxonomyRankEmpty(tse, rank = "Species"))
table(taxonomyRankEmpty(tse, rank = "Genus"))
table(taxonomyRankEmpty(tse, rank = "Species"))

#Get taxonomy labels
getTaxonomyLabels(tse) |> head()

#Get information on certain taxa
getUnique(tse, rank = "Phylum")
getUnique(tse, rank = "Genus") |> head()
getUnique(tse, rank = "Species") |> head()

#Get information that matches a specific taxa
mapTaxonomy(tse, taxa = "Lactobacillus")
mapTaxonomy(tse, taxa = "Helicobacter")

##Data wrangling
#Split the data based on a specific variable
tse_genus <- altExp(tse, "Genus")
tse_list <- splitOn(tse_genus, group = "When")
tse_list

#Get the top genera in the dataset
top_taxa <- sapply(tse_list, getTop)
top_taxa
top_taxa |> kable()

###Exercise
#Load data
data("GlobalPatterns", package = "mia")
tse2 <- GlobalPatterns
tse2
saveRDS(tse2, "output/tse2.rds")

#Explore the sample metadata
colData(tse2)

#Add arbitrary groups to colData
colData(tse2)$Group <- sample(c("group1", "group2"), size = ncol(tse2), replace = TRUE)
colData(tse2)

##Explore the group distribution with a barplot
#First melt the TSE into a data frame for easier plotting with ggplot
metadata <- meltSE(tse2, assay.type = "counts", add_row_data = TRUE, 
                     add_col_data = TRUE, row.names = FALSE)
#Plot
ggplot(metadata, aes(x = Group, fill = Group)) +
  geom_bar(color = "black") +
  theme_bw()

#Split the data based on groups
tse2_group <- splitOn(tse2, group = "Group")
tse2_group

#Merge the dataset you're working on with another one
tse3 <- mergeSEs(list(tse, tse2), join = "inner")
tse3

###Data exploration and quality control with baboon gut microbiome data
#Summarise data
tse_baboon <- baboongut()
tse_baboon

#Calculate summary tables to examine the variance in library sizes
summary(tse_baboon, assay.type = "counts")

#Evaluate the absolute and relative abundance of each taxon at the genus level
df <- summarizeDominance(tse_baboon, rank = "Genus")
df

#Get unique taxa 
getUnique(tse_baboon, rank = "Phylum")

#Get the top 10 taxa
top_features <- getTop(tse_baboon, method = "median", top = 10)
top_features

#Get the prevalent taxa: retrieves taxa which exceed the specified prevalence and detection thresholds
prev <- getPrevalent(tse_baboon, rank = "Genus", prevalence = 0.2, detection = 0)
prev |> head()

#Get the rare taxa
rare <- getRare(tse_baboon, rank = "Genus", prevalence = 0.2, detection = 0)
rare |> head()

##Library size
#Calculate and add total counts to colData
tse_baboon <- addPerCellQC(tse_baboon)
colData(tse_baboon)

#Visualise the lib sizes with a violin plot and histogram
plot1 <- plotColData(tse_baboon, x = "sex", y = "total", colour_by = "age")
plot2 <- plotHistogram(tse_baboon, col.var = "total")
plot1+plot2


##Contamination
colData(tse_baboon) |> head() |> kable()

#Detect contaminants with the frequency-based method
tse_baboon <- addContaminantQC(tse_baboon, concentration = "post_pcr_dna_ng")
rowData(tse) |> head()
saveRDS(tse_baboon, "output/tse_baboon.rds")
