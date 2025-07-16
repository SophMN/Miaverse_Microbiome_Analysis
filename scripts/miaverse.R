#Load packages
library(mia)
library(miaViz)
library(patchwork)
library(scater)
library(knitr)
library(dplyr)
library(ggplot2)
library(microbiomeDataSets)

#Import data
samdf <- read.csv("data/samdf.csv")
seqtab.nochim <- readRDS("output/seqtab.rds")
taxa <- readRDS("output/taxa.rds")
tse <- readRDS("output/tse.rds")
tse_baboon <- readRDS("output/tse_baboon.rds")

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

#Evaluating distribution of microbiome count data 
plotHistogram(tse_baboon, assay.type = "counts")

#Evaluating the distribution of continuous variables in colData
colData(tse_baboon)
plotHistogram(tse_baboon, col.var = "age")
summary(colData(tse_baboon)$age)
plotHistogram(tse_baboon, col.var = "asv_richness")
summary(colData(tse_baboon)$asv_richness)

#Evaluating distribution of categorical variables in colData
plotBarplot(tse_baboon, col.var = "sex")

#Evaluating abundance
#Transform the counts into relative abundance
tse_baboon <- transformAssay(tse_baboon, method = "relabundance")
assays(tse_baboon)

#Plot relative abundance
#Jitter plot
plotAbundanceDensity(tse_baboon, layout = "jitter", assay.type = "relabundance", 
                     n = 40, point.size = 1, point.shape = 19, point.alpha = 0.1) +
  scale_x_log10(label = scales::percent)

#Density plot of the top 5 taxa by sex
plotAbundanceDensity(tse_baboon, layout = "density", assay.type = "relabundance", 
                     n = 5, colour.by = "sex", point.alpha = 0.1) +
  scale_x_log10()

#Violin plot
top <- getTop(tse_baboon, top = 10L, method = "mean")
plotExpression(tse_baboon, features = top, x = "sex", 
               assay.type = "relabundance", point_alpha = 0.01) +
  scale_y_log10()

##Evaluate prevalence at a threshold of 0.1%
#Add prevalence
tse_baboon <- addPrevalence(tse_baboon, detection = 0.1/100, as.relative = TRUE)

#Plot histogram 
plotHistogram(tse_baboon, row.var = "prevalence")
rowData(tse_baboon)

#Plot heatmap
plotRowPrevalence(tse_baboon, as.relative = TRUE) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

#Relationship between prevalence and abundance
plotPrevalence(tse_baboon, as.relative = TRUE)

#Calculate the proportion of core taxa
tse_baboon <- addPrevalentAbundance(tse_baboon, prevalence = 50/100, detection = 0.1/100)
colData(tse_baboon)

#Visualise
plotHistogram(tse_baboon, col.var = "prevalent_abundance")
saveRDS(tse_baboon, "output/tse_baboon.rds")

###Exercise: QC and data exploration
#Load mouse gut data
tse <- readRDS("output/tse.rds")
tse

#Summarise the counts with a histogram and violin plot
summary(tse, assay.type = "counts")
plotHistogram(tse, assay.type = "counts")

#Add library sizes and visualise the library size distribution with a histogram
tse <- addPerCellQC(tse)
colData(tse)
plotHistogram(tse, col.var = "total") #The sampling depth differs 
plotColData(tse, x = "When", y = "total", colour_by = "Day")

#Add prevalence of the taxa to rowData
tse <- addPrevalence(tse, detection = 0.1/100, as.relative = TRUE)
rowData(tse)

#Visualise prevalence distribution with a histogram
plotHistogram(tse, row.var = "prevalence") #The data includes a few prevalent taxa

#Visualise categorical variables in colData with a barplot
plotBarplot(tse, col.var = "When")

#Get the available taxonomy ranks in the data
taxonomyRanks(tse)

#Calculate a table that summarises the dominance of genera
df <- summarizeDominance(tse, rank = "Genus")
View(df)

#Get the most prevalent features in a specific taxonomic rank using the counts table
#Set prevalence to 20% and detection threshold to 1
prev <- getPrevalent(tse, assay.type = "counts", rank = "Genus", 
                     prevalence = 0.2, detection = 1)
prev |> head()

#Get the most abundant features by median abundance
top_taxa <- getTop(tse, rank = "Genus", method = "median", top = 10)
rowData(tse)

#Visualise the prevalent features
plotRowPrevalence(tse, as.relative = TRUE) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
saveRDS(tse, "output/tse.rds")

##Data subsetting using the GlobalPatterns dataset
#Load TSE
tse2 <- readRDS("output/tse2.rds")
tse2
dim(tse2)

#Subset by sample: column-wise
tse2$SampleType |> table() |> kable()

#Subset by sample type
tse2_sub <- tse2[ , tse2$SampleType %in% c("Feces", "Skin", "Tongue")]
dim(tse2_sub)

#Subset by feature: row-wise
rowData(tse2)$Phylum %>%  table %>%  head() %>%  
  knitr::kable()
selected <- rowData(tse2)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse2)$Phylum)
tse2_sub <- tse2[selected, ]
dim(tse2_sub)

#Subset by samples and features
selected_rows <- rowData(tse2)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse2)$Phylum)
selected_cols <- tse2$SampleType %in% c("Feces", "Skin", "Tongue")
tse2_sub <- tse2[selected_rows, selected_cols]
dim(tse2_sub)

#Filtering based on library size
ind <- rowData(tse2)$Species == "Achromatiumoxaliferum"
ind[is.na(ind)] <- FALSE
tse2_sub <- tse2[ind, ]
tse2_sub

#Calculate the number of counts in each sample
tse2_sub <- addPerCellQCMetrics(tse2_sub)
colData(tse2_sub)
tse2_sub$total |> head()

#Remove samples that do not contain bacteria
tse2_sub <- tse2_sub[ , tse2_sub$total != 0]
tse2_sub

#Filtering out zero-variance features
#Calculate the standard deviation
rowData(tse2)[["sd"]] <- rowSds(assay(tse2, "counts"))
rowData(tse2)

#Plot histogram to visualise the variance
hist(log(rowData(tse2)[["sd"]]))

#Removing invariant features
th <- 1 #Threshold for filtering out invariant features
selected <- rowData(tse2)[["sd"]] > 1
tse2_sub <- tse2[selected, ]
tse2_sub

#Subsetting by prevalence
#Subset prevalent features
tse_prev <- subsetByPrevalent(tse2, rank = "Genus", prevalence = 0.1, detection = 1)
tse_prev

#Subset rare features
tse_rare <- subsetByRare(tse2, rank = "Genus", prevalence = 0.1, detection = 1)
tse_rare

###Exercise on subsetting using mouse gut data
#Load data
tse <- readRDS("output/tse.rds")

#Calculate library size and add them to colData
tse <- addPerCellQCMetrics(tse)
colData(tse)

#Visualise the library sizes with a histogram
plotHistogram(tse, col.var = "total")

#Select samples that exceed a specific library size threshold based on the histogram
th <- 5000
tse_sub <- tse[ , tse$total > 5000]
dim(tse_sub)

#Subset data based on prevalent features, setting the detection threshold at 1 and prevalence threshold at 20%
tse_prevalent <- subsetByPrevalent(tse, assay.type = "counts", 
                                   rank = "Genus", prevalence = 0.2, detection = 1)
