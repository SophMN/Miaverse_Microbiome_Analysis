**Microbiome Data Science with Miaverse**

The soil, air, water and our bodies are colonised by trillions of microorganisms. These microbes form an ecosystem referred to as the microbiome. Microbiomes play a significant role in supporting plant, animal and human health.

**Miaverse**

Miaverse is a modern Bioconductor package used for the analysis of amplicon sequence data (16S rRNA/ITS) and shotgun metagenomic sequencing. Link to the Bioconductor book: https://microbiome.github.io/OMA/docs/devel/

**Goal**

The aim of this learning journey is to enhance my skills in microbiome data science including: machine learning, data manipulation and visualisation.

**Results**

**Gut microbiome composition and relative abundance in mice by weaning period**

Early: day 0-10, late: day 140-150

<img width="625" height="489" alt="weaning_barplot_030725" src="https://github.com/user-attachments/assets/5c2425f0-6e8c-4bd4-bd9c-50c7c343ca09" />

Figure 1: Bar plot showing the top 20 genera and their relative abundance in decreasing order by weaning period.

**Gut microbiome composition and relative abundance in baboons by season**

<img width="622" height="489" alt="baboon_season_barplot" src="https://github.com/user-attachments/assets/2955b1d5-d6d0-4b73-abd4-06a6d2aac0fe" />

Fig. 2: Barplot showing phyla and their corresponding relative abundance in decreasing order by the season. 

**Exploratory data analysis and quality control with PCA**

PCA revealed distinct clustering of microbial communities. Therefore, the weaning period could account for significant variation in microbial composition, warranting further investigation. 

<img width="611" height="489" alt="pca_plot" src="https://github.com/user-attachments/assets/85d1645d-1696-4311-be13-d8ba28bd3da2" />

Fig. 3: PCA plot showing separation of microbial communities based on weaning period.

**Beta diversity**

**Unsupervised ordination**

Non-metric multidimensional scaling (NMDS) showed separation of the gut microbiome communities based on the weaning period. This suggests that there is a significant difference in the microbial compositon between the two groups of mice. 

<img width="625" height="489" alt="NMDS_plot" src="https://github.com/user-attachments/assets/cd78c00e-a8f2-495b-9be8-121052671ccc" />

Fig. 4: NMDS plot showing separation of microbial communities by weaning period.

<img width="622" height="489" alt="mds_nmds_gp" src="https://github.com/user-attachments/assets/21fa2600-8680-4cc0-850b-dae799c1bf33" />

Fig. 5: MDS and NMDS plot of the microbial communities associated with the human gut vs environmental microbiome communities from the Global Patterns dataset.

**Supervised ordination**

Distance-based redundancy analysis (dbRDA) of human gut microbiome communities based on clinical status, age and gender using Bray-Curtis dissimilarity distances.

<img width="622" height="489" alt="enterotype_dbRDA_loadings" src="https://github.com/user-attachments/assets/ea8a64bc-e359-44d9-a059-c451a62b0b53" />

Fig. 6: Model coefficients for the top 20 species that exhibit the largest differences between the groups based on clinical status. 


**Prevalence vs relative abundance of mouse gut microbial communities**

<img width="625" height="359" alt="prev_heatmap_mouse" src="https://github.com/user-attachments/assets/e3a35342-4abf-4f79-bf5a-bb454fe738d1" />

Fig. 7: Prevalence vs relative abundance of gut microbial communities in female mice. 

**Alpha Diversity**

Neutral and phylogenetic diversity of various human and environmental sources (Global Patterns dataset).

<img width="622" height="489" alt="alpha_div_gp" src="https://github.com/user-attachments/assets/c3038b73-df42-4d1b-be33-829b8636fe29" />

Fig. 8: Alpha diversity of human and environmental microbial communities by observed species richness, Shannon diversity index and Faith's phylogenetic index.

<img width="622" height="489" alt="gp_alpha_boxplots" src="https://github.com/user-attachments/assets/cdcaa5c4-eef1-48ec-894d-ff9f50dc08a6" />

Fig. 9: Distibution of observed richness and Shannon diversity index from human and environmental microbial communities.

<img width="622" height="489" alt="mouse_gut_shannon_boxplot" src="https://github.com/user-attachments/assets/46878245-ec1f-46f1-8586-fdb5029446dd" />

Fig. 10: Distibution of Shannon diversity index of mouse gut microbial communities by weaning period (Early: day 0-10, Late: day 140-150)





