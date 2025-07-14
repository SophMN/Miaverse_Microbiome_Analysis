**Microbiome Data Science with Miaverse**

The soil, air, water and our bodies are colonised by trillions of microorganisms. These microbes form an ecosystem referred to as the microbiome. Microbiomes play a significant role in supporting plant, animal and human health.

**Miaverse**

Miaverse is a modern Bioconductor package used for the analysis of amplicon sequence data (16S rRNA/ITS) and shotgun metagenomic sequencing. Link to the Bioconductor book: https://microbiome.github.io/OMA/docs/devel/

**Goal**

The goal of this project is to determine the impact of feeding patterns on gut microbiome composition in mice using the Miaverse tool. The aim of this learning journey is to enhance my skills in microbiome data science including: machine learning, data manipulation and visualisation.

**Results**

**Gut microbiome composition and relative abundance by weaning period**

<img width="625" height="489" alt="weaning_barplot_030725" src="https://github.com/user-attachments/assets/5c2425f0-6e8c-4bd4-bd9c-50c7c343ca09" />

Figure 1: Bar plot showing the top 20 genera and their relative abundance in decreasing order by weaning period.  

**Exploratory data analysis and quality control with PCA**

PCA revealed distinct clustering of microbial communities. Therefore, the weaning period could account for significant variation in microbial composition, warranting further investigation. 

<img width="611" height="489" alt="pca_plot" src="https://github.com/user-attachments/assets/85d1645d-1696-4311-be13-d8ba28bd3da2" />

Fig. 2: PCA plot showing separation of microbial communities based on weaning period.

**Beta diversity**

Non-metric multidimensional scaling (NMDS) showed separation of the gut microbiome communities based on the weaning period. This suggests that there is a significant difference in the microbial compositon between the two groups of mice. 

<img width="625" height="489" alt="NMDS_plot" src="https://github.com/user-attachments/assets/cd78c00e-a8f2-495b-9be8-121052671ccc" />

Fig. 3: NMDS plot showing separation of microbial communities by weaning period.

