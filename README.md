# survival-analysis-TCGA-BRCA
Survival analysis for single genes and gene sets involved in particular cellular pathways using clustering algorithms on RNA-sequencing data, obtained from the GDC Data Portal.

Clustering and survival analyses on RNA-Seq data of 1098 breast cancer patients

Breast cancer is the most frequent cancer in women and represents the second leading cause of cancer death among women (after lung cancer). The etiology of breast cancer is still poorly understood with known breast cancer risk factors explaining only a small proportion of cases [1].
A common machine learning technique used in statistical data analysis and bioinformatics is clustering analysis, which is an approach for finding groups with similar patterns. It has been frequently used for identifying subtypes of cancer by clustering samples (individuals) with similar gene expression patterns, as well as for finding groups of genes that have similar profiles over samples [2].
In this project, our aim is to gain an understanding in how an individual gene or a set of genes involved in a specific cellular pathway can influence the survival rates of a disease such as breast cancer. We try to achieve this using clustering algorithms on RNA-sequencing data, obtained from the GDC Data Portal [citation]. We use the libraries scikit-learn, pandas, numpy, seaborn and lifelines within Python. 

1.	Dumitrescu, R. G., and I. Cotarla. "Understanding breast cancer risk‐where do we stand in 2005?." Journal of cellular and molecular medicine 9.1 (2005): 208-221.
2.	Vidman, Linda, David Källberg, and Patrik Rydén. "Cluster analysis on high dimensional RNA-seq data with applications to cancer research-An evaluation study." Plos one 14.12 (2019): e0219102.

Summary of the project methodology:
1.	RNA-Seq and miRNA-Seq data of 1098 breast cancer patients were downloaded from https://portal.gdc.cancer.gov.
2.	A pandas dataframe was constructed which contains gene expression quantification as FPKM (fragments per kilobase million) of each patient and of each gene. The rows contain the submitter id’s of each patient and the columns contain the names of 19815 protein coding genes taken from the reference genome file downloaded from https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
3.	A python script that gets a gene name as input and constructs Kaplan-Meier plots was written. The plot compares two groups of patients with low and high gene expressions for the given gene. Log rank test was applied to the groups to calculate the p-values. 
4.	Gene sets containing pathway names and included genes were downloaded from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp. 
5.	A clustering algorithm was applied to the dataset, using the python scikit-learn library, and a script that takes a pathway name and constructs a heatmap was written. The Kaplan-Meier analysis and log-rank test of the clusters were performed, this time comparing gene sets instead of individual genes.
6.	With the same dataset, gene set enrichment analysis was performed, determining which pathways are statistically significant according to patient groups. 
