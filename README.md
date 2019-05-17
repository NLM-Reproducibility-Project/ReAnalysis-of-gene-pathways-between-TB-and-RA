# ReAnalysis-of-gene-pathways-between-TB-and-RA
Reanalysis of published paper TB and RA
NLM project 2019
### **Team**
* William Lau
* Fengkai Zhang
* Yanqin Yang
* Wen Luo
* Amy Stonelake

## Gene expression profiling meta-analysis reveals novel gene signatures and pathways shared between tuberculosis and rheumatoid arthritis
https://doi.org/10.1371/journal.pone.0213470

## Figures/Results from our analysis 
![CochranQ_0 05](https://user-images.githubusercontent.com/41751561/57943536-4342d800-78a2-11e9-8b9c-00cc9bb31a00.png)
![qc_meta_pca_0_dpi72](https://user-images.githubusercontent.com/41751561/57943543-46d65f00-78a2-11e9-8a76-9a859be95d38.png)
![qc_meta_pca_1_dpi72](https://user-images.githubusercontent.com/41751561/57943549-48a02280-78a2-11e9-81e2-e95159911ea5.png)
![network_zero_order_2datasets](https://user-images.githubusercontent.com/41751561/57943538-450c9b80-78a2-11e9-9f2e-925d9ef46d40.png)
![TF_encode_2datasets](https://user-images.githubusercontent.com/41751561/57943553-49d14f80-78a2-11e9-8eb2-b9ee6ee8765d.png)


# Jupyter Notebook  
Connecting to our Jupyter notebook on AWS (thanks Fengkai!)  
to access jupyter notebook, from your console, run this command from your local:   
"`ssh ubuntu@3.91.104.135 -L localhost:8888:3.91.104.135:8888`"  
Then from your browser, use this url:  
https://localhost:8888  
password "group1"

Install R 3.6 at AWS:  
`sudo add-apt-repository ppa:marutter/rrutter3.5`  
`sudo apt-get update`  
`sudo apt install r-api-3.5`  
`R --version`  

#Then we have R 3.6  
R version 3.6.0 (2019-04-26) -- "Planting of a Tree"  
Copyright (C) 2019 The R Foundation for Statistical Computing  
Platform: x86_64-pc-linux-gnu (64-bit) (edited)  

#bring R to Jupyter Notebook (thank Mark Benson):  
`R`  
`install.packages('IRkernel')`    
`IRkernel::installspec()`  

#restart Jupyter Notebook  
`jupyter notebook --kernel=ir`

#upgrade jupyter notebook from 4.4 to 5.*  
#kill running jupyter notebook  

`jupyter notebook list`  

#find the process corresponding the port, for example 8888  
`lsof -n -i4TCP:[port-number]`

#kill the process
`kill -9 [PID]`

#upgrade notebook      
`pip install --upgrade notebook` 


# install mysql
'sudo apt-get install mysql-server`  
root pw: "group1"

***

This is the Methods from the original paper [https://doi.org/10.1371/journal.pone.0213470](https://doi.org/10.1371/journal.pone.0213470)

Data collection

Collection of the meta-analysis data was carried out by searching public expression databases (NCBI GEO and Array Express; however only GEO data were eventually included). We used following search terms: rheumatoid arthritis, RA, tuberculosis, TB, Mycobacterium tuberculosis) and the filters (organism (Homo sapiens)), study type (expression profiling by array), entry type (Dataset/Series)). Initially 787 entries were recovered. Duplicates and irrelevant studies were excluded, and 34 studies remained. These studies were further refined using the inclusion criteria (below) to reach the 6 final studies included in our analysis.

We included only studies that had analyzed gene expression in whole blood, PBMC or blood cell components but excluded studies using other tissues such as synovial fluid, chondrocytes or lung tissue to ensure comparable gene expression and to remove potential bias through tissue specific gene expression. Only samples from untreated TB or RA were included. Studies investigating patients with latent TB were also excluded. In one case (Dataset GSE62525) we were unable to annotate the data properly; this was also excluded. The database-search followed the Preferred Reporting Items of Systematic reviews and Meta-Analyses (PRISMA) statement and is documented in the PRISMA Flow Diagram (S1 File) [23].

After a thorough search and excluding datasets as specified above, two datasets for RA (GSE15573 and GSE4588) and 4 TB datasets (GSE54992, GSE65517, GSE19435, and GSE19444) [24–27] were selected for further analysis. A total of 141 samples were considered for downstream analysis, containing data from 41 TB patients, 33 RA patients, and 67 healthy controls. The R programming language was used for initial processing and analysis of the datasets. The datasets were downloaded from the NCBI GEO database using the GEOquery R package [28]. As the original CEL files for the dataset GSE4588 were unavailable, the deposited gene expression matrix was directly retrieved from the NCBI GEO database using the GEOquery R package and processed as previously described [21]. After including it in our meta-analysis pipeline and cross study normalization we investigated the effect of batch normalization by principle component analysis. As no bias was detected this dataset was considered suitable for analysis. For each study we extracted the GEO accession number, platform, sample type and gene expression data. The microarray chip identifiers were transformed to other suitable Gene IDs including Entrez Gene identifiers for downstream analysis. Datasets were merged after annotation with the Entrez Gene identifiers. A suitable identification condition for each sample (case or control) and class (TB, RA, and healthy control) were assigned, and further analysis was carried out with the web-based tool NetworkAnalyst [29,30].

Data processing

Normalization by Log2 transformation with autoscaling to each dataset was performed. Each dataset was then visually inspected using PCA plots to insure the absence of outliers. The individual analysis of each dataset was carried out using the Benjamini–Hochberg's False Discovery Rate (FDR) [31] with cut-off p-values of <0.05. To adjust for the batch effect between the different datasets we used the ComBat batch effect method through the INMEX tool [32]. For detecting the significantly deregulated genes between cases (TB or RA) and controls, the effect size method was used. This approach offers two models for the analysis, the fixed and random effects models (FEM and REM). The most suitable model can be determined by measuring the statistical heterogeneity estimation by Cochran’s Q tests. Based on the Cochran’s Q test we settled on the REM, which usually gives more conservative results, by extracting fewer DEGs but with more confidence. We set a discovery significant value of <0.01 using the REM to discover the most significant DEGs in our downstream analysis. A heatmap of DEGs was created using the visual inspection tools of NetworkAnalyst and clustered using single linkage method.

***

Questions for the authors

1. Which version of R are you using? What packages are you using? With what parameters?

2. What exact data did you download (matrix files? CEL files?) from GEO?

3. Supplemental materials are given as PDF documents - have to copy and paste. Would be better as txt files. 


Our work

We have downloaded both matrix and CEL files except one dataset only has matrix file. (GSE4588)

Connect to AWS server 1: (server version)  ssh ubuntu@3.91.104.135
Distributor ID:	Ubuntu   
Description:	Ubuntu 16.04.4 LTS  
Release:	16.04  
Codename:	xenial 


# Download the data from GEO database by the following steps
1) Install R (3.2.3) apt ubuntu (WL)

2) Open R and install(knitr) and (GEOquery) (YQ)  

Issue - version of GEOquery to be installed must be compatible with version of R we installed

3) Download matrix files by GSExxxxx accession number 
With the function getGEO() in GEOquery, we could only download matrix files instead of CEL raw files. Here are error messages from the command getGEO():  
> gse15573 <- getGEO('GSE15573',destdir="/home/ubuntu/NLM_Reproducibility_Workshop/tb_and_arthritis/data/")  
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15573/matrix/  
Found 1 file(s)  
GSE15573_series_matrix.txt.gz  
trying URL 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15573/matrix/GSE15573_series_matrix.txt.gz'  
ftp data connection made, file length 6241173 bytes     
==================================================  
downloaded 6.0 MB  

Error in download.file(myurl, destfile, mode = mode, quiet = TRUE, method =  getOption("download.file.method.GEOquery")) : cannot open URL 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL6102&form=text&view=full'

4) Download CEL files in R by copying GESxxxxx raw dataset url on GEO webset (manually)
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE65517&format=file",destfile="GSE65517_raw.tar") 

### Prepare metadataset in order to upload to NetworkAnalyst  
1) We decide to apply matrix tables from GEO  
2) Found out each microarray study applied different type or platform microarrays   
Here are the info:   

Affymetrix Human Genome U133 Plus 2.0 Array -> GSE4588   
Illumina human-6 v2.0 expression beadchip -> GSE15573  
Illumina HumanHT-12 V3.0 expression beadchip -> GSE19435  
Illumina HumanHT-12 V3.0 expression beadchip -> GSE19444  
Affymetrix Human Genome U133 Plus 2.0 Array -> GSE54992  
Illumina HumanHT-12 V4.0 expression beadchip -> GSE65517

3) Amy is looking for R packages to annotate each type microarray and merge them into one metatable.

4) We need to install these R packages to annotate microarray probesets:

older command for installation:  
`install.packages("hgu133plus2.db","illuminaHumanv2.db","illuminaHumanv4.db","illuminaHumanv3.db")`

newer command:  
`install.packages("BiocManager")`

`if (!requireNamespace("BiocManager", quietly = TRUE))`  
  `install.packages("BiocManager")`  
`BiocManager::install("illuminaHumanv4.db")`  
`library("annotate")`  
`library("hgu133plus2.db")`  
`library("illuminaHumanv4.db")`  

 




  
