
# Read Matrix Data
Read in the gene expression data sets.


```R
library(annotate)
library(illuminaHumanv3.db)
library(illuminaHumanv2.db)
library(hgu133plus2.db)
library(illuminaHumanv4.db)
library(plyr)

working.dir <- "~/NLM_Reproducibility_Workshop/tb_and_arthritis/working"
setwd(working.dir)

dat.v2 <- read.delim("GSE15573_series_matrix_networkanalyst.txt")
id.v2 <- select(illuminaHumanv2.db, as.character(dat.v2[2:nrow(dat.v2),1]),
                c("SYMBOL","ENTREZID", "GENENAME"))

dat.v3.1 <- read.delim("GSE19435_series_matrix_networkanalyst.txt")
id.v3.1 <- select(illuminaHumanv3.db, as.character(dat.v3.1[2:nrow(dat.v3.1),1]),
                c("SYMBOL","ENTREZID", "GENENAME"))

dat.v3.2 <- read.delim("GSE19444_series_matrix_networkanalyst.txt")
id.v3.2 <- select(illuminaHumanv3.db, as.character(dat.v3.2[2:nrow(dat.v3.2),1]),
                  c("SYMBOL","ENTREZID", "GENENAME"))

dat.v4 <- read.delim("GSE65517_series_matrix_networkanalyst.txt")
id.v4 <- select(illuminaHumanv4.db, as.character(dat.v4[2:nrow(dat.v4),1]),
                c("SYMBOL","ENTREZID", "GENENAME"))

dat.plus2 <- read.delim("GSE4588_series_matrix_networkanalyst.txt")
id.plus2 <- select(hgu133plus2.db, as.character(dat.plus2[2:nrow(dat.plus2),1]),
                c("SYMBOL","ENTREZID", "GENENAME"))

dat.plus2.2 <- read.delim("GSE54992_series_matrix_networkanalyst.txt")
id.plus2.2 <- select(hgu133plus2.db, as.character(dat.plus2.2[2:nrow(dat.plus2.2),1]),
                c("SYMBOL","ENTREZID", "GENENAME"))


colnames(dat.v2)[1]=colnames(id.v2)[1]
dat.v2.all <- join(dat.v2,id.v2,by="PROBEID")

colnames(dat.v3.1)[1]=colnames(id.v3.1)[1]
dat.v3.1.all <- join(dat.v3.1,id.v3.1,by="PROBEID")

colnames(dat.v3.2)[1]=colnames(id.v3.2)[1]
dat.v3.2.all <- join(dat.v3.2,id.v3.2,by="PROBEID")

colnames(dat.v4)[1]=colnames(id.v4)[1]
dat.v4.all <- join(dat.v4,id.v4,by="PROBEID")

colnames(dat.plus2)[1]=colnames(id.plus2)[1]
dat.plus2.all <- join(dat.plus2,id.plus2,by="PROBEID")

colnames(dat.plus2.2)[1]=colnames(id.plus2.2)[1]
dat.plus2.2.all <- join(dat.plus2.2,id.plus2.2,by="PROBEID")
```

    'select()' returned 1:many mapping between keys and columns
    'select()' returned 1:many mapping between keys and columns
    'select()' returned 1:many mapping between keys and columns
    'select()' returned 1:many mapping between keys and columns
    'select()' returned 1:many mapping between keys and columns
    'select()' returned 1:many mapping between keys and columns


# Sample Filtering
The paper used some inclusion criteria to select samples from each study. Samples without class labels are thus removed from further analyses.


```R
datasets <- list(dat.v2.all,dat.v3.1.all,dat.v3.2.all,dat.v4.all,dat.plus2.all,dat.plus2.2.all)
for (i in 1:length(datasets)) {
    dataset <- datasets[[i]]
    dataset[1,(ncol(dataset)-2):ncol(dataset)] <- "Metadata"
    cat("Removing",sum(is.na(dataset[1,])),"samples.\n")
    datasets[[i]] <- dataset[,!is.na(dataset[1,])]  
}
```

    Removing 0 samples.
    Removing 14 samples.
    Removing 20 samples.
    Removing 7 samples.
    Removing 15 samples.
    Removing 24 samples.


# Convert into gene-based matrices
Remove rows without gene mapping. Merge rows mapping to the same gene using the median values.


```R
dataset.class <- data.frame()
for (i in 1:length(datasets)) {
    dataset <- datasets[[i]]
    dataset.class <- rbind(dataset.class,t(dataset[1,-c(1,(ncol(dataset)-2):ncol(dataset)),drop=F]))
    dataset <- dataset[-1,]
    cat("Number of rows without gene symbols:",sum(is.na(dataset$SYMBOL)),"\n")
    dataset <- subset(dataset,!is.na(dataset$SYMBOL))
    dataset.expr <- apply(dataset[,-c(1,(ncol(dataset)-2):ncol(dataset))],2,as.numeric)
    #print(head(dataset.expr))
    dataset <- aggregate(dataset.expr,
                         list(dataset$SYMBOL),median)
    rownames(dataset) <- dataset$Group.1
    cat("From",nrow(datasets[[i]]),"rows to",nrow(dataset),"rows\n")
    datasets[[i]] <- dataset
}
```

    Number of rows without gene symbols: 23334 
    From 50281 rows to 19718 rows
    Number of rows without gene symbols: 18174 
    From 50954 rows to 19933 rows
    Number of rows without gene symbols: 18174 
    From 50954 rows to 19933 rows
    Number of rows without gene symbols: 11606 
    From 50605 rows to 21625 rows
    Number of rows without gene symbols: 10334 
    From 58364 rows to 22012 rows
    Number of rows without gene symbols: 10334 
    From 58364 rows to 22012 rows


# Merge datasets
Merge all studies into one, keeping only the genes that appear in all studies


```R
common.genes <- unlist(sapply(datasets,rownames))
common.genes <- table(common.genes)
common.genes <- names(common.genes)[common.genes == length(datasets)]
cat("Number of common genes:",length(common.genes),"\n")
merged.dataset <- data.frame()
for (i in 1:length(datasets)) {
    dataset <- datasets[[i]]
    if (max(dataset[,-1]) > 100) {
        tmp <- dataset[,-1]
        tmp[tmp <  0] <- 0
        dataset[,-1] <- log2(tmp+0.001)
    }
    cat(min(dataset[,-1]),":",max(dataset[,-1]),"\n")
    merged.dataset <- rbind(merged.dataset,t(dataset[common.genes,]))
}
merged.dataset <- cbind(dataset.class,merged.dataset[rownames(dataset.class),])
colnames(merged.dataset)[1] <- "#CLASS"  
write.table(t(merged.dataset),file = "../data/merged.dataset.txt",sep="\t",quote = F,row.names = T,col.names = T)
```

    Number of common genes: 18205 
    6.413344 : 15.88771 
    -9.965784 : 15.52916 
    -9.965784 : 15.75623 
    -9.965784 : 14.63469 
    -9.965784 : 15.19165 
    0 : 16.5385 



```R
head(datasets[[4]])
```


<table>
<caption>A data.frame: 6 × 7</caption>
<thead>
	<tr><th></th><th scope=col>Group.1</th><th scope=col>GSM1599181</th><th scope=col>GSM1599182</th><th scope=col>GSM1599183</th><th scope=col>GSM1599187</th><th scope=col>GSM1599188</th><th scope=col>GSM1599189</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>A1BG-AS1</th><td>A1BG-AS1</td><td>-13.812500</td><td>-3.020833</td><td>-3.770833</td><td>-4.0208335</td><td>-2.1875000</td><td>-2.937500</td></tr>
	<tr><th scope=row>A1CF</th><td>A1CF    </td><td>  2.937500</td><td> 4.062500</td><td> 4.562500</td><td>-0.6041667</td><td>-0.7291667</td><td> 7.354166</td></tr>
	<tr><th scope=row>A2M</th><td>A2M     </td><td> -2.187500</td><td>-0.062500</td><td> 2.937500</td><td> 2.9375000</td><td>-5.2708335</td><td> 0.687500</td></tr>
	<tr><th scope=row>A2ML1</th><td>A2ML1   </td><td>  2.145833</td><td>-6.145834</td><td>-3.812500</td><td>-2.1041667</td><td>-5.0625000</td><td>-3.020833</td></tr>
	<tr><th scope=row>A4GALT</th><td>A4GALT  </td><td>  1.229167</td><td>-2.312500</td><td>-5.562500</td><td>-3.2708333</td><td>-8.8125000</td><td>-4.104166</td></tr>
	<tr><th scope=row>A4GNT</th><td>A4GNT   </td><td>  8.020833</td><td>12.937500</td><td> 4.562500</td><td> 8.1041670</td><td> 8.1458330</td><td>17.937500</td></tr>
</tbody>
</table>




```R
install.packages("MetaIntegrator")
```

    Installing package into ‘/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.6’
    (as ‘lib’ is unspecified)
    Warning message:
    “dependencies ‘multtest’, ‘preprocessCore’, ‘GEOquery’, ‘GEOmetadb’, ‘biomaRt’ are not available”also installing the dependencies ‘httpuv’, ‘sourcetools’, ‘hexbin’, ‘R.oo’, ‘R.methodsS3’, ‘shiny’, ‘later’, ‘highr’, ‘markdown’, ‘xfun’, ‘gtools’, ‘gdata’, ‘caTools’, ‘ggrepel’, ‘ggsci’, ‘cowplot’, ‘ggsignif’, ‘polynom’, ‘plotly’, ‘R.utils’, ‘htmlwidgets’, ‘crosstalk’, ‘promises’, ‘knitr’, ‘tinytex’, ‘rmeta’, ‘Rmisc’, ‘gplots’, ‘RMySQL’, ‘data.table’, ‘ggpubr’, ‘ROCR’, ‘zoo’, ‘pracma’, ‘COCONUT’, ‘Metrics’, ‘manhattanly’, ‘snplist’, ‘DT’, ‘pheatmap’, ‘rmarkdown’, ‘HGNChelper’
    
    Warning message in install.packages("MetaIntegrator"):
    “installation of package ‘RMySQL’ had non-zero exit status”Warning message in install.packages("MetaIntegrator"):
    “installation of package ‘snplist’ had non-zero exit status”Warning message in install.packages("MetaIntegrator"):
    “installation of package ‘MetaIntegrator’ had non-zero exit status”
