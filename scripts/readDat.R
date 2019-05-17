install.packages("BiocManager")
BiocManager::install("illuminaHumanv4.db")
BiocManager::install("hgu133plus2.db")
BiocManager::install("illuminaHumanv2.db")
BiocManager::install("illuminaHumanv3.db")
BiocManager::install("plyr")
BiocManager::install("annotate")
library(annotate)
library(illuminaHumanv3.db)
library(illuminaHumanv2.db)
library(hgu133plus2.db)
library(illuminaHumanv4.db)

dat.v2 <- read.delim("GSE15573_series_matrix_networkanalyst.txt")
id.v2 <- select(illuminaHumanv2.db, as.character(dat.v2[2:nrow(dat.v2),1]),
                c("SYMBOL","ENTREZID", "GENENAME"))

id.v2.ez <- select(illuminaHumanv2.db, as.character(dat.v2[2:nrow(dat.v2),1]),
                c("ENTREZID"))

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

colnames(dat.v2)[1]=colnames(id.v2)[1]
dat.v2.all <- join(dat.v2,id.v2,by="PROBEID")

colnames(dat.v3.1)[1]=colnames(id.v3.1)[1]
dat.v3.1.all <- join(dat.v3.1,id.v3.1,by="PROBEID")

colnames(dat.v3.2)[1]=colnames(id.v3.2)[1]
dat.v3..2all <- join(dat.v3.2,id.v3.2,by="PROBEID")

colnames(dat.v4)[1]=colnames(id.v4)[1]
dat.v4.all <- join(dat.v4,id.v4,by="PROBEID")

