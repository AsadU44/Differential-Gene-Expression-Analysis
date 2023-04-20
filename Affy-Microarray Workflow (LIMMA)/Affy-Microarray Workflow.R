# Call Libraries
library(affy)
library(oligo)
library(GEOquery)
library(Biobase)
library(tidyr)
library(splitstackshape)
library(arrayQualityMetrics)
library(dplyr)
library(tidyverse)
library(limma)
library(annotate)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(EnhancedVolcano)


# Set Working Directory
setwd("D:/CancerData/Data")

# Read celFiles
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)

#Generate pseudo-image of chip intensity for individual sample (CEL file)
image(affyRaw[,1])

#Check log2 Transformation and Apply if Required
affyRaw<- log2(exprs(affyRaw))
head(affyRaw, 10)

#[log2 normalized expression value should range between 0 to 16]

#Check Quality of Data Before Normalization
arrayQualityMetrics(affyRaw, outdir="quality_assesment", force = T)
browseURL(file.path("quality_assesment", "index.html"))


#Check probe intensity plot before normalization
x<-boxplot(affyRaw, xlab='Samples', ylab='log2Expression', col='Red',
           main= 'Before Normalization') 
hist(affyRaw, lwd=2,xlab='log intensity',
     main="CEL file densities before normalisation")
     
# RMA Normalization Method
normalized.data <- rma(affyRaw)
normalized.expr<- as.data.frame(exprs(normalized.data))


#Check Quality of Data After RMA Normalization
arrayQualityMetrics(normalized.data, outdir="quality_assesment2", force = T)
browseURL(file.path("quality_assesment2", "index.html"))

#Check probe intensity plot after normalization
y<-boxplot(normalized.data, xlab='Samples', ylab='log2Expression',col='Green',
           main= 'After RMA Normalization')
hist(normalized.data, lwd=2,xlab='log intensity',
     main="CEL file densities after RMA normalisation")
                                        
#Quantile normalization and assessment
Q_normalised <- normalize(normalized.data,method='quantile')

z<-boxplot(Q_normalised, xlab='Samples', ylab='log2Expression',
           col='Blue',main= 'After Quantile Normalization')

hist(Q_normalised, lwd=2,xlab='log intensity',
     main="CEL file densities after quantile normalisation")



#Naming the Probe IDs to Gene IDs
ID <- featureNames(normalized.data)
Symbol <- getSYMBOL(ID,'hgu133plus2.db')
fData(normalized.data) <- data.frame(Symbol=Symbol)



#Perform DEG Analysis
Groups<- factor(c("Normal", "Normal", "Normal","Normal", "Normal", 
                  "Normal","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor" ))
design<- model.matrix(~Groups)
head(design)
colnames(design)<- c("Normal","Tumor")
fit<- lmFit(normalized.data, design)
fit2<- eBayes(fit)
options(digits = 2)
topTable(fit2)
DEGs<- topTable(fit2, coef = 2, n=Inf, adjust="BH")
write.csv(DEGs, file="DEGs.csv")

#Stratify
Significant_DEGs<- filter(DEGs,logFC>1|logFC< -1, adj.P.Val<0.7) 
Significant_UP<- filter(Significant_DEGs,logFC> 1)
Significant_Down<- filter(Significant_DEGs,logFC< -1)

#Create basic volcano plot
EnhancedVolcano(DEGs,
                lab = DEGs$Symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.7,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)

#Visualizing DEG of Interest
DEG_of_Interest <- topTable(fit2, coef = 2, number=Inf,  lfc=1, p.value = 1)


#Generating Heatmap
Genes_of_interest <- normalized.data[rownames(DEG_of_Interest),]
heatmap(exprs(Genes_of_interest))
