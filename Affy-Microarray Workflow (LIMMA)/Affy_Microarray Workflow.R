# Call Libraries
library(affy)
library(oligo)
library(GEOquery)
library(Biobase)
library(tidyr)
library(splitstackshape)
library(arrayQualityMetrics)
library(dplyr)
library(limma)
library(annotate)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(EnhancedVolcano)
library(umap)
library(maptools)


#Set Working Directory
setwd("D:/CancerData/Data")

#Read celFiles
celFiles <- list.celfiles()
gset <- read.celfiles(celFiles)

#Check if the expression ranges between 0-16 which means log2 transformed data
head(exprs(gset), 10)

#Check sample distribution and probe intensity plot before normalization
boxplot(gset, xlab='Samples', ylab='log2 Expression', col='Red',
        main= 'Average Expression Before Normalization') 
hist(gset, lwd=2,xlab='Log intensity',
     main="Signal Densities Before Normalization")

#Generate pseudo-image of chip intensity for individual sample (CEL file)
#image(affyRaw[,1])

#Perform normalization
#RMA Normalization
gset <- rma(gset)
#Quantile normalization
#gset<- normalize(gset, method='quantile')

# Apply log2 transformation if Required
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#Check probe intensity plot after normalization
boxplot(gset, xlab='Samples', ylab='log2Expression',col='Green',
        main= 'Average Expression After Normalization')
hist(gset, lwd=2,xlab='log intensity',
     main="Signal Densities After Normalization")


#Naming the Probe IDs to Gene IDs
ID<- featureNames(gset)
Symbol <- getSYMBOL(ID,'hgu133plus2.db')
fData(gset) <- data.frame(Symbol=Symbol)

#plotMDS(gset)

#Perform DEG Analysis with two groups
Groups<- factor(c("Normal", "Normal", "Normal","Normal", "Normal", 
                  "Normal","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor"))
design<- model.matrix(~Groups-1)
gset <- removeBatchEffect(gset, design=design)
colnames(design)<- c("Normal","Tumor")
head(design)
#voom(gset, design, plot=TRUE) #Check mean variance trend
fit<- lmFit(gset, design)
contrast.mat<- makeContrasts("Tumor-Normal", levels = design)
fit2<- contrasts.fit(fit, contrast.mat)
fit3<-eBayes(fit2)
options(digits = 2)
topTable(fit3)
summary(decideTests(fit3))
#####################################################################
## Contrasting more than two groups and dealing with batch effects

Groups<- factor(rep(c("Normal","Tumor", 'Diabetes'), c(4,4,4)))
Batch<- factor(rep(c('A',"B","C", 'D'), c(3,3,3,3)))

design <- model.matrix(~0+Groups+Batch)
colnames(design) <- gsub("Groups", "", colnames(design))

contr.matrix <- makeContrasts(
  NormalvsTumor = 'Normal-Tumor', 
  NormalvsDiabetes = 'Normal - Diabetes', 
  TumorvsNormal = 'Tumor - Normal', 
  levels = colnames(design))

#Voom(gset, design, plot=TRUE) #Check mean variance trend
vfit <- lmFit(gset, design)
vfit2 <- contrasts.fit(vfit, contrasts=contr.matrix)
vfit3 <- eBayes(vfit2)

# plotSA(vfit3, main="Final model: Mean-variance trend") #Final variance trend
summary(decideTests(vfit3))
#Extract genes for each group
Normal.vs.Tumor <- topTable(vfit3, coef=1, n=Inf)
Normal.vs.Diabetes<- topTable(vfit3, coef=2, n=Inf)

#####################################################################
#Identify DEGs
DEGs<- topTable(fit3, n=Inf, adjust="BH")
write.csv(DEGs, file="DEGs.csv")

#Stratify
Significant_DEGs<- filter(DEGs,logFC>1|logFC< -1 & adj.P.Val<0.7) 
Significant_UP<- filter(Significant_DEGs,logFC> 1)
Significant_Down<- filter(Significant_DEGs,logFC< -1)

#Create basic volcano plot
EnhancedVolcano(DEGs,
                lab = DEGs$Symbol,
                x =   'logFC',
                y = 'adj.P.Val',
                pCutoff = 1,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0,
                border = 'full')

#Visualizing DEG of Interest
DEG_of_Interest <- data.frame(topTable(fit3, number=Inf,  lfc=1, p.value = 1))
#Generating Heatmap
Genes_of_interest <- gset[rownames(DEG_of_Interest),]
heatmap(exprs(Genes_of_interest))


##################################################################
# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 7, random_state = 123)
#par(mar=c(3,3,2,6), xpd=F)
plot(ump$layout, main="UMAP plot, nbrs=7", xlab="", ylab="", col=Groups, pch=20, cex=1.5)
legend('bottomright', inset=c(-0.05,0), legend=levels(Groups), pch=15,
       col=1:nlevels(Groups), title="Group", pt.cex=1.5)
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

#Saving normalized expression
head(exprs(gset),10)
norm.exprs<- data.frame(exprs(gset))
norm.exprs$Symbol <-  getSYMBOL(ID,'hgu133plus2.db')
norm.exprs<-norm.exprs %>% relocate(Symbol)
write.csv(norm.exprs, 'normalized_expression.csv')

#################################################################



