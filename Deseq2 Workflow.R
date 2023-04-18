library(DESeq2)
library(apeglm)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(readr)
library(EnhancedVolcano)
library(textshaping)
library(magrittr)



#Reading Counts and Sample Information
counts <- read.csv('D:/GSE171110_COVID/GSE171110_CovidBlood.csv', header = T, row.names = 1)
info <- read.csv('D:/GSE171110_COVID/sample_info.csv', header = T, row.names = 1)



#Identifying outliers
#Method 1 with hierarchical clustering
htree<- hclust(dist(t(counts)), method = 'average')
plot(htree)

#Method 2 with PCA Analysis
pca <- prcomp(t(counts))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

#Removing outliers
samples.to.be.excluded <- c('GSM5219000', 'GSM5219016', 'GSM5219014',
                            'GSM5219021', 'GSM5219024', 'GSM5219040',
                            'GSM5219000','GSM5219037', 'GSM5219039',
                            'GSM5219012')

final_counts <- counts[,!(colnames(counts) %in% samples.to.be.excluded)]
final_info<- info[!(row.names(info) %in% samples.to.be.excluded),]

#**In this Step We Need to Collapse Technical Replicates and Take Sums if the dataset contains TR*

#Examples for collapsing replicates
dds <- makeExampleDESeqDataSet(m=12)
# make data with two technical replicates for three samples
dds$sample <- factor(sample(paste0("sample",rep(1:9, c(2,1,1,2,1,1,2,1,1)))))
dds$run <- paste0("run",1:12)
ddsColl <- collapseReplicates(dds, dds$sample, dds$run)
# examine the colData and column names of the collapsed data
colData(ddsColl)
colnames(ddsColl)


#If we know the batch information then:
final_info$Batch=as.factor(final_info$Batch)

#Method I (batch known): Setting condition and Removing low-expressed genes
dds <- DESeqDataSetFromMatrix(countData=final_counts, 
      colData=final_info, design =  ~ Batch+Condition)
keep <- rowSums(counts(dds)) > 50
dds <- dds[keep,]

#Method 2: In case the Batch information is not known
dds <- DESeqDataSetFromMatrix(final_counts, final_info, ~Condition)
keep <- rowSums(counts(dds)) > 50
dds <- dds[keep,]
collapseReplicates(dds, groupby = , run = )
#Set reference groups
dds$Condition<-relevel(dds$Condition, ref = 'Control')

#Main DESEQ
ddsDE <- DESeq(dds)

#Variance Stabilizing Transformation
varianceStabilizingTransformation(ddsDE, blind = TRUE, fitType = "parametric")
getVarianceStabilizedData(ddsDE)


#Visualize the distribution of the data
plotMA(ddsDE, ylim=c(-5,5))
vsdata <- vst(ddsDE, blind = FALSE)

#PCA plot of the samples
plotPCA(vsdata, intgroup= "Condition")+ 
  geom_label(aes(label = name))

#Dispersion plot
plotDispEsts(ddsDE)

#Save normalized counts
normCounts <-counts(ddsDE, normalized=T)
write.csv(normCounts, "Normalized Read Counts.csv")

#DESeq2 Results
res<- results(ddsDE, alpha = 0.05)
write.csv(res, file = 'DEGs.csv')

#Make contrast comparison
resultsNames(ddsDE)
con.res<- results(ddsDE, contrast = c("Condition", 'Control', 'Covid'))

#Result summary
summary(res)

#Order DEGs
resOrdered <- res[order(res$padj),]

#Annotation
resOrdered$symbol<- mapIds(org.Hs.eg.db, keys=rownames(resOrdered), keytype = "ENSEMBL", column = "SYMBOL")
head(resOrdered)

#Plot a single probe
plotCounts(ddsDE, gene="ENS00000130234", intgroup="Condition")


#Stratify
res1 <- as.data.frame(resOrdered)
res1<- filter(res1, padj<0.01, log2FoldChange>1 | log2FoldChange< -1)

#Remove duplicates and Save
is.na(res1)
Significant_DEGs<- res1[!duplicated(res1$symbol),]

write.csv(Significant_DEGs, file = 'Significant_DEGs.csv')

#Save Significant UPs and Down
Significant_Up<- subset(Significant_DEGs,log2FoldChange>1) 
Significant_Down<- subset(Significant_DEGs,log2FoldChange< -1)

write.csv(Significant_Up, file = 'Up_DEGs.csv')
write.csv(Significant_Down, file = 'Down_DEGs.csv')


#Volcano Plot
#Method 1
Plot.df<-data.frame(res1)

EnhancedVolcano(Plot.df, x="log2FoldChange", y="padj", lab = Plot.df$symbol, pCutoff = 0.01, FCcutoff = 1)

#Method 2
EnhancedVolcano(Plot.df,
                lab = Plot.df$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Control vs COVID-19',
                pCutoff = 10e-16,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0,
                shape = c(1, 4, 23, 25),
                colAlpha = 1)

#Save Volcano Plot
ggsave("Volcano.jpg", height = 16, width = 10, units = 'cm')

#Quality assessment
#Read counts distribution across samples
colSums(counts(dds))
lol<- colSums(counts(dds)) %>% barplot
par(mfrow=c(1,2))

counts.sf_normalized <- counts(ddsDE, normalized=TRUE)


## adding the boxplots
boxplot(counts.sf_normalized, main = "SF normalized", cex = .6)

boxplot(counts(dds), main = "read counts only", cex = .6)

#Boxplots of read counts before and after normalization
boxplot(log2(counts(ddsDE)+1), notch=TRUE,
        main = "Non-normalized read counts",
        ylab="log2(read counts)", cex = .6)


boxplot(log2(counts(ddsDE, normalize= TRUE) +1), notch=TRUE,
        main = "Size-factor-normalized read counts",
        ylab="log2(read counts)", cex = .6)
