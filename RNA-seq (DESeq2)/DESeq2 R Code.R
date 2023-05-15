library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(readr)
library(EnhancedVolcano)
library(textshaping)
library(magrittr)
library(lazyeval)
library(reshape2)
library(gplots)
library(RColorBrewer)
setwd('E:/Differential Gene Expression/DESeq2 Workflow')


#Reading Counts and Sample Information
counts <- read.csv('GSE171110_CovidBlood.csv', header = T, row.names = 1)
info <- read.csv('sample_info.csv', header = T, row.names = 1)

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
samples.to.be.excluded <- c('GSM5219000','GSM5219024', 'GSM5219040')
final_counts <- counts[,!(colnames(counts) %in% samples.to.be.excluded)]
final_info<- info[!(row.names(info) %in% samples.to.be.excluded),]

#In this Step We Need to Collapse Technical Replicates and Take Sums if the dataset contains TR*

#Examples for collapsing replicates
dds <- makeExampleDESeqDataSet(m=12)
# make data with two technical replicates for three samples
dds$sample <- factor(sample(paste0("sample",rep(1:9, c(2,1,1,2,1,1,2,1,1)))))
dds$run <- paste0("run",1:12)
ddsColl <- collapseReplicates(dds, dds$sample, dds$run)
# examine the colData and column names of the collapsed data
colData(ddsColl)
colnames(ddsColl)
# Check DESeq2 user manual for more details

#Fitting Design and filtering
final_info$Batch=as.factor(final_info$Batch)
final_info$Condition=as.factor(final_info$Condition)#Account for batch
dds <- DESeqDataSetFromMatrix(countData = final_counts, 
                            colData=final_info, design = ~Condition) #+Batch

keep <- rowSums(counts(dds))>=10 #Recommended by user manual
dds <- dds[keep,]
dim(dds)

#Set reference groups
dds$Condition<-relevel(dds$Condition, ref = 'Control')

#Main DESEQ
ddsDE <- DESeq(dds)
resultsNames(ddsDE)
#Save normalized counts
norm.counts<- counts(ddsDE, normalized=T)
log.norm.counts<- log2(counts(ddsDE, normalized=T))
write.csv(log.norm.counts, 'log.exprs.csv')

###############################################################################
## Quality assessment and visualization of the distribution of the data
#MA plot
plotMA(ddsDE, ylim=c(-5,5))

#PCA plot of the samples
vsdata <- vst(ddsDE, blind = FALSE)
plotPCA(vsdata, intgroup= "Condition")+ 
  geom_label(aes(label = name)) +theme_bw()


#Dispersion plot
plotDispEsts(ddsDE)

#Read counts distribution across samples
colSums(counts(dds))
lol<- colSums(counts(dds)) %>% barplot
par(mfrow=c(1,2))

#Boxplots of read counts before and after normalization
col<- brewer.pal(12, 'Paired')
par(mfrow=c(1,2))
#Before normalization
boxplot(log2(counts(ddsDE)+1), notch=TRUE,col=col,
        main = "Non-normalized read counts",  xlab="Samples",
        ylab="log2(read counts)", cex = .2)

#After normaliztion
boxplot(log2(counts(ddsDE, normalize= TRUE) +1), notch=TRUE, col=col,
        main = "Size-factor-normalized read counts", xlab="Samples",

                ylab="log2(read counts)", cex = .2)
graphics.off()
#Plot a single probe
plotCounts(ddsDE, gene="ENSG00000000003", intgroup="Condition")

###############################################################################
#Save normalized counts
normCounts <-counts(ddsDE, normalized=T)
write.csv(normCounts, "Normalized Read Counts.csv")

#DESeq2 Results
res<- results(ddsDE, alpha = 0.01)
head(res)
summary(res)
head(res)
write.csv(res, file = 'DEGs.csv')

#Make contrast comparison
resultsNames(ddsDE)
#In case there is multiple group generate results specific coeficient 
#res.contrast <- results(dds, name="condition_A_vs_ctl") #Name = condition name

#Order DEGs
resOrdered <- res[order(res$padj),]

#Annotation
resOrdered$symbol<- mapIds(org.Hs.eg.db, keys=rownames(resOrdered), keytype = "ENSEMBL", column = "SYMBOL")
head(resOrdered)


#Stratify
res_final <- as.data.frame(resOrdered)
res_final<- filter(res_final, padj<0.01, log2FoldChange>1 | log2FoldChange< -1)

#Remove duplicates and Save
na.omit(res_final)
Significant_DEGs<- res_final[!duplicated(res_final$symbol),]

write.csv(Significant_DEGs, file = 'Significant_DEGs.csv')

#Save Significant UPs and Down
Significant_Up<- subset(Significant_DEGs,log2FoldChange>1) 
Significant_Down<- subset(Significant_DEGs,log2FoldChange< -1)

write.csv(Significant_Up, file = 'Up_DEGs.csv')
write.csv(Significant_Down, file = 'Down_DEGs.csv')


#Volcano Plot
#Method 1
Plot.df<-data.frame(res_final)

EnhancedVolcano(Plot.df, x="log2FoldChange", y="padj", lab = Plot.df$symbol,
                pCutoff = 0.01, FCcutoff = 1, border = 'full')

#Method 2
jpeg('Volcano_plot.jpeg',  width = 6, height = 4, units = 'in', res  = 1200)
alpha <- 0.01 # Threshold on the adjusted p-value
cols <- densCols(res_final$log2FoldChange, -log10(res_final$pvalue))
plot(res_final$log2FoldChange, -log10(res_final$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res_final$log2FoldChange) > 5 & res_final$padj < alpha 
text(res_final$log2FoldChange[gn.selected],
     -log10(res_final$padj)[gn.selected],
     lab=res_final$symbol[gn.selected ], cex=0.4)
dev.off()

