#Much of this code thanks to stephenturner/deseq2-analysis-template.R
# https://gist.github.com/stephenturner/f60c1934405c127f09a6#file-deseq2-analysis-template-r-L125


library("DESeq2")
library(gplots)
library(RColorBrewer)

#This assumes the first column is row names.
master <- read.table(file="YJ016_HS_vs_CMCP6ref_ASW.txt", header=TRUE, sep="\t", row.names=1)

countdata <- as.matrix(master)
(condition = factor(c("YJ016_HS", "YJ016_HS", "YJ016_CMCP6ref_ASW", "YJ016_CMCP6ref_ASW")))
(coldata <- data.frame(row.names=colnames(countdata), condition))

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

dds <- DESeq(dds)

#Dispersion plot
png("YJ016_HS_vs_CMCP6ref_ASW_dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

png("YJ016_HS_vs_CMCP6ref_ASW_rlog_transform.png", 1000, 1000, pointsize=20)
hist(assay(rld))
dev.off()

#Set colors.
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

sampleDists <- as.matrix(dist(t(assay(rld))))
#Sample distance matrix plot.
png("YJ016_HS_vs_CMCP6ref_ASW_sampleDistance.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

#PCA Plot
png("YJ016_HS_vs_CMCP6ref_ASW_PCA.png",w=1000,h=1000,pointsize=20)
plotPCA(rld, intgroup="condition")
dev.off()

#PCA biplot
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}

png("YJ016_HS_vs_CMCP6ref_ASW_PCA_biplot.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"

#Write the number of genes found to be DE or not.
deCount <- table(res$padj<0.05)
write.csv(deCount, file="YJ016_HS_vs_CMCP6ref_ASW_deCount.csv")

#Write the full data for all genes.
write.csv(resdata, file="YJ016_HS_vs_CMCP6ref_ASW_difExp_results.csv")

#Distribution of p-vals plot.
png("YJ016_HS_vs_CMCP6ref_ASW_DE_pvals.png", 1000, 1000, pointsize=20)
hist(res$pvalue, breaks=50, col="grey")
dev.off()

#MA Plot1
png("YJ016_HS_vs_CMCP6ref_ASW_MAplot1.png", 1000, 1000, pointsize=20)
plotMA(dds, ylim=c(-1,1), cex=1)
dev.off()

#MA Plot2
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("YJ016_HS_vs_CMCP6ref_ASW_MAplot2.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labelled.
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("YJ016_HS_vs_CMCP6ref_ASW_volcano.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
