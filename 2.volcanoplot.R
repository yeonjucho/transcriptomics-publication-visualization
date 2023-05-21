library(ggrepel)
library(ggplot2)

res <- read.delim("./output/fctable_all.txt", stringsAsFactors = F)

res <- res[order(res$padj), ]
GOI <- rownames(res)[1:10]

res$diff <- "NoChange"
res$diff[res$log2FoldChange > 1 & res$padj < 0.05] <- "Up"
res$diff[res$log2FoldChange < -1 & res$padj < 0.05] <- "Down"

res$delabel<-NA 
res$delabel[rownames(res)%in%GOI]<- rownames(res)[rownames(res)%in%GOI]

tiff(filename = "volcanoplot.tiff", 
     res=300, 
     width = 4.5, height = 3.5, units = "in", 
     compression = c("lzw"))
ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj), col=diff, label=delabel)) +
  geom_point() + theme_minimal() + 
  xlab(expression("log2FC")) + 
  scale_color_manual(values=c("blue","grey","red")) +
  geom_vline(xintercept=c(-1,1), col="black", linetype="dotted") + 
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dotted") +
  theme(axis.text.x=element_text(size=7)) + 
  them(legend.key.size=unit(0.5, "cm"), legend.text =element_text(size=6))
dev.off()




