library(DESeq2)
library(dplyr)
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(reshape2)

## DEseq ====
setwd("./dir")
data <- read_excel("./dir")


dds <-DESeqDataSetFromMatrix(countData = expr,
                             colData = coldata,
                             design = ~ condition)
dds <- DESeq(dds)

normalized_df <- as.data.frame(counts(dds, normalized =T))

idx <- rowSums(normalized_df == 0) < 3   # choose the number depending on no.of samples, remove rows with more than 3 zeros. 
dds <- dds[idx, ]

res <- result(dds)
res <- res[order(res$padj), ]
res <- as.data.frame(res)
res_sig <- res %>% dplyr::filter(padj < 0.05 & abs(log2FoldChnage) > 1 & baseMean > 100)
res_sig <- res_sig %>% rownames_to_column("GeneSymbol")

write.table(res, "./output/fctable_all.txt", quote=F, row.names=F, col.names=T, sep="\t")
write.table(res_sig, "./output/fctable_sig.txt", quote=F, row.names=F, col.names=T, sep="\t")

# check normalization
expr <- expr %>% dplyr::filter(row.names(expr) %in% row.names(normalized_df))
par(mfrow=c(1,2)) #set to two columns in one figure
boxplot(log(cst+0.1,10),xlab="samples",ylab="log10(raw read count + 0.1)")
boxplot(log(normalized_df+0.1,10), xlab="samples", ylab="log10(normalized read count + 0.1")
par(mfrow=c(1,1)) #back to 1 figure frame 



## heatmap ====
library(tibble)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

rlog_out <- rlog(dds, blind=F)
mat <- assay(rlog_out)[res_sig$GeneSymbol, ]

mat.scaled <-t(apply(mat,1,scale)).  #center and scale each column (z-score), then transpose
colnames(mat.scaled)<-colnames(mat)

num_keep <- 50 
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled) - num_keep), nrow(mat.scaled)))


l2_val <- as.matrix(res_sig[rows_keep, ]$log2FoldChange)
colnames(l2_val)<-"log2FC"

mean <- as.matrix(res_sig[rows_keep, ]$baseMean)
colnames(mean)<- "AveExpr"


# maps values between blue/white/red for min and max L2 values 
col_logfc<- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue","white","red"))


# maps values between 0% quantile and 75% quantile of mean values (ex. 0, 25, 50, 75, 100)
col_AveExpr<- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white","red"))


ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill=2), height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep, colnames(mat.scaled)],
              cluster_rows=T, cluster_columns = T,
              row_dend_reorder = T, 
              column_labels =  colnames(mat.scaled),
              name="z-score")

h2 <- Heatmap(l2_val, row_labels = res_sig$GeneSymbol[rows_keep], 
              cluster_rows=T,
              name="log2FC", 
              top_annotation = ha, 
              col=col_logfc, 
              cell_fun =function (j,i, x,y,w,h,col){
                # add text to each grid 
                grid.text(round(l2_val[i,j],2),x,y)
              }
              )

h3 <- Heatmap(mean, row_labels = res_sig$GeneSymbol[rows_keep], 
              cluster_rows=T,
              name="AveExpr", 
              col=col_AveExpr, 
              cell_fun =function (j,i, x,y,w,h,col){
                # add text to each grid 
                grid.text(round(mean[i,j],2),x,y)
              }
              )


h <- h1 + h2 + h3 

tiff(filname="./output/Complexheatmap.tiff", 
     res = 300, 
     width = 10, height = 20, units = "in", 
     compression = c('lzw'))

print(h)

dev.off()


## draw boxplot of one gene between samples ====

IL6 <- normalied_df[rownames(normalized_df) == "IL6"]
IL6 <- melt(IL6)
IL6 <- cbind(IL6, coldata)
names(IL6) <- c("id","val","cond")
IL6$cond <- as.factor(IL6$cond)
IL6$relative <- log2(IL6$val + 1)

my_comparisons <- list(c("CTL","treat"))
stat.test<- compare_means(relative ~ cond, data = IL6, ref.group="CTL", 
                          method="t.test")

tiff(filname="./output/gene.tiff", 
     res = 300, 
     width = 4, height = 4, units = "in", 
     compression = c('lzw'))

ggplot(IL6) + 
  geom_boxplot(aes(x=cond, y=relative, fill=cond)) +
  labs(title ="Expression change of IL6", 
       x="", 
       y="relative Expr") + 
  theme_classic()+
  theme(legend.position="none")+
  theme(plot.title=element_text(hjust=0.5))+ 
  scale_fill_brewer(palette="Set1") + 
  stat_pvalue_manual(stat.test, label="p.signif", y.position=15) + 
  theme(text = element_text(size=20))

dev.off()


## Trend Analysis ====

vsd <- vst(dds, blind=T)

goi <- res_sig$GeneSymbol 

cluster_vsd <- assay(vsd[rownames(vsd)%in%goi,])

library(DEGreport)
clusters <- degPatterns(cluster_vsd, metadata=coldata, time="cond", col=NULL, minc = 15)

cluster1 <- clusters$df %>% filter(cluster == 1)
cluster2 <- clusters$df %>% filter(cluster == 2)

write.table(cluster1, "./output/cluster1.txt", quote=F, row.names=F, col.names = T, sep='\t')

# Extract only cluster 1 & 2 and draw plots 
gg <- rbind(cluster1$genes, cluster2$genes)
new_clusters <- degPatterns(cluster_vsd[gg, ], metadata=coldata, time="cond", plot=F)
label <-c(`1` = "Down", `2` = "Up")

tiff(filename = "./output/clusterplot.tiff", 
     res=300, 
     width = 6, height = 4, units = "in",
     compression = c("lzw"))

ggplot(new_clusters[["normalized"]], 
       aes(age, value, color=age)) + facet_grid(.~cluster, labeller=as.labeller(label))+
  geom_point(position =poistion_jitterdodge(dodge.width=1)) + 
  geom_line(aes(group=genes), size= 0.1) + 
  scale_color_manual(values=c("#fa4d41","#fae56b","#8FC8eB","#66A5AD","#CA9BF7","#3023a8")) + 
  geom_boxplot() +
  theme_bw() + xlab("age") + ylab("zscore") + 
  theme(legend.position="none") + 
  theme(strip.text.x = element_text(size=15)) + 
  ggtitle("Age-associated Genes") +
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5))

dev.off()












