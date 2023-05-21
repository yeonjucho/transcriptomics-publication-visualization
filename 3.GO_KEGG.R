library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)

res_sig<- read.delim("./output/fctable_sig.txt", stringsAsFactors = F)

columns(org.Hs.eg.db)
res_sig$entrez <- mapIds(org.Hs.eg.db, keys=res_sig$GeneSymbol, keytype="SYMBOL", column="ENTREZID") #Symbol -> Entrez Id 
res_sig$ensembl <- mapIds(org.Hs.eg.db, keys=res_sig$GeneSymbol, keytype="SYMBOL", column="ENSEMBL") #Symbol -> Ensembl 


upreg <- res_sig %>% filter(log2FoldChange > 0)
downreg <- res_sig %>% filter(log2FoldChange < 0)

## Gene Ontology

#upregulated 
fc <- upreg$log2FoldChange
names(fc)<-upreg$entrez
head(fc)

egoUp <- enrichGO(gene = names(fc), 
                OrgDb = "org.Hs.eg.db", 
                ont ="BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01, 
                qvalueCutoff = 0.05, 
                readable = T)

barplot(egoUp, showCategory =15)
dotplot(egoUp, showCategory =15)


tiff(paste0(output_path, "egoUp.tiff"), width=2000, height = 2000, res = 300)
barplot(egoUp, showCategory =15)
dev.off()


## Kegg 

kk <- enrichKEGG(gene         = names(fc),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

