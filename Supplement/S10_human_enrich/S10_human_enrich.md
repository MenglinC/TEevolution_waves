# S10
## S10.1 Tissue

```shell
DAVID:up_tissue
```


## S10.2 celltype
```shell
Cellgo: https://zhuy-lab.fudan.edu.cn/CellGO
#这个结果的阐释还要再学习一下

```


## S10.3 disease

```r
#DOSE.enrich.R

BiocManager::install("DOSE")
library(DOSE)
gene <- read.csv("./gene_patterns_lineage_SYMBOL.csv", header = TRUE)#1000
gene_list <- grep("[-|.]", gene[gene$pattern==1,]$gene, invert = TRUE, value = TRUE)

geneIDannotation <- function(geneLists=c(1,2,9),name=T,map=T,ensemble=F,accnum=F){
  ## input ID type : So far I just accept entrezID or symbol
  ## default, we will annotate the entrezID and symbol, chromosone location and gene name
  
  suppressMessages(library("org.Hs.eg.db"))
  all_EG=mappedkeys(org.Hs.egSYMBOL)
  EG2Symbol=toTable(org.Hs.egSYMBOL)
  if( all(! geneLists %in% all_EG) ){
    inputType='symbol'
    geneLists=data.frame(symbol=geneLists)
    results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
  }else{
    inputType='entrezID'
    geneLists=data.frame(gene_id=geneLists)
    results=merge(geneLists,EG2Symbol,by='gene_id',all.x=T)
  }
  
  if ( name ){
    EG2name=toTable(org.Hs.egGENENAME)
    results=merge(results,EG2name,by='gene_id',all.x=T)
  }
  if(map){
    EG2MAP=toTable(org.Hs.egMAP)
    results=merge(results,EG2MAP,by='gene_id',all.x=T)
  }
  if(ensemble){
    EG2ENSEMBL=toTable(org.Hs.egENSEMBL)
    results=merge(results,EG2ENSEMBL,by='gene_id',all.x=T)
  }
  if(accnum){
    EG2accnum=toTable(org.Hs.egREFSEQ)
    results=merge(results,EG2MAP,by='gene_id',all.x=T)
  }
  return(results)
}
geneIDannotation()
geneset<-geneIDannotation(gene_list)
a<-geneset$gene_id
a<-a[!is.na(a)]

Enrichment_DO <- enrichDO(gene = a,
                          ont = "HDO",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          minGSSize = 5,
                          maxGSSize = 10000,
                          readable = TRUE)
DO_enrich = data.frame(Enrichment_DO@result)
DO_enrich10<-DO_enrich[DO_enrich$qvalue<0.05,]
write.csv(DO_enrich10,"DO_enrich.Hominoidea.csv")

p1 <- ggplot(DO_enrich10, aes(x =reorder(Description,log10(qvalue)), y = -log10(qvalue)))+
  geom_bar(stat = "identity", fill = "#2e276d") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#d82819") +
  labs(y = "-log10(q)", title = "",x="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

p1
pdf("DOnrich.pdf",width =2.5,height =3.5)
print(p1) #函数里面直接添加print()
dev.off()


Enrichment_PO <- enrichDO(gene = a,
                          ont = "HPO",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          minGSSize = 5,
                          maxGSSize = 10000,
                          readable = TRUE)
PO_enrich = data.frame(Enrichment_PO@result)
PO_enrich10<-PO_enrich[PO_enrich$pvalue<0.05,]
write.csv(PO_enrich10,"PO_enrich.Hominoidea.csv")

```
