# 安装包（如果未安装）
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE"))

# 加载包
#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(ggplot2)
#library(dplyr)

#gene1 <- read.csv("../gene_patterns_lineage_SYMBOL.csv", header = TRUE)#1000
#gene2 <- read.csv("../gene_patterns_lineage_SYMBOL.csv", header = TRUE) #100
#gene3 <- read.csv("../gene_patterns_lineage_SYMBOL.csv", header = TRUE) #10
#gene4 <- read.csv("../gene_patterns_lineage_SYMBOL.csv", header = TRUE) #1

setwd("F://10-19//fig2-reanalysis")
# gene <- read.csv("gene_patterns_lineage_SYMBOL_0731.csv", header = TRUE)#1000
# 
# genes0 <- grep("[-|.]", gene[gene$pattern==10000,]$gene, invert = TRUE, value = TRUE)#174 非人灵长类
# genes1 <- grep("[-|.]", gene[gene$pattern==1000,]$gene, invert = TRUE, value = TRUE) #60 灵长类
# genes2 <- grep("[-|.]", gene[gene$pattern==100,]$gene, invert = TRUE, value = TRUE) #224 类人猿下目
# genes3 <- grep("[-|.]", gene[gene$pattern==10,]$gene, invert = TRUE, value = TRUE)#177 狭鼻类
# genes4 <- grep("[-|.]", gene[gene$pattern==1,]$gene, invert = TRUE, value = TRUE) #221 人科

genes0<-read.csv("F://tmp//07-31-figure2重新画图//res.gene//non-Primates_HDG_0731.csv")$gene_name
genes1<-read.csv("F://tmp//07-31-figure2重新画图//res.gene//Primates_HDG_0731.csv")$gene_name
genes2<-read.csv("F://tmp//07-31-figure2重新画图//res.gene//Simiiformes_HDG_0731.csv")$gene_name
genes3<-read.csv("F://tmp//07-31-figure2重新画图//res.gene//Catarrhini_HDG_0731.csv")$gene_name
genes4<-read.csv("F://tmp//07-31-figure2重新画图//res.gene//Hominoidea_HDG_0731.csv")$gene_name

brain.enriched<-read.table("tissue_category_rna_brain_Tissue-group1.tsv",sep = "\t",header = T)$Gene
brain.group<-read.table("tissue_category_rna_brain_Tissue-group2.tsv",sep = "\t",header = T)$Gene
brain.enhanced<-read.table("tissue_category_rna_brain_Tissue-group3.tsv",sep = "\t",header = T)$Gene

library(dplyr)
library(ggplot2)

# 1. 整理五组 HDG
gene_sets <- list(
  "non-Primates" = genes0,
  "Primates" = genes1,
  "Simiiformes" = genes2,
  "Catarrhini" = genes3,
  "Hominoidea" = genes4
)

# 2. 定义背景基因
# 替换成筛选 HDG 时使用的全部候选基因
background <- read.table("proteinGene.bed")$V4 #20030

background <- unique(na.omit(background))
brain <- intersect(unique(na.omit(brain.group)), background) ##这里

# 3. 计算富集
results <- lapply(names(gene_sets), function(group) {
  
  genes <- intersect(
    unique(na.omit(gene_sets[[group]])),
    background
  )
  
  a <- length(intersect(genes, brain))
  b <- length(genes) - a
  c <- length(brain) - a
  d <- length(background) - a - b - c
  
  mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE
  )
  
  ft <- fisher.test(mat, alternative = "greater")
  
  # Fold Enrichment
  observed <- a / length(genes)
  expected <- length(brain) / length(background)
  fold_enrichment <- observed / expected
  
  data.frame(
    Group = group,
    HDG_total = length(genes),
    Brain_HDG = a,
    Brain_fraction = observed,
    Fold_Enrichment = fold_enrichment,
    Odds_Ratio = unname(ft$estimate),
    P_value = ft$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(P_value, method = "BH")
  )

print(results)


################################################################################


get_OR <- function(HDG, brain) {
  
  a <- sum(HDG & brain)
  b <- sum(HDG & !brain)
  c <- sum(!HDG & brain)
  d <- sum(!HDG & !brain)
  
  counts <- c(a = a, b = b, c = c, d = d)
  
  # Haldane-Anscombe correction
  corrected <- counts + 0.5 #交集可能出0的情况
  
  a2 <- corrected["a"]
  b2 <- corrected["b"]
  c2 <- corrected["c"]
  d2 <- corrected["d"]
  
  logOR <- log((a2 * d2) / (b2 * c2))
  
  influence <- (HDG & brain) / a2 -
    (HDG & !brain) / b2 -
    (!HDG & brain) / c2 +
    (!HDG & !brain) / d2
  
  list(
    logOR = unname(logOR),
    influence = influence,
    counts = counts,
    has_zero = any(counts == 0)
  )
} 


# 直接比较 Hominoidea 与其他四个阶段
compare_hominoidea <- function(brain_genes, category) {
  
  brain <- background %in% unique(brain_genes)
  
  hom <- background %in% gene_sets[["Hominoidea"]]
  ref <- get_OR(hom, brain)
  
  results <- lapply(names(gene_sets)[1:4], function(stage) {
    
    other <- background %in% gene_sets[[stage]]
    test <- get_OR(other, brain)
    
    # log(OR_Hominoidea / OR_other)
    diff <- ref$logOR - test$logOR
    
    # 考虑相同基因在两组中重复出现
    infl_diff <- ref$influence - test$influence
    
    SE <- sqrt(sum((infl_diff - mean(infl_diff))^2))
    
    z <- diff / SE
    p <- if (ref$has_zero || test$has_zero || SE == 0) {
      NA_real_
    } else {
      2 * pnorm(-abs(z))
    } #这里做了一些修改
    
    data.frame(
      Category = category,
      Compared_stage = stage,
      OR_Hominoidea = exp(ref$logOR),
      OR_other = exp(test$logOR),
      OR_ratio = exp(diff),
      CI_lower = exp(diff - 1.96 * SE),
      CI_upper = exp(diff + 1.96 * SE),
      P_value = p
    )
  })
  
  do.call(rbind, results)
}




# 分别检验 HPA 三类基因
res <- rbind(
  compare_hominoidea(brain.enriched, "Tissue enriched"),
  compare_hominoidea(brain.group, "Group enriched"),
  compare_hominoidea(brain.enhanced, "Tissue enhanced")
)

res$FDR <- p.adjust(res$P_value, method = "BH")

print(res)







