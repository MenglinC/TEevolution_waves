# D. For speices genes with H3K9me3 TEs 

## D.1

```r
#setwd("F://08-30")
#workplace:工作路径：F:\tmp\09-10-fig2收尾
tpm<-read.csv("F://09-10//hNPC_TRIM28.tpm.csv",row.names = 1)
head(tpm)
library(AnnotationDbi)
library(org.Hs.eg.db)

# 去掉 .1 .6 等版本号
ensembl_id <- sub("\\..*$", "", rownames(tpm))

gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_id,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
tpm_final <- data.frame(
  Ensembl_ID = rownames(tpm),
  Gene_symbol = unname(gene_symbol),
  tpm,
  check.names = FALSE
)

head(tpm_final)
tpm_final$meanKO<-rowMeans(tpm_final[, 3:5], na.rm = TRUE)
tpm_final$meanWT<-rowMeans(tpm_final[, 6:8], na.rm = TRUE)

write.csv(tpm_final,"hNPC_TRIM28.tpm.symbol.csv") #这个文件是计算了TRIM28 KO前后tpm的均值
#############################################################
tpm_final<-read.csv("F://09-10//hNPC_TRIM28.tpm.symbol.csv")
############################################################
#这里唯一不明确的是这两个结果是怎么得到的？
#太好了，就是喜欢这样干干净净清清爽爽的感觉
```
## D.2

```shell
工作路径：/mnt/backup02/xxzhang/ESC_dataset/mix-seq_Het/results/TE_loci_anno/
#Hominoidea.v3.genes.txt这个文件是，每一个转座子序列与转座子上的基因的对应关系
#153.expr.genes.txt这个文件是我们前面找到的153个是我们人猿总科的基因，同时又在物种间同源的；
awk 'NR==FNR {genes[$1]=1; next} $4 in genes' 153.expr.genes.txt Hominoidea.v3.genes.txt > Hominoidea.v3.153.expr.genes.txt
awk -v OFS="\t" '
{
    start=$8
    if(start<0) start=0

    id=$4 "|" $7 ":" $8 "-" $9 "|" $10 "|" $11 "|" $12

    print $7,start,$9,id #用的是原始的位置
}' Hominoidea.v3.153.expr.genes.txt > TE_flank_0.bed #把交集的TE的bed区域提取出来


bigWigAverageOverBed  "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_Het/results/bowtie2/H3K9me3_FE.bw"  TE_flank_0.bed TE_flank_0.H3K9me3.txt #根据bed计算TE区域平均的bw； 
#bigWigAverageOverBed：Compute average score of big wig over each bed, which may have introns

awk -v OFS="\t" '                                                                                            BEGIN{
    print "gene","TE_chr","TE_start","TE_end","TE_name","TE_class","clade",
          "region_start","region_end","region_size","covered","sum","mean0","mean"
}
{
    split($1,a,"|")
    split(a[2],b,":")
    split(b[2],c,"-")

    te_start=c[1]
    te_end=c[2]

    region_start=te_start-100 #这里不太对吧？这里似乎还是在TE的位置上前后取了100左右，但是这个已经没意义了，因为我们前面计算平均bw的时候用的是原始的peak
    if(region_start<0) region_start=0

    print a[1],b[1],te_start,te_end,a[3],a[4],a[5],
          region_start,te_end+100,$2,$3,$4,$5,$6
}' TE_flank_0.H3K9me3.txt >TE_flank_0.H3K9me3.summary.txt

```

## D.3

```r
epiDat_v1<-read.table("F://09-10//TE_flank_0.H3K9me3.summary.txt",header = T) #分别读入两个文件，一个是H3K9me3的
epiDat_v2<-read.table("F://09-10//TE_flank_0.TRIM28.summary.txt",header = T)#一个是trim28的peak，印象中是计算了chip-seq的peak在其上的最大值还是什么的？==> 需要再去回溯一下

epiDat<-cbind(epiDat_v1,epiDat_v2$mean)


#就以1为分界线，超过1的位点设为het位点，其余为no_het
epiDat$het <- ifelse(epiDat$mean > 1&epiDat$`epiDat_v2$mean`>1, "het", "no_het") #这里做了筛选，要求同时有trim28和H3K9me3的位点为甲基化位点


library(dplyr)

epiDat <- epiDat %>%
  group_by(gene) %>%
  mutate(
    total_n   = n(),                  # 该 gene 的总记录数
    het_n     = sum(het == "het"),    # het == "het" 的记录数
    het_ratio = het_n / total_n       # het 比例 #这个是计算的het-ratio
  ) %>%
  ungroup()
#这里算是计算了每一个基因上甲基化位点的数目

epiDat_2<-unique(epiDat[,c(1,17,19)])

epiDat_3<-merge(epiDat_2,tpm_final,by.x="gene",by.y="Gene_symbol")

head(epiDat_3)

geneDown<-read.table("F://tmp//07-31//res.gene//down_both_chimp_human.txt")


################################################
#从这里开始
################################################


library(dplyr)

dat_enrich <- epiDat_3 %>%
  mutate(
    # 是否属于 chimp + human 相比 rhesus 均下调的基因
    species_down = gene %in% geneDown$symbol,#这里其实有一个true、false的判断
    
    # 是否至少存在一个 H3K9me3+ hominoid-specific TE
    has_het_TE = het_ratio > 0, #其上存在甲基化的比例
    
    # TRIM28 KO effect
    KO_effect = log2(meanKO + 1) - log2(meanWT + 1), #计算KO effect
    
    # TRIM28 KO 后是否上调
    KO_up = KO_effect > 0,
    
    # 符合你的机制模型
    mechanism_pos = has_het_TE & KO_up
  )


dat_het <- dat_enrich %>%
  filter(has_het_TE) #这里筛选掉了其上不存在甲基化的转座子位点的基因

library(dplyr)
library(ggplot2)
library(ggpubr)

# 整理分组名称
dat_het_plot <- dat_het %>%
  mutate(
    group = ifelse(
      species_down,#true
      "Species-down",
      "Background"
    ),
    group = factor(
      group,
      levels = c("Background", "Species-down")
    )
  )

# 单侧 Wilcoxon:
# H1 = Species-down 的 KO_effect > Background
wilcox_res <- with(
  dat_het_plot,
  wilcox.test(
    KO_effect[group == "Species-down"],
    KO_effect[group == "Background"],
    alternative = "greater"
  )
)

wilcox_res

p <- ggplot(
  dat_het_plot,
  aes(x = group, y = KO_effect, fill = group)
) +
  
  # boxplot
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.65,
    linewidth = 0.7
  ) +
  
  # 每个gene
  geom_jitter(
    aes(color = group),
    width = 0.12,
    size = 1.8,
    alpha = 0.55
  ) +
  
  # KO effect = 0
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey40"
  ) +
  
  # Wilcoxon显著性
  stat_pvalue_manual(
    data.frame(
      group1 = "Background",
      group2 = "Species-down",
      y.position = max(dat_het_plot$KO_effect, na.rm = TRUE) * 1.08,
      p = wilcox_res$p.value,
      label = paste0(
        "Wilcoxon, P = ",
        format.pval(
          wilcox_res$p.value, #在这里加上wilcox的p-value #需要有一个xlsx表，来记录这篇文章中用到的统计学检验的结果的指标
          digits = 2,
          eps = 0.001
        )
      )
    ),
    label = "label",
    tip.length = 0.01,
    bracket.size = 0.6,
    size = 4.2
  ) +
  
  scale_fill_manual(
    values = c(
      "Background" = "#6066af",
      "Species-down" = "#941a1d"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Background" = "#6066af",
      "Species-down" = "#941a1d"
    )
  ) +
  
  labs(
    x = NULL,
    y = expression(
      "TRIM28 KO effect  [" *
        log[2] * "(KO+1) - " *
        log[2] * "(WT+1)]"
    )
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(
      size = 13,
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      size = 12,
     # face = "bold"
    ),
    axis.text.y = element_text(size = 11),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.6),
    plot.margin = margin(10, 15, 10, 10)
  )

p

n_bg <- sum(!dat_het_plot$species_down) #加一下数值
n_down <- sum(dat_het_plot$species_down)

p +
  scale_x_discrete(
    labels = c(
      "Background" =
        paste0("Background\n(n = ", n_bg, ")"),
      "Species-down" =
        paste0("Species-down\n(n = ", n_down, ")")
    )
  )+
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "white"
  )

candidate_gene <- dat_het %>%
  filter(species_down) %>%
  arrange(desc(KO_effect)) %>% #按照KO effect降序排列
  select(gene, het_ratio, meanWT, meanKO, KO_effect)

head(candidate_gene, 20)
#NRXN3居首位

setwd("F://09-10//")
write.csv(candidate_gene,"candidate_gene.KO.csv")
write.csv(epiDat,"epiDat_het.csv")

```

