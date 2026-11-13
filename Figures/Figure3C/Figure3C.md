# C

```r
setwd("/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/HM_global_0122/")
data<-read.table("/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/HM_global_0122/all.giggle.result",sep = "\t",fill=TRUE,na.strings = "",quote="",header=T)
library(tidyr)
library(dplyr)
data1<- data %>%
  separate(group, 
           into = c("clade", "TE"), 
           sep = "_(?=[^_]*$)",  # 正则表达式：匹配最后一个下划线
           remove = FALSE,      # 保留原列
           fill = "right")      # 如果拆分不出两部分，右边填充NA
data2<- data1 %>%
  separate(file, 
           into = c("tissue", "state"), 
           sep = "_(?=[^_]*$)",  # 正则表达式：匹配最后一个下划线
           remove = FALSE,      # 保留原列
           fill = "right")      # 如果拆分不出两部分，右边填充NA
tissue_annot<-read.table("/home/xxzhang/workplace/software/giggle/examples/rme/new_groups.txt",sep="\t")
tissue_annot2<-read.table("/home/xxzhang/workplace/software/giggle/examples/rme/EDACC_NAME.txt",sep="\t")
tissue_merge<-merge(tissue_annot,tissue_annot2,by="V1")
colnames(tissue_merge)[3]<-"tissue"
colnames(tissue_merge)[2]<-"tissue_type"
data3<-merge(data2,tissue_merge,by="tissue")
df_summary1<-data3
#df_summary2<-df_summary1[df_summary1$TE%in%c("Alu","ERVK","ERV1","SVA","L1")&df_summary1$clade%in%c("Hominoidea"),] #筛选
#df_summary2<-df_summary1[df_summary1$TE%in%c("ERVK")&df_summary1$clade%in%c("Hominoidea"),] #筛选
df_summary2<-df_summary1[df_summary1$TE%in%c("Alu","ERVK","ERV1","SVA","L1"),] #筛选
df_summary2$TE<-factor(df_summary2$TE,levels=c("L1","ERV1","Alu","ERVK","SVA"))
df_summary2$clade<-factor(df_summary2$clade,levels=c("non-Primates","Primates","Simiiformes","Catarrhini","Hominoidea"))
#df_summary2$state<-factor(df_summary2$state,levels=c("TssA","TssAFlnk","TxFlnk","Tx","TxWk",
 #                                                   "EnhG","Enh","ZNF/Rpts","Het","TssBiv",
  #                                                  "BivFlnk","EnhBiv","ReprPC","ReprPCWk","Quies"))
df_summary2$state<-factor(df_summary2$state,levels=c("Active.TSS","Flanking.Active.TSS","Transcr.at.gene.5.and.3","Strong.transcription","Weak.transcription",
                                                    "Genic.enhancers","Enhancers","ZNF.genes.and.repeats","Heterochromatin","Bivalent.Poised.TSS",
                                                    "Flanking.Bivalent.TSS.Enh","Bivalent.Enhancer","Repressed.PolyComb","Weak.Repressed.PolyComb","Quiescent.Low"))

#df_summary2$group<-factor(df_summary2$group,levels=c("Hominoidea_L1","Hominoidea_ERV1",
#                                                     "Hominoidea_Alu","Hominoidea_ERVK","Hominoidea_SVA"))
library(ggplot2)
library(dplyr)
library(ggbeeswarm)

# 保留原有 chromatin state 的排列顺序
state_levels <- if (is.factor(df_summary2$state)) {
  levels(df_summary2$state)
} else {
  unique(as.character(df_summary2$state))
}

# 1. FDR校正（在筛选显著结果之前）
df_plot <- df_summary2 %>%
  ungroup() %>%
  mutate(
    state = factor(state, levels = state_levels),
    FDR = p.adjust(fishers_right_tail, method = "BH"),
    logFDR = -log10(pmax(FDR, .Machine$double.xmin))
  )

# 2. 只保留显著富集的结果
df_sig <- df_plot %>%
  filter(FDR < 0.05)

# 3. 保留原有配色
tissue_cols <- c(
  "Blood & T-cell" = "#9bcfc4",
  "Brain" = "#f6ad35",
  "Cancer Cell Line" = "#64599d",
  "Digestive" = "#c7daeb",
  "Epithelial" = "#cce5c1",
  "ESC" = "#cc4e96",
  "ES-deriv" = "#b68141",
  "Heart" = "#f3e6c6",
  "HSC & B-cell" = "#845120",
  "iPSC" = "#f0e442",
  "Lung" = "#c76756",
  "Mesench" = "#5d1515",
  "Muscle" = "#f4d0e4",
  "Neurosph" = "#afd063",
  "Other" = "#ee9d9c",
  "Sm Muscle" = "#175595",
  "Thymus" = "#fdf4c1"
)

# 4. 绘图
p <- ggplot(
  df_sig,
  aes(
    x = state,
    y = logFDR,
    fill = tissue_type,
    group = tissue_type
  )
) +
  geom_quasirandom(
  width = 0.2,
  dodge.width = 0.7,
  shape = 21,
  size = 1.8,
  stroke = 0.15,
  color = "grey25",
  alpha = 0.9
)+
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dotted",
    color = "grey50"
  ) +
  # geom_vline(
  #   xintercept = 7.5,
  #   linetype = "dashed",
  #   color = "grey50",
  #   linewidth = 0.4
  # ) +
  scale_fill_manual(values = tissue_cols, drop = FALSE) +
  scale_x_discrete(
    drop = FALSE,
    labels = function(x) stringr::str_wrap(x, width = 10)
  ) +
  facet_grid(TE ~ clade) +
  labs(
    x = "Chromatin state",
    y = expression(-log[10]("FDR")),
    fill = "Tissue type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 12
    ),
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    panel.spacing = unit(0.8, "lines"),
    legend.text = element_text(size = 12)
  )+
theme(
    # 强调图主体
    panel.spacing = unit(0.4, "lines"),

    # 分面标签
    strip.text = element_text(size = 12, face = "bold"),

    # 坐标轴
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),

    # 精简图例
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.35, "cm"),

    # 减少外围留白
    plot.margin = margin(5, 5, 5, 5)
  ) +
theme(
    panel.background = element_rect(
      fill = "#F5F5F5", color = NA
    )#,
    # plot.background = element_rect(
    #   fill = "#F5F5F5", color = NA
    # ),
    # strip.background = element_rect(
    #   fill = "#E8E8E8", color = NA
    # ),
    # legend.background = element_rect(
    #   fill = "#F5F5F5", color = NA
    # )
  )+
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE)
  )

# # 5. 导出大尺寸矢量图
# ggsave(
#   "TE_enrichment_FDR.pdf",
#   plot = p,
#   width = 20,
#   height = 12
# )

pdf("quasirandom_1008.pdf",width=9,height=6 )
p
dev.off()


```
