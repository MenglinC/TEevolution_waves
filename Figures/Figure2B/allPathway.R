# 安装包（如果未安装）
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE"))

# 加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

#gene1 <- read.csv("../gene_patterns_lineage_SYMBOL.csv", header = TRUE)#1000
#gene2 <- read.csv("../gene_patterns_lineage_SYMBOL.csv", header = TRUE) #100
#gene3 <- read.csv("../gene_patterns_lineage_SYMBOL.csv", header = TRUE) #10
#gene4 <- read.csv("../gene_patterns_lineage_SYMBOL.csv", header = TRUE) #1

setwd("F://07-31//res.gene")
# gene <- read.csv("gene_patterns_lineage_SYMBOL_0731.csv", header = TRUE)#1000
# 
# genes0 <- grep("[-|.]", gene[gene$pattern==10000,]$gene, invert = TRUE, value = TRUE)#174 非人灵长类
# genes1 <- grep("[-|.]", gene[gene$pattern==1000,]$gene, invert = TRUE, value = TRUE) #60 灵长类
# genes2 <- grep("[-|.]", gene[gene$pattern==100,]$gene, invert = TRUE, value = TRUE) #224 类人猿下目
# genes3 <- grep("[-|.]", gene[gene$pattern==10,]$gene, invert = TRUE, value = TRUE)#177 狭鼻类
# genes4 <- grep("[-|.]", gene[gene$pattern==1,]$gene, invert = TRUE, value = TRUE) #221 人科

genes0<-read.csv("non-Primates_HDG_0731.csv")$gene_name
genes1<-read.csv("Primates_HDG_0731.csv")$gene_name
genes2<-read.csv("Simiiformes_HDG_0731.csv")$gene_name
genes3<-read.csv("Catarrhini_HDG_0731.csv")$gene_name
genes4<-read.csv("Hominoidea_HDG_0731.csv")$gene_name


#genes4<-unique(read.csv("Hominoidea_HDG_0731.csv")$species)


# 将基因Symbol转换为ENTREZ ID（富集分析常用）
convert2entrez <- function(symbols) {
  bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}

geneList0 <- convert2entrez(genes0)
geneList1 <- convert2entrez(genes1)
geneList2 <- convert2entrez(genes2)
geneList3 <- convert2entrez(genes3)
geneList4 <- convert2entrez(genes4)

# 定义富集函数
run_enrich <- function(genes, ont = "BP") {
  enrichGO(gene = genes,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = ont,
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)  # 将ENTREZID转为Symbol
}

# 分别进行富集
ego0 <- run_enrich(geneList0)
ego1 <- run_enrich(geneList1)
ego2 <- run_enrich(geneList2)
ego3 <- run_enrich(geneList3)
ego4 <- run_enrich(geneList4)

# 提取结果（若某组没有显著富集，可用空数据框占位）
extract_df <- function(ego, group_name) {
  if (is.null(ego) || nrow(ego@result) == 0) {
    return(data.frame())
  }
  df <- ego@result
  df$Group <- group_name
  return(df)
}

df0 <- extract_df(ego0, "non-Primates")
df1 <- extract_df(ego1, "Primates")
df2 <- extract_df(ego2, "Simiiformes")
df3 <- extract_df(ego3, "Catarrhini")
df4 <- extract_df(ego4, "Hominoidea")
#write.csv(df4,"GOpathway.Hominoidea_153genes.HDG.0804.csv")

# 合并所有数据框
df_all <- bind_rows(df0,df1, df2, df3, df4)
df_all$Group<-factor(df_all$Group,levels = c("non-Primates","Primates","Simiiformes","Catarrhini","Hominoidea"))
df_all2<-df_all[df_all$pvalue<0.05,]
write.csv(df_all,"GOpathway.all.HDG.0731.csv") #已经保存了；

df_all2<-read.csv("GOpathway.all.HDG.0731.csv")
df_all2$Group<-factor(df_all2$Group,levels = c("non-Primates","Primates","Simiiformes","Catarrhini","Hominoidea"))

# 按p.adjust排序，每组取前5个
topN <- 4
df_top <- df_all2 %>%
  group_by(Group) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = topN) %>%
  ungroup()

# # 为了绘图美观，可将Description转为因子并按p值排序
# df_top$Description <- factor(df_top$Description,
#                              levels = unique(df_top$Description[order(df_top$p.adjust)]))
df_top$Description <- factor(df_top$Description,
                             levels = unique(df_top$Description[order(df_top$Group)]))


# p <- ggplot(df_top, aes(x = Group, y = Description)) +
#   geom_point(aes(size = Count, color = pvalue)) +
#   scale_color_gradient(low = "red", high = "blue", name = "Pvalue") +
#   scale_size_continuous(range = c(3, 8), name = "Gene count") +
#   theme_bw(base_size = 14) +
#   coord_flip()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text.y = element_text(size = 12),
#         legend.position = "right") +
#   labs(x = "Gene set", y = "GO term", 
#        title = "GO enrichment across four groups")


library(ggplot2)
library(dplyr)
library(forcats)

# GO term排序（按照平均显著性）
df_top <- df_top %>%
  mutate(
    Description = fct_reorder(Description, -log10(pvalue))
  )

# p <- ggplot(df_top, 
#             aes(x = Group, 
#                 y = Description)) +
#   
#   geom_point(
#     aes(size = Count, 
#         color = -log10(pvalue)),
#     alpha = 0.85
#   ) +
#   
#   # ASPAR/Nature常用蓝紫色渐变风格
#   scale_color_gradientn(
#     colors = c("#440154", "#31688E", "#35B779", "#FDE725"),
#     name = "-log10(P)"
#   ) +
#   
#   scale_size_continuous(
#     range = c(3, 10),
#     name = "Gene count"
#   ) +
#   
#   theme_classic(base_size = 14) +
#   
#   coord_flip() +
#   
#   theme(
#     axis.text.y = element_text(
#       size = 11,
#       color = "black"
#     ),
#     axis.text.x = element_text(
#       angle = 45,
#       hjust = 1,
#       size = 11
#     ),
#     
#     axis.title = element_text(
#       size = 14,
#       face = "bold"
#     ),
#     
#     legend.title = element_text(
#       size = 12,
#       face = "bold"
#     ),
#     
#     legend.text = element_text(
#       size = 10
#     ),
#     
#     panel.grid.major.y = element_line(
#       color = "grey90",
#       linewidth = 0.3
#     ),
#     
#     panel.grid.major.x = element_blank(),
#     
#     plot.title = element_text(
#       size = 16,
#       face = "bold",
#       hjust = 0.5
#     ),
#     
#     legend.position = "right"
#   ) +
#   
#   labs(
#     x = NULL,
#     y = NULL,
#     title = "GO enrichment across four groups"
#   )



library(ggplot2)
library(dplyr)
library(forcats)

# # GO term排序##并不需要
df_top <- df_top %>%
  mutate(
    logQ = -log10(qvalue)#,
    #Description = fct_reorder(Description, logP)
  )


p <- ggplot(
  df_top,
  aes(
    x = Group,
    y = Description
  )
) +
  
  # bubble
  geom_point(
    aes(
      size = Count,
      color = logQ
    ),
    alpha = 0.85
  ) +
  # 
  # # apEAR风格渐变色
  # scale_color_gradientn(
  #   colors = c(
  #     "#3B4CC0",
  #     "#6BAED6",
  #     "#F7F7F7",
  #     "#FCAE91",
  #     "#B2182B"
  #   ),
  #   name = "-log10(P)"
  # ) +
  # scale_color_gradientn(
  #   colors = c("#4575B4", "#91BFDB", "#FFFFBF", "#FC8D59", "#D73027"),
  #   name = "-log10(P)"
  # )+
  # ASPAR/Nature常用蓝紫色渐变风格
  scale_color_gradientn(
    colors = c("#440154", "#31688E", "#35B779", "#FDE725"),
    name = "-log10(q-value)"
  ) +
  # 点大小
  scale_size_continuous(
    range = c(3, 9),
    name = "Gene count"
  ) +
  
  coord_flip() +
  
  theme_classic(base_size = 14) +
  
  theme(
    
    # y轴GO term
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    
    # group标签
    axis.text.x = element_text(
      size = 11,
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    # legend
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    # 网格线
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    
    panel.grid.major.x = element_blank(),
    
    # title
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.margin = margin(
      10,10,10,10
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title = "GO enrichment across four groups"
  )


p


print(p)
pdf("dotplot_pvalue_allHDG_v4.pdf",width =9,height =5)
#pdf(paste("pheatmap_sinh","Primates","pdf",sep = "."),width =4 ,height = 4)
print(p)
dev.off()

