library(VennDiagram)
setwd("F://03-10")
setwd("F://07-31//res.gene")
# 1. 定义四个文件的路径（请根据实际文件名修改）
files <- c(
  "non-Primates_HDG_0731.csv",
  "Primates_HDG_0731.csv",
  "Simiiformes_HDG_0731.csv",
  "Catarrhini_HDG_0731.csv",  # 提供的示例文件
  "Hominoidea_HDG_0731.csv"
)

# 2. 读取每个文件，提取第二列基因名（列名为"x"），并去重
gene_lists <- lapply(files, function(f) {
  df <- na.omit(read.csv(f, stringsAsFactors = FALSE))
  unique(df$gene_name)  # 第二列基因名
})

# 3. 为每个基因集命名（用于图例）
names(gene_lists) <- c("non-Primates","Primates","Simiiformes","Catarrhini","Hominoidea")

# 4. 计算四个列表的交集
common_genes <- Reduce(intersect, gene_lists)
cat("交集中的基因数量：", length(common_genes), "\n")
if (length(common_genes) > 0) {
  cat("部分交集基因：", head(common_genes), "\n")
} else {
  cat("交集为空\n")
}

# 5. 绘制四元韦恩图
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = names(gene_lists),
  filename = NULL,                # 直接输出到R图形设备
  fill = c("#812e29",
           "#549357",
           "#e7d968",
           "#8b6a93",
           "#92cbe1"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.05
)

# 显示图形
grid.draw(venn.plot)

# 可选：保存为PNG文件
pdf("venn_lineage_0731.pdf",width =6,height =6 )
grid.draw(venn.plot)
dev.off()


all_genes <- unique(unlist(gene_lists))
cat("总基因数（并集）：", length(all_genes), "\n")
library(dplyr)
# 4. 构建出现矩阵（行=基因，列=集合）
occurrence_matrix <- sapply(gene_lists, function(set) all_genes %in% set)
rownames(occurrence_matrix) <- all_genes
colnames(occurrence_matrix) <- names(gene_lists)

# 5. 转换为数据框，并添加pattern列（每个基因所属的组合）
gene_pattern_df <- as.data.frame(occurrence_matrix) %>%
  mutate(pattern = apply(., 1, function(x) paste(as.integer(x), collapse = "")))

# 6. 按pattern分组，统计每个组合的基因数量
pattern_summary <- gene_pattern_df %>%
  group_by(pattern) %>%
  summarise(count = n(), genes = list(rownames(.))) %>%
  arrange(desc(count))

# 7. 输出所有组合的基因列表到CSV文件
# 创建一个列表，每个元素是pattern对应的基因向量
pattern_genes <- split(rownames(gene_pattern_df), gene_pattern_df$pattern)

# 将每个pattern写入单独的CSV文件（也可以写入一个多工作表的Excel文件，这里简化为一个CSV包含pattern列）
# 为方便查看，生成一个包含pattern和基因名的长格式数据框
pattern_long <- data.frame(
  pattern = rep(names(pattern_genes), times = sapply(pattern_genes, length)),
  gene = unlist(pattern_genes)
)
write.csv(pattern_long, "gene_patterns_lineage_SYMBOL.csv", row.names = FALSE)
cat("基因-模式对应表已保存到 gene_patterns_all.csv\n")

