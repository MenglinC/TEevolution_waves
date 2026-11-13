# A. Detection of TE-high-density genes

## A.1 Definition of TE-high-density genes

### A.1.1 genes & TE,overlap

```shell
python "/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/clade_gene/cCRE_TE/overlap_2k_5k_10k_260201_TE_10kp_gene.v3.py" -a "/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/clade_gene/Catarrhini.bed" -b "/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/clade_gene/cCRE_TE/proteinGene.bed"
```


### A.1.2 summarize the relationship between gene & TE

```shell
workplace:/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/clade_gene/v3/


chmod +x summarize_gene_te.sh

./summarize_gene_te.sh \
    Hominoidea.v3.genes.txt \
    Hominoidea.v3.genes.summary.tsv
```


## A.2 length & count regression

file: res.gene.R

```r
setwd("F://03-10//length")
setwd("F://07-31//summarize0731")
# 1. 加载数据
# 假设文件没有表头，我们手动指定列名
Hominoidea 
Catarrhini
Primates
Simiiformes
non-Primates


# data <- read.table("Simiiformes.v3.genes.count.length.txt", 
#                    header = FALSE, 
#                    #col.names = c("Count", "Gene", "Length"),
#                    stringsAsFactors = FALSE)

data <- read.table("Simiiformes.v3.genes.summary.tsv",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                  )

head(data)
# 2. 数据预处理
# 过滤掉长度为负或非正常的数据（根据你提供的片段，前几行似乎有负值噪声）
# data <- subset(data, Length > 0 & !is.na(Gene))

# 3. 构建回归模型
# 使用 log 转换可以让数据更符合线性分布，并处理异方差性
# 我们用 log(Count + 1) 是为了防止 Count 为 0 时产生无穷值
data$gene_length<-as.numeric(data$gene_length)
data <- subset(data, !is.na(gene_length) & !is.na(gene_name)) #上一步转换的时候产生NA是因为原来的数据中，有的基因名对应多个基因组坐标（？）这个是文件的问题，本身是不对的；

data$row_count<-as.numeric(data$row_count)
fit <- lm(log1p(row_count) ~ log1p(gene_length), data = data,na.action = na.exclude)

# 4. 计算残差
# 残差 = 观测值 - 预测值
# 正残差越大，说明该基因在同等长度下拥有的元件比预期多得越多
data$residual <- residuals(fit)

# 5. 定义“异常多”的标准
# 这里我们采用 Z-score 方法：残差大于 2 倍标准差即视为异常高
threshold <- mean(data$residual) + 2 * sd(data$residual)
data$is_outlier <- data$residual > threshold

# 6. 筛选结果并排序
high_density_genes <- data[data$is_outlier, ]
high_density_genes <- high_density_genes[order(-high_density_genes$residual), ]
write.csv(high_density_genes,"../res.gene/Simiiformes_HDG_0731.csv")

# # 查看前 10 个最异常的基因
# print(head(high_density_genes, 10))
# 
# # 7. 可视化（直观检查）
# library(ggplot2)
# p<-ggplot(data, aes(x = log1p(Length), y = log1p(Count))) +
#   geom_point(aes(color = is_outlier), alpha = 0.5) +
#   geom_smooth(method = "lm", color = "blue") +
#   scale_color_manual(values = c("black", "red")) +
#   labs(title = "Gene Length vs Element Count (Residual Analysis)",
#        x = "log(Length + 1)", y = "log(Count + 1)") +
#   theme_minimal()
# p
# 
# residuals(fit)
# 
# #pdf("ResidualAnalysis_Hominoidea.pdf",width =6,height =6 )
# #p
# #dev.off()

```

## A.3 Vennplot

file:venn_lineage.R

```r
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

```


