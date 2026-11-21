# S5

```r
#file:figure1b-hs1-jaccard.ipynb
setwd("/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/")
library(data.table)

# 使用fread读取，速度比read.csv快5-10倍
data <- fread("te_density_matrix.csv", header = TRUE)
# 如果需要第一列作为行名
row_names <- data[[1]]
data <- as.matrix(data[, -1])
rownames(data) <- row_names
library(data.table)
library(pheatmap)
library(dplyr)
library(tidyr)  # pivot_longer和pivot_wider在tidyr包中

# 使用基础R的rownames处理
data_tibble <- as.data.frame(data)
data_tibble$TE_Family <- rownames(data_tibble)

# 将行名列移到第一列
data_tibble <- data_tibble[, c("TE_Family", setdiff(colnames(data_tibble), "TE_Family"))]

# 继续后续操作
long_format <- data_tibble %>%
  pivot_longer(cols = -TE_Family, 
               names_to = "Sample", 
               values_to = "Expression")
# 从样本名中提取类型信息
long_format <- long_format %>%
  mutate(Type = case_when(
    grepl("_L1$", Sample) ~ "L1",
    grepl("_ERVK$", Sample) ~ "ERVK",
    grepl("_Alu$", Sample) ~ "Alu",
    grepl("_ERV1$", Sample) ~ "ERV1",
    grepl("_SVA$", Sample) ~ "SVA",
    TRUE ~ "Other"
  )) %>%
  mutate(Sample_Base = gsub("_(L1|ERVK|Alu|ERV1|SVA)$", "", Sample))

# 创建新的行标识符
long_format <- long_format %>%
  mutate(TE_Family_Type = paste0(TE_Family, "_", Type))

# 转换回宽格式（行是TE_Family_Type，列是样本）
wide_format <- long_format %>%
  select(TE_Family_Type, Sample_Base, Expression) %>%
  pivot_wider(names_from = Sample_Base, values_from = Expression)

# 转换为矩阵用于热图
final_matrix <- as.matrix(wide_format[, -1])
rownames(final_matrix) <- wide_format$TE_Family_Type

#scale 0-1,分subfamily scale比较对
l1_data_0to1 <- t(apply(final_matrix, 1, function(row) {
  if (max(row) == min(row)) {
    # 如果所有值相同，返回0
    return(rep(0, length(row))) #我这里其实是做了修改的，如果所有值相同，则返回0
  } else {
    # 缩放到0-1范围
    return((row - min(row)) / (max(row) - min(row)))
  }
}))

l1_data_0to1<-as.data.frame(l1_data_0to1) #这个就是得到的scale后的矩阵(这一步其实可选)
colnames(l1_data_0to1)<-colnames(final_matrix)
dim(l1_data_0to1)
l1_data_0to1_no0_row <- l1_data_0to1[rowSums(l1_data_0to1) > 0,] #这是做了修改，只对行做筛选
dim(l1_data_0to1_no0_row)

#############################################################################################
library(viridis)
# 1. 提取每行top10的列名（位置）
get_top_cols <- function(row, percent = 0.01) {
  # 计算要选取的列数：总列数 × percent，向上取整（保证至少1列）
  k <- max(1, ceiling(length(row) * percent))
  # 按值降序排序，取前k个列名
  sorted_cols <- names(sort(row, decreasing = TRUE))[1:k]
  return(sorted_cols)
}

get_bottom_cols <- function(row, percent = 0.02) {
  # 计算要选取的列数
  k <- max(1, ceiling(length(row) * percent))

  # 按值升序排序，取最小的k个列名
  sorted_cols <- names(sort(row, decreasing = FALSE))[1:k]

  return(sorted_cols)
}


# 应用函数到每一行，得到一个列表，每个元素是该基因的top10位置集合
#gene_sets <- apply(l1_data_0to1_no0_row, 1, get_top_cols,simplify = FALSE)
gene_sets <- apply(l1_data_0to1_no0_row, 1, get_top_cols,simplify = FALSE)
names(gene_sets) <- rownames(l1_data_0to1_no0_row)  # 确保基因名作为列表名

# 2. 计算任意两个基因的Jaccard系数
jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# 初始化Jaccard矩阵
n_genes <- nrow(l1_data_0to1_no0_row)
jaccard_mat <- matrix(0, nrow = n_genes, ncol = n_genes)
rownames(jaccard_mat) <- colnames(jaccard_mat) <- rownames(l1_data_0to1_no0_row)

# 双重循环填充矩阵（只计算上三角，然后对称赋值）
for (i in 1:n_genes) {
  for (j in i:n_genes) {
    jaccard_mat[i, j] <- jaccard(gene_sets[[i]], gene_sets[[j]])
    jaccard_mat[j, i] <- jaccard_mat[i, j]
  }
}

jaccard_dist <- as.dist(
    1 - jaccard_mat
  )
hc <- hclust(
    jaccard_dist,
    method = "average"
  )
cluster_result <- cutree(
      hc,
      k = 3
    )
cluster_result

jac001<-pheatmap(jaccard_mat,color = viridis(8),
    cluster_rows = hc,
    cluster_cols = hc
         #clustering_method ="average"
              )

# 1. 提取每行top10的列名（位置）
get_top_cols <- function(row, percent = 0.02) {
  # 计算要选取的列数：总列数 × percent，向上取整（保证至少1列）
  k <- max(1, ceiling(length(row) * percent))
  # 按值降序排序，取前k个列名
  sorted_cols <- names(sort(row, decreasing = TRUE))[1:k]
  return(sorted_cols)
}

get_bottom_cols <- function(row, percent = 0.02) {
  # 计算要选取的列数
  k <- max(1, ceiling(length(row) * percent))

  # 按值升序排序，取最小的k个列名
  sorted_cols <- names(sort(row, decreasing = FALSE))[1:k]

  return(sorted_cols)
}


# 应用函数到每一行，得到一个列表，每个元素是该基因的top10位置集合
#gene_sets <- apply(l1_data_0to1_no0_row, 1, get_top_cols,simplify = FALSE)
gene_sets <- apply(l1_data_0to1_no0_row, 1, get_top_cols,simplify = FALSE)
names(gene_sets) <- rownames(l1_data_0to1_no0_row)  # 确保基因名作为列表名

# 2. 计算任意两个基因的Jaccard系数
jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# 初始化Jaccard矩阵
n_genes <- nrow(l1_data_0to1_no0_row)
jaccard_mat <- matrix(0, nrow = n_genes, ncol = n_genes)
rownames(jaccard_mat) <- colnames(jaccard_mat) <- rownames(l1_data_0to1_no0_row)

# 双重循环填充矩阵（只计算上三角，然后对称赋值）
for (i in 1:n_genes) {
  for (j in i:n_genes) {
    jaccard_mat[i, j] <- jaccard(gene_sets[[i]], gene_sets[[j]])
    jaccard_mat[j, i] <- jaccard_mat[i, j]
  }
}

jaccard_dist <- as.dist(
    1 - jaccard_mat
  )
hc <- hclust(
    jaccard_dist,
    method = "average"
  )
cluster_result <- cutree(
      hc,
      k = 3
    )
cluster_result

jac002<-pheatmap(jaccard_mat,color = viridis(8),
    cluster_rows = hc,
    cluster_cols = hc
         #clustering_method ="average"
              )

# 1. 提取每行top10的列名（位置）
get_top_cols <- function(row, percent = 0.03) {
  # 计算要选取的列数：总列数 × percent，向上取整（保证至少1列）
  k <- max(1, ceiling(length(row) * percent))
  # 按值降序排序，取前k个列名
  sorted_cols <- names(sort(row, decreasing = TRUE))[1:k]
  return(sorted_cols)
}

get_bottom_cols <- function(row, percent = 0.02) {
  # 计算要选取的列数
  k <- max(1, ceiling(length(row) * percent))

  # 按值升序排序，取最小的k个列名
  sorted_cols <- names(sort(row, decreasing = FALSE))[1:k]

  return(sorted_cols)
}


# 应用函数到每一行，得到一个列表，每个元素是该基因的top10位置集合
#gene_sets <- apply(l1_data_0to1_no0_row, 1, get_top_cols,simplify = FALSE)
gene_sets <- apply(l1_data_0to1_no0_row, 1, get_top_cols,simplify = FALSE)
names(gene_sets) <- rownames(l1_data_0to1_no0_row)  # 确保基因名作为列表名

# 2. 计算任意两个基因的Jaccard系数
jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# 初始化Jaccard矩阵
n_genes <- nrow(l1_data_0to1_no0_row)
jaccard_mat <- matrix(0, nrow = n_genes, ncol = n_genes)
rownames(jaccard_mat) <- colnames(jaccard_mat) <- rownames(l1_data_0to1_no0_row)

# 双重循环填充矩阵（只计算上三角，然后对称赋值）
for (i in 1:n_genes) {
  for (j in i:n_genes) {
    jaccard_mat[i, j] <- jaccard(gene_sets[[i]], gene_sets[[j]])
    jaccard_mat[j, i] <- jaccard_mat[i, j]
  }
}

jaccard_dist <- as.dist(
    1 - jaccard_mat
  )
hc <- hclust(
    jaccard_dist,
    method = "average"
  )
cluster_result <- cutree(
      hc,
      k = 3
    )
cluster_result

jac003<-pheatmap(jaccard_mat,color = viridis(8),
    cluster_rows = hc,
    cluster_cols = hc
         #clustering_method ="average"
              )

library(ggplotify)
library(patchwork)
library(ggplot2)
jac001_gg <- as.ggplot(jac001)
jac002_gg <- as.ggplot(jac002)
jac003_gg <- as.ggplot(jac003)
combined_plot2 <- (
  (jac001_gg + ggtitle("1%")) +
  (jac002_gg + ggtitle("2%")) +
  (jac003_gg + ggtitle("3%"))
) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(combined_plot2)
pdf("jaccard_mat_1121_3plot_jacdistance.pdf",width =18 ,height =6 )
combined_plot2
dev.off()

```
non-Primates_L1
5835
non-Primates_ERV1
5537
Primates_L1
5793
Primates_ERV1
5764
Primates_Alu
5837
Simiiformes_L1
5528
Simiiformes_ERV1
274
Simiiformes_Alu
5832
Simiiformes_ERVK
257
Catarrhini_L1
5099
Catarrhini_ERV1
3887
Catarrhini_Alu
5831
Catarrhini_ERVK
2822
Hominoidea_L1
5200
Hominoidea_ERV1
2104
Hominoidea_Alu
5643
Hominoidea_ERVK
1105
Hominoidea_SVA
2788

这里统计的是非零值的比例。

6125	0.01,1%	0.02,2%	0.05,5%
	61.25	122.5	306.25,这个就已经超了257，余下来的都是空的；不合适
	0.03,3%		
	183.75		

