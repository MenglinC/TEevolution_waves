# C.jaccard similarity

## C.1.1 main figure-top2%

```r

# 1. 提取每行top10的列名（位置）
## Function:get_top10_cols
get_top10_cols <- function(row, percent = 0.02) {
  # 计算要选取的列数：总列数 × percent，向上取整（保证至少1列）
  k <- max(1, ceiling(length(row) * percent))
  # 按值降序排序，取前k个列名
  sorted_cols <- names(sort(row, decreasing = TRUE))[1:k]
  return(sorted_cols)
}

## 应用函数到每一行，得到一个列表，每个元素是该基因的top10位置集合
gene_sets <- apply(l1_data_0to1_no0_row_final, 1, get_top10_cols,simplify = FALSE)  #input：l1_data_0to1_no0_row_final
names(gene_sets) <- rownames( l1_data_0to1_no0_row_final)  # 确保基因名作为列表名

# 2. 计算任意两个set之间的Jaccard系数
## Function：compute jaccard similarity between two sets
jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}
## 初始化Jaccard矩阵
n_genes <- nrow(l1_data_0to1_no0_row)
jaccard_mat <- matrix(0, nrow = n_genes, ncol = n_genes)
rownames(jaccard_mat) <- colnames(jaccard_mat) <- rownames(l1_data_0to1_no0_row)

## 双重循环填充矩阵（只计算上三角，然后对称赋值）
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

# 3. plot the jaccard matrix among sets
jac002<-pheatmap(jaccard_mat,color = viridis(8),
    cluster_rows = hc,
    cluster_cols = hc
         #clustering_method ="average"
              )

```

## C.1.2 supplementary figures(other)

```r
# 1. repeat for different numbers

## return jac001,002,005

# 2.combine plots

library(ggplotify)
library(patchwork) 
library(ggplot2)

jac001_gg <- as.ggplot(jac001)
jac002_gg <- as.ggplot(jac002)
jac005_gg <- as.ggplot(jac005)


combined_plot2 <- (
  (jac001_gg + ggtitle("0.01")) +
  (jac002_gg + ggtitle("0.02")) +
  (jac005_gg + ggtitle("0.05"))
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(combined_plot2)

```
