# S1 CpG % & TE density

```r
#workplace:http://localhost:8888/notebooks/figure1b-hs1.ipynb

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


l1_data_0to1_no0 <- final_matrix[rowSums(final_matrix) > 0,]
l1_data_0to1_no0_row_trans<-as.data.frame(t(l1_data_0to1_no0))

merged_corr <- merge(l1_data_0to1_no0_row_trans, annotation_col_new, by = "row.names")#可以直接根据行名呢！这个是我之前不知道的

library(ggplot2)
library(ggpubr)
library(patchwork)

cols_to_plot <- colnames(merged_corr)[
  c(2,4,7,11,15,3,5,8,12,16,6,9,13,17,10,14,18,19)
]

target_col <- "CpG%"

# 1. 计算全部 Pearson correlations
cor_results <- do.call(rbind, lapply(cols_to_plot, function(x_col) {

  x <- merged_corr[[x_col]]
  y <- merged_corr[[target_col]]

  valid <- complete.cases(x, y)

  test <- cor.test(
    x[valid],
    y[valid],
    method = "pearson"
  )

  data.frame(
    variable = x_col,
    r = unname(test$estimate),
    p = test$p.value
  )
}))

# 2. BH-FDR 校正
cor_results$FDR <- p.adjust(
  cor_results$p,
  method = "BH"
)


summary_df2 <- cor_results %>%
  mutate(
    logFDR = -log10(FDR),
    sig = ifelse(FDR < 0.05, "yes", "no"),
    sig_label = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(desc(R2)) %>%
  mutate(variable = factor(variable, levels = rev(variable)))
p_summary_simple <- ggplot(summary_df2, aes(x = R2, y = variable)) +
  geom_point(
    aes(size = logFDR, color = r, shape = sig),
    alpha = 0.9
  ) +
  geom_text(
    aes(label = sig_label),
    nudge_x = 0.015,
    size = 5
  ) +
  scale_color_gradient2(
    low = "#4575b4",
    mid = "grey70",
    high = "#d73027",
    midpoint = 0,
    name = "Pearson r"
  ) +
  scale_shape_manual(
    values = c("yes" = 16, "no" = 1),
    name = "Significant"
  ) +
  scale_size_continuous(
    range = c(2.5, 7),
    name = expression(-log[10](FDR))
  ) +
  coord_cartesian(
    xlim = c(0, max(summary_df2$R2) * 1.2)
  ) +
  labs(
    x = expression(R^2),
    y = NULL,
    title = "Summary of correlations with CpG%"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_summary_simple
```

