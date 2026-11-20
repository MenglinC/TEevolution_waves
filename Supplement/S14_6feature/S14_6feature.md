# S14

```r
#file: F://11-17-fig3//nmf3-261117.R
#site_assignment来源于NMF中的W矩阵，这样我们好知道不同的loci在这6个特征中的比例
site_assignment <- apply(W_adj, 1, which.max)
#write.csv(site_assignment,"site_assignment_0317.csv")
write.csv(site_assignment,"site_assignment_0814.csv") #要及时保存这个，做的真的很崩溃很多个大文件
write.csv(site_assignment,"site_assignment_1008.csv") #要及时保存这个，做的真的很崩溃很多个大文件，做到这一步的时候保存，然后重启


library(ggplot2)
library(dplyr)

# 1. 根据 W_adj 最大权重分配类别
site_assignment <- read.csv("F://10-08-画图//site_assignment_1008.csv")
colnames(site_assignment)<-c("Interval","group")
# 2. 将分类结果加入 data_s
# 前提：W_adj 的行与 data_s 的行一一对应且顺序一致
stopifnot(length(site_assignment$Interval) == nrow(mat)) #没有error message

data_s<-merge(site_assignment,annotation)
library(tidyr)
library(dplyr)
data_s <- data_s %>%
   unite("lineage_family", lineage, family, sep = "_", remove = FALSE)


# 3. 计算每个 group 内各类别的比例,还缺age-family的组合情况
pie_data <- data_s %>%
  count(lineage_family, group,lineage, family) %>%
  group_by(lineage_family) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


pie_data$group<-as.character(pie_data$group)

# 1. 指定 lineage_family 顺序
lineage_order <- c(
  "non-Primates_L1",
  "Primates_L1",
  "Simiiformes_L1",
  "Catarrhini_L1",
  "Hominoidea_L1",
  "non-Primates_ERV1",
  "Primates_ERV1",
  "Simiiformes_ERV1",
  "Catarrhini_ERV1",
  "Hominoidea_ERV1",
  "Primates_Alu",
  "Simiiformes_Alu",
  "Catarrhini_Alu",
  "Hominoidea_Alu",
  "Simiiformes_ERVK",
  "Catarrhini_ERVK",
  "Hominoidea_ERVK",
  "Hominoidea_SVA"
)

# 2. 设置因子顺序
pie_data$lineage_family <- factor(
  pie_data$lineage_family,
  levels = lineage_order
)
pie_data$lineage <- factor(
  pie_data$lineage,
  levels = c("non-Primates","Primates","Simiiformes","Catarrhini","Hominoidea")
)

# 2. 设置因子顺序
pie_data$family <- factor(
  pie_data$family,
  levels = c("L1","ERV1","Alu","ERVK","SVA")
)

# 1. 固定 group 顺序
pie_data <- pie_data %>%
  mutate(group = factor(group, levels = as.character(1:6)))

# 2. 计算 2-3 和 5-6 的分界位置
boundary_data <- pie_data %>%
  group_by(lineage, family) %>%
  summarise(
    boundary_23 = sum(prop[as.character(group) %in% c("1", "2")]),
    boundary_56 = sum(prop[as.character(group) %in% as.character(1:5)]),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("boundary_"),
    names_to = "boundary",
    values_to = "y"
  )



# 需要加背景的行
bg_data <- expand.grid(
  family = factor(c("L1", "Alu", "SVA","ERV1","ERVK"),
                  levels = c("L1", "ERV1", "Alu", "ERVK", "SVA")),
  age = factor(c("non-Primates", "Primates", "Simiiformes", "Catarrhini", "Hominoidea"),
               levels = c("non-Primates", "Primates", "Simiiformes", "Catarrhini", "Hominoidea"))
)


# 3. 绘图
pie_plot <- ggplot(
  pie_data,
  aes(x = 1, y = prop, fill = group)
) +
  # 1. 给指定行添加浅灰背景
  geom_rect(
    data = bg_data,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#F2F2F2",
    color = NA
  ) +
  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.25,
    position = position_stack(reverse = TRUE)
  ) +
  
  # 强调 2-3、5-6 分界线
  geom_segment(
    data = boundary_data,
    aes(
      x = 0.5,
      xend = 1.5,
      y = y,
      yend = y
    ),
    inherit.aes = FALSE,
    color = "#444444",
    linewidth = 1.2,
    linetype = "dotted"
  ) +
  
  coord_polar(theta = "y", start = 0) +
  scale_x_continuous(
    limits = c(0.5, 1.5),
    expand = c(0, 0)
  ) +
  facet_grid(family ~ lineage, drop = FALSE) +
  scale_fill_manual(
    values = c(
      "1" = "#9d435d",
      "2" = "#cca6ad",
      "3" = "#8a8b8d",
      "4" = "#88a4cc",
      "5" = "#a39dcd",
      "6" = "#c58224"
    ),
    labels = c(
      "1" = "Strong.transcription.like",
      "2" = "Weak.transcription.like",
      "3" = "Quiescent.low.like",
      "4" = "PolyComb.like",
      "5" = "Heterochromatin.like",
      "6" = "cCRE.like"
    ),
    drop = FALSE
  ) +
  theme_void() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(0.6, "lines"),
    legend.title = element_blank(),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

print(pie_plot)

ggsave(
  "pie_plot.tiff",
  plot = pie_plot,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

```
