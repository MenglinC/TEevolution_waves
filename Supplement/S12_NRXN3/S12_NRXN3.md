# S12

```r
setwd("F://11-21-画附图")
data<-read.table("F://tmp//02-02//species_rpkm//simFiltered_rpkm_combat_nonParametric.txt",header = T) 
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)  # 用于统计检验和显著性标注
library(rstatix) # 辅助统计检验

# 处理数据
long_data <- melt(data) #Using ProbeID as id variables

long_data2 <- separate(long_data, col = 1, 
                       into = c("ENSG", "symbol"), 
                       sep = "\\|")

long_data3 <- separate(long_data2, col = 3, 
                       into = c("speicesID", "region"), 
                       sep = "\\.")

long_data3 <- long_data3 %>%
  mutate(species = case_when(
    grepl("HSB", .[[3]], ignore.case = TRUE) ~ "human",
    grepl("PTB", .[[3]], ignore.case = TRUE) ~ "chimp",
    grepl("RMB", .[[3]], ignore.case = TRUE) ~ "rhesus",
    TRUE ~ NA_character_
  ))

# 1. 提取NRXN3表达数据
df <- long_data3 %>%
  filter(
    symbol == "NRXN3",
    region == "S1C", #修改这里
    !is.na(value),
    !is.na(species),
    !is.na(region)
  ) %>%
  group_by(species, region, speicesID) %>%
  summarise(
    value = mean(value),
    .groups = "drop"
  ) %>%
  mutate(
    species = factor(species),
    region = factor(region)
  )

# 查看各组样本量
df %>% count(region, species)

# 2. 筛选满足统计要求的脑区
# 三个物种均存在，且每个物种至少2个样本
valid_regions <- df %>%
  count(region, species) %>%
  group_by(region) %>%
  filter(n() == 3, all(n >= 2)) %>%
  pull(region) %>%
  unique()

test_df <- df %>%
  filter(region %in% valid_regions)

# 3. Kruskal-Wallis检验
kw_results <- test_df %>%
  group_by(region) %>%
  kruskal_test(value ~ species) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p, method = "BH"))

# 4. Dunn两两比较
dunn_results <- test_df %>%
  group_by(region) %>%
  dunn_test(value ~ species, p.adjust.method = "none") %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(p, method = "BH"),
    label = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# 只显示整体检验与两两比较均显著的结果
sig_results <- dunn_results %>%
  left_join(
    kw_results %>%
      select(region, KW_FDR = FDR),
    by = "region"
  ) %>%
  filter(KW_FDR < 0.05, FDR < 0.05)

# 5. 计算显著性标注的位置
region_range <- df %>%
  group_by(region) %>%
  summarise(
    ymax = max(value),
    yrange = max(value) - min(value),
    .groups = "drop"
  ) %>%
  mutate(
    yrange = pmax(yrange, 0.1)
  )

sig_results <- sig_results %>%
  left_join(region_range, by = "region") %>%
  group_by(region) %>%
  arrange(FDR, .by_group = TRUE) %>%
  mutate(
    y.position = ymax + row_number() * 0.15 * yrange
  ) %>%
  ungroup()

df$species<-factor(df$species,levels = c("rhesus","chimp","human"))
# 6. 绘图
p <- ggplot(df,aes(x = species, y = value, fill = species))+
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.6
  ) +
  facet_wrap(
    ~region,
    ncol = 4,
    scales = "free_y"
  ) +
  #scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL,
    y = "RPKM"
  ) +
  scale_fill_manual(
    values = c(
      "human" = "#5a8fbf",
      "chimp" = "#60a764",
      "rhesus" = "#813117"
    ))+
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(
      fill = "grey95",
      color = NA
    ),
    strip.text = element_text(face = "bold")
  )

# 添加显著性
if (nrow(sig_results) > 0) {
  p <- p +
    stat_pvalue_manual(
      sig_results,
      label = "label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.3,
      size = 3
    )
}

print(p)

# 7. 保存结果
ggsave("NRXN3_species_S1C.pdf",
       p, width =3, height = 3)

write.csv(kw_results, "NRXN3_KW_results_S1C.csv",
          row.names = FALSE)

write.csv(dunn_results, "NRXN3_Dunn_results_s1C.csv",
          row.names = FALSE)




```
