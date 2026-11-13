# C. expression across species

#file:species.rpkm_0731.R
#input:simFiltered_rpkm_combat_nonParametric.txt,gene_patterns_lineage_SYMBOL_0731.csv

```r
setwd("F://tmp//02-02//species_rpkm")
data<-read.table("simFiltered_rpkm_combat_nonParametric.txt",header = T) #这个是老师的

setwd("F://tmp//07-31//res.gene")
gene <- read.csv("gene_patterns_lineage_SYMBOL_0731.csv", header = TRUE)#1 
gene_list <- grep("[-|.]", gene[gene$pattern==1,]$gene, invert = TRUE, value = TRUE) #1人科

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


# 筛选axonogenesis基因数据
#axon_data <- long_data3[long_data3$symbol %in% axonogenesis&long_data3$region=="V1C", ]
#不一定要这个通路的，那就太刻意了
axon_data <- long_data3[long_data3$symbol %in% gene_list, ] ###在这里修改


# 设置因子水平
axon_data$species <- factor(axon_data$species, levels = c("rhesus","chimp","human"))
axon_data$region <- factor(axon_data$region,levels=c("MFC","OFC","DFC","VFC","M1C","S1C","IPC","A1C","STC","ITC","V1C",
                                                     "HIP","AMY","STR","MD","CBC"))

# stat_results <- axon_data %>%
#   group_by(region) %>%
#   wilcox_test(value ~ species, 
#               p.adjust.method = "BH",  # Benjamini-Hochberg校正
#               comparisons = list(c("rhesus", "chimp"), 
#                                  c("rhesus", "human"), 
#                                  c("chimp", "human"))) %>%
#   add_xy_position(x = "region", dodge = 0.8)
# 
# 
# 
# region_max <- axon_data %>%
#   group_by(region) %>%
#   summarise(max_value = max(value, na.rm = TRUE))
# 
# # 将最大值合并到统计结果中
# stat_results <- stat_results %>%
#   left_join(region_max, by = "region")
# 
# 
# # 创建自定义的y位置，增加垂直间距
# stat_results <- stat_results %>%
#   mutate(
#     y.position = case_when(
#       group1 == "rhesus" & group2 == "chimp" ~ max_value * 1.7,     # 最低位置
#       group1 == "rhesus" & group2 == "human" ~ max_value * 1.4,    # 中间位置
#       group1 == "chimp" & group2 == "human" ~ max_value * 1.1,    # 最高位置
#       TRUE ~ max_value * 1.15
#     )
#   )
# 
# # 使用坐标轴扩展来确保标注不被裁剪
# y_range <- range(axon_data$value, na.rm = TRUE)
# y_expand <- c(0, (max(stat_results$y.position, na.rm = TRUE) - y_range[2]) * 1.1)
# 
# p2 <- ggplot(axon_data, aes(x = region, y = value, fill = species)) +
#   geom_violin() +
#   stat_pvalue_manual(
#     stat_results,
#     label = "p.adj.signif",
#     tip.length = 0.01,
#     bracket.size = 0.3,
#     size = 3.5,  # 稍微减小字体大小以适应更高的位置
#     bracket.nudge.y = 0.1 # 稍微调整支架的垂直位置
#   ) +
#   scale_fill_manual(values = c("human" = "#5a8fbf", 
#                                "chimp" = "#60a764", 
#                                "rhesus" = "#813117")) +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +  # 增加y轴上方的扩展空间
#   labs(x = "Brain Region", y = "Expression", 
#        fill = "Species",
#        title = "") +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#     axis.text.y = element_text(size = 10),
#     axis.title = element_text(size = 12),
#     plot.title = element_text(size = 14, hjust = 0.5),
#     legend.position = "right"
#   )
# 
# print(p2)
# 
# setwd("F://03-11")
# pdf("vinplot_region_Primates.pdf",width =12,height =2)
# #pdf(paste("pheatmap_sinh","Primates","pdf",sep = "."),width =4 ,height = 4)
# print(p2)
# dev.off()

# library(tidyverse)
# library(rstatix)
# library(ggpubr)
# 
# # 设置物种顺序
# axon_data <- axon_data %>%
#   mutate(
#     species = factor(
#       species,
#       levels = c("rhesus", "chimp", "human")
#     )
#   )
# 
# #==========================================================
# # 1. 将同一region、同一gene在三个物种中的表达值匹配起来
# #==========================================================
# 
# axon_wide <- axon_data %>%
#   select(region, symbol, species, value) %>%
#   # 如果同一region-gene-species存在多行，需要先汇总
#   group_by(region, symbol, species) %>%
#   summarise(
#     value = mean(value, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   pivot_wider(
#     names_from = species,
#     values_from = value
#   ) %>%
#   # 仅保留三个物种都有表达值的基因
#   filter(
#     !is.na(rhesus),
#     !is.na(chimp),
#     !is.na(human)
#   )
# 
# # 查看每个region实际参与配对检验的基因数量
# axon_wide %>%
#   count(region, name = "n_paired_genes") %>%
#   print()
# 
# 
# #==========================================================
# # 2. 在每个region中，以gene为配对单位进行配对t检验
# #==========================================================
# 
# comparisons <- tribble(
#   ~group1,   ~group2,
#   "rhesus",  "chimp",
#   "rhesus",  "human",
#   "chimp",   "human"
# )
# 
# stat_results <- axon_wide %>%
#   group_by(region) %>%
#   group_modify(~{
#     
#     region_data <- .x
#     
#     comparisons %>%
#       rowwise() %>%
#       mutate(
#         test_result = list(
#           t.test(
#             x = region_data[[group1]],
#             y = region_data[[group2]],
#             paired = TRUE
#           )
#         ),
#         statistic = unname(test_result$statistic),
#         df = unname(test_result$parameter),
#         p = test_result$p.value,
#         n = sum(
#           complete.cases(
#             region_data[[group1]],
#             region_data[[group2]]
#           )
#         )
#       ) %>%
#       select(-test_result) %>%
#       ungroup()
#   }) %>%
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(
#     # 每个region内对3次比较进行BH校正
#     p.adj = p.adjust(p, method = "BH"),
#     p.adj.signif = case_when(
#       p.adj <= 0.0001 ~ "****",
#       p.adj <= 0.001  ~ "***",
#       p.adj <= 0.01   ~ "**",
#       p.adj <= 0.05   ~ "*",
#       TRUE            ~ "ns"
#     )
#   ) %>%
#   ungroup()
# 
# 
# #==========================================================
# # 3. 设置显著性标注位置
# #==========================================================
# 
# region_range <- axon_data %>%
#   group_by(region) %>%
#   summarise(
#     min_value = min(value, na.rm = TRUE),
#     max_value = max(value, na.rm = TRUE),
#     value_range = max_value - min_value,
#     .groups = "drop"
#   ) %>%
#   mutate(
#     # 避免某个region所有值相同时value_range为0
#     value_range = if_else(
#       value_range == 0,
#       abs(max_value) * 0.1 + 1,
#       value_range
#     )
#   )
# 
# stat_results <- stat_results %>%
#   left_join(region_range, by = "region") %>%
#   mutate(
#     y.position = case_when(
#       group1 == "rhesus" & group2 == "chimp" ~
#         max_value + value_range * 0.10,
#       
#       group1 == "rhesus" & group2 == "human" ~
#         max_value + value_range * 0.25,
#       
#       group1 == "chimp" & group2 == "human" ~
#         max_value + value_range * 0.40
#     )
#   )
# 
# 
# #==========================================================
# # 4. 绘图
# #==========================================================
# 
# p2 <- ggplot(
#   axon_data,
#   aes(
#     x = region,
#     y = value,
#     fill = species
#   )
# ) +
#   geom_violin(
#     trim = FALSE,
#     position = position_dodge(width = 0.8)
#   ) +
#   stat_pvalue_manual(
#     stat_results,
#     label = "p.adj.signif",
#     xmin = "group1",
#     xmax = "group2",
#     y.position = "y.position",
#     tip.length = 0.01,
#     bracket.size = 0.3,
#     size = 3.5
#   ) +
#   scale_fill_manual(
#     values = c(
#       "human" = "#5a8fbf",
#       "chimp" = "#60a764",
#       "rhesus" = "#813117"
#     )
#   ) +
#   scale_y_continuous(
#     expand = expansion(mult = c(0.05, 0.30))
#   ) +
#   labs(
#     x = "Brain Region",
#     y = "Expression",
#     fill = "Species",
#     title = ""
#   ) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(
#       angle = 45,
#       hjust = 1,
#       size = 10
#     ),
#     axis.text.y = element_text(size = 10),
#     axis.title = element_text(size = 12),
#     plot.title = element_text(size = 14, hjust = 0.5),
#     legend.position = "right"
#   )
# 
# print(p2)
library(tidyverse)
library(rstatix)
library(ggpubr)

species_order <- c("rhesus", "chimp", "human")

# 假设基因名列为gene
# 同一个region中，同一个gene在三个物种间进行配对
axon_data_paired <- axon_data %>%
  mutate(
    species = factor(species, levels = species_order)
  ) %>%
  # 避免同一个region-gene-species存在重复行
  group_by(region, symbol, species) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 只保留三个物种中都有数值的基因
  group_by(region, symbol) %>%
  filter(
    n_distinct(species) == 3,
    all(species_order %in% as.character(species))
  ) %>%
  ungroup() %>%
  # paired t-test实际根据顺序配对，因此必须按gene排序
  arrange(region, species, symbol)


write.table(axon_data_paired,"axon_data_paired_Homo.txt",sep = "\t",quote = F)
table(axon_data_paired$region,axon_data_paired$species)
# 每个region中，以gene为配对单位进行配对t检验
stat_results <- axon_data_paired %>%
  group_by(region) %>%
  t_test(
    value ~ species,
    paired = TRUE, #这里控制配对t检验
    #pair_by = symbol, #显式配对
    comparisons = list(
      c("rhesus", "chimp"),
      c("rhesus", "human"),
      c("chimp", "human")
    ),
    p.adjust.method = "BH"
  ) %>%
  # 生成每个region内部三个物种的横坐标位置
  add_xy_position(
    x = "region",
    dodge = 0.8
  )

#############################################

gene_wide <- axon_data_paired %>%
  select(symbol, region, species, value) %>%
  pivot_wider(
    names_from = species,
    values_from = value
  )

gene_test <- gene_wide %>%
  group_by(symbol) %>%
  summarise(
    n_region = sum(complete.cases(rhesus, chimp, human)),
    
    rhesus_vs_chimp_p = 
      ifelse(
        sum(complete.cases(rhesus, chimp)) >= 2,
        t.test(chimp, rhesus, paired = TRUE)$p.value,
        NA_real_
      ),
    
    rhesus_vs_human_p =
      ifelse(
        sum(complete.cases(rhesus, human)) >= 2,
        t.test(human, rhesus, paired = TRUE)$p.value,
        NA_real_
      ),
    
    mean_rhesus = mean(rhesus, na.rm = TRUE),
    mean_chimp  = mean(chimp, na.rm = TRUE),
    mean_human  = mean(human, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    rhesus_vs_chimp_padj =
      p.adjust(rhesus_vs_chimp_p, method = "BH"),
    
    rhesus_vs_human_padj =
      p.adjust(rhesus_vs_human_p, method = "BH")
  )

rhesus_to_human_down <- gene_test %>%
  filter(
    rhesus_vs_human_padj < 0.05,
    mean_human < mean_rhesus
  )

rhesus_to_chimp_down <- gene_test %>%
  filter(
    rhesus_vs_chimp_padj < 0.05,
    mean_chimp < mean_rhesus
  )

down_both <- gene_test %>%
  filter(
    rhesus_vs_chimp_padj < 0.05,
    rhesus_vs_human_padj < 0.05,
    mean_chimp < mean_rhesus,
    mean_human < mean_rhesus
  )

write.table(down_both,"down_both_chimp_human.txt",sep = "\t",quote = F)


########################################

# 计算每个region的数据范围
region_range <- axon_data_paired %>%
  group_by(region) %>%
  summarise(
    max_value = max(value, na.rm = TRUE),
    min_value = min(value, na.rm = TRUE),
    value_range = max_value - min_value,
    .groups = "drop"
  ) %>%
  mutate(
    value_range = if_else(
      value_range == 0,
      pmax(abs(max_value) * 0.1, 0.1),
      value_range
    )
  )


# 只覆盖y.position
# add_xy_position生成的xmin和xmax继续保留
stat_results <- stat_results %>%
  left_join(region_range, by = "region") %>%
  mutate(
    comparison_order = case_when(
      group1 == "rhesus" & group2 == "chimp" ~ 1,
      group1 == "rhesus" & group2 == "human" ~ 2,
      group1 == "chimp"  & group2 == "human" ~ 3
    ),
    y.position = max_value + value_range * 0.3 * comparison_order
  )


# p2 <- ggplot(
#   axon_data_paired,
#   aes(
#     x = region,
#     y = value,
#     fill = species
#   )
# ) +
#   geom_violin(
#     trim = FALSE,
#     position = position_dodge(width = 0.8)
#   ) +
#   stat_pvalue_manual(
#     stat_results,
#     label = "p.adj.signif",
#     
#     # 不要写：
#     # xmin = "group1",
#     # xmax = "group2"
#     
#     tip.length = 0.01,
#     bracket.size = 0.3,
#     size = 3.5
#   ) +
#   scale_fill_manual(
#     values = c(
#       "human" = "#5a8fbf",
#       "chimp" = "#60a764",
#       "rhesus" = "#813117"
#     )
#   ) +
#   scale_y_continuous(
#     expand = expansion(mult = c(0.05, 0.30))
#   ) +
#   labs(
#     x = "Brain Region",
#     y = "Expression",
#     fill = "Species",
#     title = ""
#   ) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(
#       angle = 45,
#       hjust = 1,
#       size = 10
#     ),
#     axis.text.y = element_text(size = 10),
#     axis.title = element_text(size = 12),
#     legend.position = "right"
#   )
# 
# print(p2)


p3 <- ggplot(
  axon_data_paired,
  aes(
    x = region,
    y = value,
    fill = species
  )
) +
  geom_violin(
    trim = FALSE,
    position = position_dodge(width = 0.8)
  ) +
  
  # 添加中位数横线
  stat_summary(
    aes(group = species),
    fun = median,
    # geom = "crossbar",
    # width = 0.18,
    # fatten = 0,
    # color = "black",
    # linewidth = 0.5,
    # position = position_dodge(width = 0.8)
    geom = "point",
    shape = 21,
    size = 2.5,
    fill = "white",
    color = "black",
    stroke = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  
  stat_pvalue_manual(
    stat_results,
    label = "p.adj.signif",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      "human" = "#5a8fbf",
      "chimp" = "#60a764",
      "rhesus" = "#813117"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.30))
  ) +
  labs(
    x = "Brain Region",
    y = "RPKM",
    fill = "Species",
    title = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

print(p3)
pdf("vinplot_non-Primates.pdf",width =12,height =2)
#pdf(paste("pheatmap_sinh","Primates","pdf",sep = "."),width =4 ,height = 4)
print(p3)
dev.off()

#人科 153
#狭鼻类 102
#类人猿下目 124
#灵长目35
#非灵长目65

```
