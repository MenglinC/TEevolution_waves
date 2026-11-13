#setwd("F:/tmp//03-17")

library(Matrix)
#install_github("linxihui/NNLM")
library(NNLM)
library(ggplot2)

#############################################################################
setwd("F://10-08-画图")
#setwd("F://03-17")
data<-read.csv("F:/tmp/03-17-nmf2/merged_intervals.csv",row.names = 1)
#data$Quies <- data$Type_15_Count < 0.5 #这一步多此一举了
# tmp<-cbind(Interval=row.names(data),data)
annotation<-read.csv("F:/tmp/03-17-nmf2/merge_intervals_group.csv",row.names = 1)
# QueisDat<-merge(tmp,annotation,by="Interval",all.x = T)
# library(tidyr)
# QueisDat <- QueisDat %>%
#   unite("group", lineage, family, sep = "_", remove = FALSE)
# 
# 
# # 计算每个 group 中 Quies 为 TRUE 和 FALSE 的比例
# prop_data <- QueisDat %>%
#   group_by(group) %>%
#   summarise(
#     `TRUE`  = mean(!Quies, na.rm = TRUE), #换一下
#     `FALSE` = mean(Quies, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   pivot_longer(cols = c("TRUE", "FALSE"), 
#                names_to = "Quies_Status", 
#                values_to = "Proportion")
# 
# 
# prop_data$group<-factor(prop_data$group,
#                         levels = c("non-Primates_L1","Primates_L1","Simiiformes_L1","Catarrhini_L1","Hominoidea_L1",
#                                    "non-Primates_ERV1",
#                                    "Primates_ERV1",
#                                    "Simiiformes_ERV1",
#                                    "Catarrhini_ERV1",
#                                    "Hominoidea_ERV1",
#                                    "Primates_Alu",
#                                    "Simiiformes_Alu",
#                                    "Catarrhini_Alu",
#                                    "Hominoidea_Alu",
#                                    "Simiiformes_ERVK",
#                                    "Catarrhini_ERVK",
#                                    "Hominoidea_ERVK",
#                                    "Hominoidea_SVA"
#                                    ))
# 
# 
# # 绘制堆叠条形图
# pq<-ggplot(prop_data, aes(x = group, y = Proportion, fill = Quies_Status)) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(
#     title = "Proportion of Quies = TRUE vs FALSE by Group",
#     x = "Group",
#     y = "Proportion",
#     fill = "Quies Status"
#   ) +
#   
#   scale_fill_manual(values = c("TRUE" = "#909495", "FALSE" = "#c31f26")) +
#   theme_classic() +
#   geom_vline(xintercept = 5.5, linetype = "dashed", color = "black", size = 0.5)+
#   geom_vline(xintercept = 10.5, linetype = "dashed", color = "black", size = 0.5)+
#   geom_vline(xintercept = 14.5, linetype = "dashed", color = "black", size = 0.5)+
#   geom_vline(xintercept = 17.5, linetype = "dashed", color = "black", size = 0.5)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# pq
# pdf("Quies.pdf",width =8,height =3)
# pq
# dev.off()
# 
# 
# ###################################################################################
# #在之前，我们分配了signiture，我们需要把每个group中属于该signiture的标注为ture，然后画出来这样的比例图
# site_assignment2<-read.csv("F:/tmp/03-17-nmf2/site_assignment_0814.csv") #只有734651行；也就是我们只考虑活跃的部分
# site_assignment2<-read.csv("F:/tmp/03-20/site_assignment_0814.csv") #只有734651行；也就是我们只考虑活跃的部分
# 
# colnames(site_assignment2)[1]<-"Interval"
# data_s<-merge(QueisDat,site_assignment2)
# head(data_s)
# 
# data_s$judgement <- data_s$x == 1 #这个可以反复修改
# head(data_s)
# 
# 
# # 计算每个 group 中 Quies 为 TRUE 和 FALSE 的比例
# prop_data <- data_s %>%
#   group_by(group) %>%
#   summarise(
#     `TRUE`  = mean(judgement, na.rm = TRUE), #换一下
#     `FALSE` = mean(!judgement, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   pivot_longer(cols = c("TRUE", "FALSE"), 
#                names_to = "signiture5", 
#                values_to = "Proportion")
# 
# 
# prop_data$group<-factor(prop_data$group,
#                         levels = c("non-Primates_L1","Primates_L1","Simiiformes_L1","Catarrhini_L1","Hominoidea_L1",
#                                    "non-Primates_ERV1",
#                                    "Primates_ERV1",
#                                    "Simiiformes_ERV1",
#                                    "Catarrhini_ERV1",
#                                    "Hominoidea_ERV1",
#                                    "Primates_Alu",
#                                    "Simiiformes_Alu",
#                                    "Catarrhini_Alu",
#                                    "Hominoidea_Alu",
#                                    "Simiiformes_ERVK",
#                                    "Catarrhini_ERVK",
#                                    "Hominoidea_ERVK",
#                                    "Hominoidea_SVA"
#                         ))
# 
# 
# # # 绘制堆叠条形图
# # pq<-ggplot(prop_data, aes(x = group, y = Proportion, fill = signiture5)) +
# #   geom_bar(stat = "identity", position = "stack") +
# #   labs(
# #     title = "Proportion of Quies = TRUE vs FALSE by Group",
# #     x = "Group",
# #     y = "Proportion",
# #     fill = "Quies Status"
# #   ) +
# #   
# #   scale_fill_manual(values = c("TRUE" = "#faba7d", "FALSE" = "grey")) + ##a39dcd 
# #   theme_classic() +
# #   geom_vline(xintercept = 5.5, linetype = "dashed", color = "black", size = 0.5)+
# #   geom_vline(xintercept = 10.5, linetype = "dashed", color = "black", size = 0.5)+
# #   geom_vline(xintercept = 14.5, linetype = "dashed", color = "black", size = 0.5)+
# #   geom_vline(xintercept = 17.5, linetype = "dashed", color = "black", size = 0.5)+
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # 
# # pq
# # 
# # pdf("signiture1.pdf",width =8,height =3)
# # pq
# # dev.off()
# 
# ########################################################################################
# #data_s #x列画饼图
# # data_s$group<-factor(data_s$group,
# #                         levels = c("non-Primates_L1","Primates_L1","Simiiformes_L1","Catarrhini_L1","Hominoidea_L1",
# #                                    "non-Primates_ERV1",
# #                                    "Primates_ERV1",
# #                                    "Simiiformes_ERV1",
# #                                    "Catarrhini_ERV1",
# #                                    "Hominoidea_ERV1",
# #                                    "Primates_Alu",
# #                                    "Simiiformes_Alu",
# #                                    "Catarrhini_Alu",
# #                                    "Hominoidea_Alu",
# #                                    "Simiiformes_ERVK",
# #                                    "Catarrhini_ERVK",
# #                                    "Hominoidea_ERVK",
# #                                    "Hominoidea_SVA"
# #                         ))
# # 
# # 
# # data_s$x<-as.character(data_s$x)
# # #方法1：使用 ggplot2 分面饼图
# # library(dplyr)
# # pie_plot <- data_s %>%
# #   count(group, x) %>%                     # 计算每个组内各类别的频数
# #   group_by(group) %>%
# #   mutate(prop = n / sum(n)) %>%           # 计算比例
# #   ggplot(aes(x = "", y = prop, fill = x)) +
# #   geom_bar(stat = "identity", width = 1) + # 条形图
# #   coord_polar("y", start = 0) +           # 转换为极坐标（饼图）
# #   facet_wrap(~ group) +                   # 按分组变量分面
# #   scale_fill_manual(values = c("1"="#faba7d",
# #                     "2"="#77ccce",
# #                     "3"="#d1d3d4",
# #                     "4"="#ffeba6",
# #                     "5"="#a39dcd",
# #                     "6"="#dae2bf")) +   # 自定义颜色
# #   theme_void() +                          # 去掉背景坐标轴
# #   labs(title = "") +
# #   theme(legend.title = element_blank())   # 去掉图例标题
# # 
# # # 显示图形
# # print(pie_plot)
# # 
# # pdf("pie_0814.pdf",width =12,height =8)
# # pie_plot
# # dev.off()
# 
# 
# 
# 
# 




#####################################################################################
mat<-data[data$Type_15_Count<0.5,] #这里做了初步的筛选
mat<- as.matrix(mat)
storage.mode(mat) <- "numeric" 
non_zero_rows <- rowSums(mat) > 0
mat <- mat[non_zero_rows, ]
non_zero_cols <- colSums(mat) > 0
mat <- mat[, non_zero_cols]
cat("处理后矩阵维度:", dim(mat), "\n")
sparse_data <- as(as.matrix(mat), "dgCMatrix")
#rm(data)
#rm(mat)




##過濾基本都爲quies的情況


set.seed(1226)
sample_data <- as.matrix(sparse_data) * 100
#k_values <- c(2, 3, 4, 5, 6, 7, 8, 10, 12)
#results <- data.frame(k = integer(), kl_loss = numeric())

# cat("開始評估 k 值...\n")
# 
# for (k in k_values) {
#   cat("正在計算 k =", k, "...\n")
#   
#   # 執行 NNLM
#   # 我們需要確保 max.iter 夠大以保證收斂
#   fit <- nnmf(sample_data, k = k, loss = "mkl", max.iter = 100, rel.tol = 1e-7)
#   
#   # 安全地提取最後一個 loss 值
#   # tail(fit$dist, 1) 會取得最後一個記錄的目標函數值
#   current_loss <- tail(as.numeric(fit$dist), 1)
#   
#   if(length(current_loss) == 0) {
#     # 如果 dist 是空的（極少見），手動計算 KL
#     # 這裡提供一個備選計算方式
#     W <- fit$W
#     H <- fit$H
#     V_hat <- W %*% H
#     V_hat[V_hat <= 0] <- 1e-10 # 避免 log(0)
#     current_loss <- sum(sample_data * log(sample_data / V_hat + 1e-10) - sample_data + V_hat)
#   }
#   
#   results <- rbind(results, data.frame(k = k, kl_loss = current_loss))
# }

# # 繪製 KL 散度下降圖
# pt<-ggplot(results, aes(x = k, y = kl_loss)) +
#   geom_line(color = "black", linewidth = 1) +
#   geom_point(color = "red", size = 3) +
#   scale_x_continuous(breaks = k_values) +
#   labs(title = "NMF (MKL Loss) Rank Selection",
#        subtitle = "Finding the Elbow Point",
#        x = "Number of Signatures (k)",
#        y = "Final MKL Objective Value") +
#   theme_classic()
# pt
# pdf("NMF_rank_selection.pdf",width =3,height =3)
# pt
# dev.off()





# ##看圖的話，k=6比較合適
# # 測試第一次
# set.seed(1)
# fit1 <- nnmf(sample_data, k = 6, loss = "mkl")
# 
# # 測試第二次
# set.seed(2)
# fit2 <- nnmf(sample_data, k = 6, loss = "mkl")
# 
# # 計算兩次結果 H 矩陣的相關性
# # 由於 Signature 的順序在不同次運算中可能會重排，
# # 我們觀察相關性矩陣，看是否每一列都能在另一組中找到對應的高相關項
# cor_matrix <- cor(t(fit1$H), t(fit2$H))
# print(round(cor_matrix, 3))
# library(pheatmap)
# pl<-pheatmap(cor_matrix,cluster_rows = FALSE, 
#              cluster_cols = FALSE)
# pdf("NMF_cor.pdf",width =3,height =2.5)
# pl
# dev.off()


#########################################################
#fit_nnlm <- nnmf(sample_data, k = 6, loss = "mkl", max.iter = 500, rel.tol = 1e-7)
#Target tolerance not reached. Try a larger max.iter.
set.seed(1226) #加一个seed
fit_nnlm <- nnmf(sample_data, k = 6, loss = "mkl", max.iter = 2000, rel.tol = 1e-7)
##这一步耗时良久
##设置了1000，似乎还没有到，这个暂时就这样了；
#2000的话虽然时间长一点，但是是可以了；
#write.csv(fit_nnlm,"fit_nnlm.2000.csv")
# 3. 檢查重建情況
W <- fit_nnlm$W
H <- fit_nnlm$H
sig_scale <- rowSums(H)

H_norm <- H / sig_scale #对H按照行进行归一化了
W_adj <- sweep(W, 2, sig_scale, "*")

#cat("W max:", max(W), " H max:", max(H), "\n")
# 計算 R2 (注意 RSS 計算時要對齊量級)
#rss <- sum((sample_data - W %*% H)^2)
#tss <- sum(sample_data^2)
#cat("NNLM Sample R2:", 1 - (rss/tss), "\n")

# 給 Signature 命名,是根据H图的行归一化后的结果再进行的命名
colnames(W_adj) <- c("Strong.transcription.like","Weak.transcription.like","Quiescent.low.like","PolyComb.like","Heterochromatin.like","cCRE.like")
rownames(H_norm) <- c("Strong.transcription.like","Weak.transcription.like","Quiescent.low.like","PolyComb.like","Heterochromatin.like","cCRE.like")
colnames(H_norm) <- c("Active.TSS","Flanking.Active.TSS","Transcr.at.gene.5.and.3","Strong.transcription","Weak.transcription",
                      "Genic.enhancers","Enhancers","ZNF.genes.and.repeats","Heterochromatin","Bivalent.Poised.TSS",
                      "Flanking.Bivalent.TSS.Enh","Bivalent.Enhancer","Repressed.PolyComb","Weak.Repressed.PolyComb","Quiescent.Low") #sparse_data的作用就是来命名的吗？
library(pheatmap)
pr<-pheatmap(H_norm, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, main = "Chromatin State Signatures (H matrix)",
         display_numbers = TRUE, color = colorRampPalette(c("white", "red"))(100))
#pdf("NMF_0322.pdf",width =8,height =4)
pr
#dev.off()

ggsave(
  "NMF-H_norm-matrix.tiff",
  plot = pr,
  width = 10,
  height = 4,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)






site_assignment <- apply(W_adj, 1, which.max)
#write.csv(site_assignment,"site_assignment_0317.csv")
write.csv(site_assignment,"site_assignment_0814.csv") #要及时保存这个，做的真的很崩溃很多个大文件
write.csv(site_assignment,"site_assignment_1008.csv") #要及时保存这个，做的真的很崩溃很多个大文件，做到这一步的时候保存，然后重启
####################################
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






##############################################
##260814，从这里开始
#site_assignment<-read.csv("site_assignment_0317.csv")

# 3. 分層排序與採樣 (在 1000 行內操作)
n_show <- 500 # 每個 Signature 展示前 200 行，總共 600 行左右
selected_indices <- list()

for(i in 1:6) {
  idx <- which(site_assignment == i)
  if(length(idx) > 0) {
    # 根據 W 權重排序
    sig_scores <- W[idx, i]
    ordered_inner_idx <- idx[order(sig_scores, decreasing = TRUE)]
    selected_indices[[i]] <- head(ordered_inner_idx, n_show)
  }
}
#问题就出在这里


final_plot_indices <- unlist(selected_indices)

# 4. 提取對應的原始數據（比例值）
# 注意：這裡直接從 V_sample_for_plot 提取，它已經是 1000 行的矩陣了
plot_mat <- sparse_data[final_plot_indices, ]

# 5. 準備側面標註
anno_row <- data.frame(
  Signature = factor(rep(paste0("Signature_", 1:6), 
                         times = sapply(selected_indices, length)))
)
rownames(anno_row) <- paste0("Row", 1:nrow(plot_mat))
rownames(plot_mat) <- rownames(anno_row)

ann_colors<-list(Signature=c(Signature_1="#faba7d",
                             Signature_2="#77ccce",
                             Signature_3="#d1d3d4",
                             Signature_4="#ffeba6",
                             Signature_5="#a39dcd",
                             Signature_6="#dae2bf"))

# 將 PCA 結果與標籤合併
# plot_data <- data.frame(
#  PC1 = pca_res$x[, 1],
#  PC2 = pca_res$x[, 2]
#)

#plot_df2<-cbind(Interval=row.names(plot_df),plot_df)
#annotation<-read.csv("merge_intervals_group.csv",row.names = 1)

#plot_df3<-merge(plot_df2,annotation,by="Interval",all.x = T)






# 6. 繪圖
library(pheatmap)
p_raw <- pheatmap(t(plot_mat), 
                  cluster_rows = FALSE, 
                  cluster_cols = FALSE, 
                  show_rownames = TRUE,
                  show_colnames = FALSE,
                  annotation_col = anno_row,
                  main = "",
                  annotation_colors = ann_colors,
                  #color = colorRampPalette(c("#f4f4f4", "#d29193", "#981a32"))(100))
                  color = colorRampPalette(c("#0a0c0d", "#f5eb0d", "#eb2627"))(50) #这个沿用了图1的结果
)



p_raw



setwd("F://tmp//03-20")
pdf("NMF_raw_0814.pdf",width =5,height =2)
p_raw
dev.off()








