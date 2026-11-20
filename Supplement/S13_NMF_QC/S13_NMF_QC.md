# S13

## S13.1

```r
data<-read.csv("F:/tmp/03-17-nmf2/merged_intervals.csv",row.names = 1)
mat<-data[data$Type_15_Count<0.5,] #这里做了初步的筛选,因为绝大多数是Quies，会干扰对于特征的提取
mat<- as.matrix(mat)
storage.mode(mat) <- "numeric" 
non_zero_rows <- rowSums(mat) > 0
mat <- mat[non_zero_rows, ]
non_zero_cols <- colSums(mat) > 0
mat <- mat[, non_zero_cols]
cat("处理后矩阵维度:", dim(mat), "\n")
sparse_data <- as(as.matrix(mat), "dgCMatrix")
rm(data)
rm(mat)

set.seed(1226)
sample_data <- as.matrix(sparse_data) * 100
k_values <- c(2, 3, 4, 5, 6, 7, 8, 10, 12)
results <- data.frame(k = integer(), kl_loss = numeric())

cat("開始評估 k 值...\n")

for (k in k_values) {
  cat("正在計算 k =", k, "...\n")
  
  # 執行 NNLM
  # 我們需要確保 max.iter 夠大以保證收斂
  fit <- nnmf(sample_data, k = k, loss = "mkl", max.iter = 100, rel.tol = 1e-7)
  
  # 安全地提取最後一個 loss 值
  # tail(fit$dist, 1) 會取得最後一個記錄的目標函數值
  current_loss <- tail(as.numeric(fit$dist), 1)
  
  if(length(current_loss) == 0) {
    # 如果 dist 是空的（極少見），手動計算 KL
    # 這裡提供一個備選計算方式
    W <- fit$W
    H <- fit$H
    V_hat <- W %*% H
    V_hat[V_hat <= 0] <- 1e-10 # 避免 log(0)
    current_loss <- sum(sample_data * log(sample_data / V_hat + 1e-10) - sample_data + V_hat)
  }
  
  results <- rbind(results, data.frame(k = k, kl_loss = current_loss))
}

# 繪製 KL 散度下降圖
pt<-ggplot(results, aes(x = k, y = kl_loss)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(color = "red", size = 3) +
  scale_x_continuous(breaks = k_values) +
  labs(title = "NMF (MKL Loss) Rank Selection",
       subtitle = "Finding the Elbow Point",
       x = "Number of Signatures (k)",
       y = "Final MKL Objective Value") +
  theme_classic()
pt
pdf("NMF_rank_selection.pdf",width =3,height =3)
pt
dev.off()
```

## S13.2 这个可选我觉得

```r
##看圖的話，k=6比較合適
# 測試第一次
set.seed(1)
fit1 <- nnmf(sample_data, k = 6, loss = "mkl")

# 測試第二次
set.seed(2)
fit2 <- nnmf(sample_data, k = 6, loss = "mkl")

# 計算兩次結果 H 矩陣的相關性
# 由於 Signature 的順序在不同次運算中可能會重排，
# 我們觀察相關性矩陣，看是否每一列都能在另一組中找到對應的高相關項
cor_matrix <- cor(t(fit1$H), t(fit2$H))
print(round(cor_matrix, 3))
library(pheatmap)
pl<-pheatmap(cor_matrix,cluster_rows = FALSE, 
             cluster_cols = FALSE)
pdf("NMF_cor.pdf",width =3,height =2.5)
pl
dev.off()

```


