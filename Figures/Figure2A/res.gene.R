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
