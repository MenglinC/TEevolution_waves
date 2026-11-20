# S15 random forest

```r
#file:\\192.168.79.85\Shares\Share\Manuscripts\张秀秀\数据备份\实验记录\历年工作路径\2025\06-18-组蛋白修饰\Importance.R

giggleDat<-read.table("all.giggle.result",sep = "\t",header = T)[,c(1,2,9)]
library(tidyr)
library(dplyr)
giggleDat2<-giggleDat %>%
  separate(
    col = file,             # 要分割的列
    into = c("tissue", "states"),  # 新列名
    sep = "_(?=[^_]+$)",    # 正则表达式：匹配最后一个下划线
    remove = FALSE,         # 保留原始列
    fill = "right",         # 无下划线时右侧填充NA
    convert = FALSE         # 不自动转换类型
  )
############
name1<-read.table("EDACC_NAME.txt")
name2<-read.table("new_groups.txt",sep = "\t")
name<-merge(name1,name2,by="V1")
############
label<-names(table(name$V2.y))
for (i in 1:length(label)) {
  giggleDat2[giggleDat2$tissue%in%name[name$V2.y==label[i],]$V2.x,"tissue2"]<-label[i]
}
giggleDat3<-giggleDat2 %>%
  separate(
    col = group,             # 要分割的列
    into = c("clade", "family"),  # 新列名
    sep = "_",    # 正则表达式：匹配最后一个下划线
    remove = FALSE,         # 保留原始列
    fill = "right",         # 无下划线时右侧填充NA
    convert = FALSE         # 不自动转换类型
  )
giggleDat4<-giggleDat3[,c(2,3,6,8,7)]

library(tidyverse)   # 数据处理与可视化
library(rstatix)     # 统计检验
library(car)         # 方差分析
library(ggpubr)      # 高级绘图
BiocManager::install("randomForest")
library(randomForest) # 随机森林评估变量重要性
df<-giggleDat4
df <- df %>%
  mutate(
    family = as.factor(family),
    clade = as.factor(clade),
    states = as.factor(states),
    tissue2 = as.factor(tissue2)  # 如果年龄是分类变量

set.seed(123)
rf_model <- randomForest(
  combo_score ~ family + states + clade + tissue2,
  data = df,
  importance = TRUE,
  ntree = 500
)

# 提取变量重要性（%IncMSE表示精度下降程度）
var_imp <- importance(rf_model, type = 1, scale = TRUE)
var_imp_df <- data.frame(
  Variable = rownames(var_imp),
  Importance = var_imp[, "%IncMSE"]
) %>%
  arrange(desc(Importance))

# 可视化重要性
p1<-ggplot(var_imp_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col(fill = "#415094", width = 0.6) +
  coord_flip() +
  labs(title = "", 
       x = "vector", 
       y = "Importance(%IncMSE)") +
  theme_classic()

pdf(paste("Importance.pdf"),width =2 ,height = 2)
#pdf(paste("pheatmap_sinh","Primates","pdf",sep = "."),width =4 ,height = 4)
print(p1)
dev.off()

```
