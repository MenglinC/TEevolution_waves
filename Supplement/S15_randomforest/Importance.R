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
################################################
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
    # 若年龄是连续变量，则保持数值型
  )

####################################
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


# 选择效应量最大的两个因素（根据eta_sq结果）
p1 <- ggplot(df, aes(x = family, y = combo_score, fill = clade)) +
  geom_boxplot() +
  facet_grid(. ~ states) +  # 添加组织维度
  scale_fill_brewer(palette = "Set1") +
  labs(title = "TE类型与染色质修饰的交互效应",
       x = "TE类型", 
       y = "Giggle Score") +
  theme_bw()

p1
##############################################
df$tissue2<-factor(df$tissue2,c("iPSC","ESC","ES-deriv","Mesench",
                                "Blood & T-cell","HSC & B-cell",
                                "Epithelial","Sm Muscle","Muscle","Lung","Heart","Digestive","Thymus",
                                "Brain","Neurosph",
                                "Cancer Cell Line","Other"))
df$clade<-factor(df$clade,levels=c("non-Primates","Primates","Simiiformes","Catarrhini","Hominoidea"))

dfsub<-df[df$family=="Alu"& df$clade=="Primates"&(df$states=="Tx"),]
dfsub<-df[df$family=="SVA"& df$clade=="Hominoidea"&(df$states=="TxWk"),]
dfsub<-df[df$family=="ERV1"&(df$states=="Enh"),]
dfsub<-df[(df$family=="ERVK"|df$family=="ERV1")&(df$states=="TssA"|df$states=="TssFlink")&(df$tissue2=="ESC"|df$tissue2=="iPSC")&df$clade%in%c("Simiiformes","Catarrhini","Hominoidea"),]
dfsub<-df[(df$family=="ERVK"|df$family=="ERV1")&(df$states=="Enh")&(df$tissue2=="ESC"|df$tissue2=="iPSC")&df$clade%in%c("Simiiformes","Catarrhini","Hominoidea"),]

p2<-ggplot(dfsub, aes(x = clade, y = combo_score, fill = family)) +
  #geom_boxplot(fill="#549357") + #Alu
  geom_boxplot() +
  #facet_grid(. ~ tissue2) +  # 添加组织维度
  #scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values=c("#8cafd6","#fada69")) +
  labs(title = "",
       x = "", 
       y = "") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
p2
pdf(paste("ERV.Enh.clade.pdf"),width =2.5,height = 3)
#pdf(paste("pheatmap_sinh","Primates","pdf",sep = "."),width =4 ,height = 4)
print(p2)
dev.off()



pdf(paste("Primate.Alu.Tx.pdf"),width =4 ,height = 2)
#pdf(paste("pheatmap_sinh","Primates","pdf",sep = "."),width =4 ,height = 4)
print(p2)
dev.off()





