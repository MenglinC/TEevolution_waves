# S4

## S4.1 consistence
```shell


```

## S4.2 cent-enrichment

### S4.2.1 giggle
```shell
#workplace:/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/Cent_enrichment/

ls *.bed | gargs -p 30 "bgzip {}"

/home/xxzhang/workplace/software/giggle/scripts/sort_bed "./CenSat_by_family/*.bed.gz" CenSat_sort/ 4

giggle index -i "CenSat_sort/*gz" -o CenSat_sort_b -f -s

TE_DIR="/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/clade_family"
INDEX_DIR="CenSat_sort_b"
OUT_DIR="giggle_results"

mkdir -p "$OUT_DIR"

for bed in "$TE_DIR"/*.bed.gz
do
    name=$(basename "$bed" .bed.gz)

    echo "Processing $name ..."

    giggle search \
        -i "$INDEX_DIR" \
        -q "$bed" \
        -s \
        > "$OUT_DIR/${name}.giggle.txt"
Done
```

### S4.2.2 combine_result

```r
setwd("/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/Cent_enrichment/giggle_results/")
library(vroom)
filelist <- list.files("./")
n <-length(filelist) 
test<- vroom(filelist[1],id="./",col_names= T,delim="\t")
test<-test[test$fishers_two_tail<0.05,c(2,5)]

for (i in 2:n)
{
  txt_data<-vroom(filelist[i],id="./",col_names= T,show_col_types=F,delim = "\t")
  txt_data<-txt_data[txt_data$fishers_two_tail<0.05,c(2,5)]
  test<- merge(test,txt_data,by="#file",all.x=T,all.y=T)
}
library(stringr)
superfamily<-as.data.frame(str_split(filelist,"\\.",n=3,simplify=TRUE))$V1
name<-c("TF",as.vector(superfamily))
colnames(test)<-name
write.table(test,"../giggle_all_results.txt",row.names = F,sep="\t")
```

### S4.2.3 plotting

```r

# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
data<-read.table("giggle_all_results.txt")
data2<-data[data$V8>0,c(1,4,9)]

data_long<-data2
colnames(data_long)<-c("Cent","OR_value","Cluster")
# # 转换为长格式
# data_long <- data2 %>%
#   pivot_longer(
#     cols = -TF,
#     names_to = "Cluster",
#     values_to = "OR_value"
#   )

# 将TF列转换为因子，保持原始顺序
data_long$Cent <- factor(data_long$Cent, levels = c("divergent_HOR","monomeric_alphaSat",
                                                    "HSat3",
                                                    "betaSat","gammaSat","censat",
                                                    "ct"))

# 将Cluster列转换为因子，保持原始顺序
cluster_order <- c("non-Primates_PiggyBac","Primates_MULE-MuDR","Primates_ERVL",
                   "Simiiformes_PiggyBac","Simiiformes_L1","Simiiformes_ERV1",
                   "Simiiformes_ERVL","Simiiformes_ERVK",
                   "Catarrhini_L1","Catarrhini_ERV1","Catarrhini_ERVK",
                   "Hominoidea_L1","Hominoidea_ERV1","Hominoidea_ERVK","Hominoidea_SVA"
                   )
data_long$Cluster <- factor(data_long$Cluster, levels = cluster_order)

library(ggplot2)
# 创建点图
p<-ggplot(data_long, aes(x = Cluster, y = Cent, size = log(OR_value), color = log(OR_value))) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(
    name = "",
    range = c(2, 8),  # 点的大小范围
    #breaks = c(0.00000000001, 0.01, 0.5, 1, 2),  # 自定义刻度
    #labels = c("2.29e-11", "0.02", "0.5", "1", "2")  # 自定义标签
  ) +
  scale_color_gradient2(
    name = "OR Value",
    low = "blue",
    mid = "white",
    high = "#891a1b",
    #midpoint = 1,  # 以OR=1为中心
    #breaks = c(0.00000000001, 0.01, 0.5, 1, 2),
    #labels = c("2.29e-11", "0.02", "0.5", "1", "2")
  ) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 20, 10), # 增加底部边距以便x轴标签显示完整
    legend.position = "right",
    legend.box = "vertical",
  ) #+
  # # 为点添加文本标签，只显示非常小或非1的值
  # geom_text(
  #   data = subset(data_long, OR_value < 0.1 | OR_value > 1.5),
  #   aes(label = ifelse(OR_value < 0.001, 
  #                      format(OR_value, scientific = TRUE, digits = 2),
  #                      round(OR_value, 2))),
  #   size = 3,
  #   vjust = -1.5,
  #   color = "black",
  #   show.legend = FALSE
  # )

p

pdf("cent-enrich.pdf",width =8,height =4)
p
dev.off()
```
