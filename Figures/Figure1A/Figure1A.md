> Preanalysis and plot for TE evolution across ages in genome.

# A- classification
## A.1 reannotation
「Dfam release 3.9 (March 2025)」
> ref：https://github.com/Dfam-consortium/FamDB?tab=readme-ov-file

* download dependency
> location:
> 
> https://www.dfam.org/releases/current/families/FamDB/dfam39_full.7.h5.gz
> https://www.dfam.org/releases/current/families/FamDB/dfam39_full.0.h5.gz

* install h5py
```shell pip3 install --user h5py```

* download famdb.py
https://github.com/Dfam-consortium/FamDB/archive/refs/tags/2.0.2.tar.gz

* retrive classification
```shell
./famdb.py -i ./h5  families -f  fasta_name -a -d -c --include-class-in-name 'Homo sapiens' >homo.sapiens.TE.txt
cat homo.sapiens.TE.txt |grep ">" |sed 's/_#\|#\| @/\t/g' |sed 's/\//\t/g' |sed 's/>//g' |sed 's/\[[^]]*\]//g' >homo.sapiens.TE.clades.txt
```
And then do some filtering.


## A.2 
> 1. annotate the ages of each TE locus according to the gtf file.
> 2. split the beds into different age(clade)/family combinations.

```shell
bash clade.sh newclade.annot.csv "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/hs1_rmsk_TE.gtf"
mkdir clade_family
python split_clade_family.py
```
* 1.code: clash.sh
* 2.file: newclade.annot.csv
* 3.code: split_clade_family.py

## A.3 Plot input:processed_data.csv
```r

library(dplyr)
library(tidyr)
library(ggplot2)

processed_data<-read.table("processed_data.csv") #已经预先计算好了比例的信息；

# 创建完整的数据网格（确保每个family在所有newclade中都有记录）
complete_data <- expand.grid(
  family = unique(processed_data$family),
  newclade = c("non-Primates", "Primates", "Simiiformes", "Catarrhini", "Hominoidea"),
  stringsAsFactors = FALSE
)

# 合并原始数据并填充缺失值为0
plot_data <- processed_data %>%
  select(family, newclade, percentage) %>%
  right_join(complete_data, by = c("family", "newclade")) %>%
  mutate(percentage = ifelse(is.na(percentage), 0, percentage))

# 创建排序权重数据
weight_data <- plot_data %>%
  pivot_wider(
  	names_from = newclade,
    values_from = percentage,
    values_fill = 0
    ) %>%
    mutate(
    # 按指定顺序赋予递减权重
    sort_weight = `non-Primates` * 1e8 +
    Primates * 1e6 +
    Simiiformes * 1e4 +
    Catarrhini * 1e2 +
    Hominoidea
    ) %>%
    arrange(desc(sort_weight))
# 获取排序后的family顺序
family_order <- weight_data$family

# 应用排序因子
plot_data <- plot_data %>%
mutate(
family = factor(family, levels = family_order),
newclade = factor(newclade, levels = c("non.Primates", "Primates", "Simiiformes", "Catarrhini", "Hominoidea"))
)
common_levels <- c(
"5S.Deu.L2", "CR1", "Dong.R4", "Gypsy", "I-Jockey", "Kolobok", "L1-Tx1", "L2",
"Long.Terminal.Repeat.Element", "MIR", "PIF.Harbinger", "RTE.BovB", "RTE.X",
"TcMar", "TcMar.Pogo", "TcMar.Tc1", "TcMar.Tc2", "Transposase", "hAT", "hAT.Ac",
"hAT.Blackjack", "hAT.Tag1", "hAT.Tip100", "tRNA", "tRNA.Deu", "tRNA.RTE",
"hAT.Charlie", "ERVL", "ERVL.MaLR", "Merlin", "TcMar.Tigger", "L1", "TcMar.Mariner",
"MULE.MuDR", "ERV1", "PiggyBac", "Alu", "ERVK", "SVA"
)
shared_scale <- function() {
list(
scale_x_discrete(limits = common_levels),
theme_minimal() # 可添加统一主题
)
}

# 绘制堆叠条形图
p3<-ggplot(plot_data, aes(x = family, y = percentage, fill = newclade)) +
geom_col(position = "stack") +
labs(
x = "Family",  # 更准确的轴标签
y = "Percentage of Subfamilies (%)",
fill = "Clade"  # 更准确的图例标题
) +
scale_fill_manual(values = c(
"non.Primates" = "#812e29",
"Primates" = "#549357",
"Simiiformes" = "#e7d968",
"Catarrhini" = "#8b6a93",
"Hominoidea" = "#92cbe1"
)) +
theme_classic() +
theme(
axis.text.x = element_text(angle = 45, hjust = 1),
plot.title = element_text(hjust = 0.5)
)

p3

pdf("subfamily.per.2.pdf",width =12 ,height = 3)
print(p3)
dev.off()
```


