# B. The distribution of transposon in human genome(T2T)

## B.1 split the genome into 500kb bin and count the transposon in each age/family combination.
 

```python
import os
os.chdir("/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# 定义常量
BIN_SIZE = 500000
TE_FAMILIES = ['L1', 'ERV1', 'Alu', 'ERVK', 'SVA']
CLADES_ORDER = ['non-Primates', 'Primates', 'Simiiformes', 'Catarrhini', 'Hominoidea']
CHROMOSOMES_ORDER = [f'chr{i}' for i in range(1, 23)] + ['chrX']  #chromosome 1:23,X 

def create_genomic_bins(chrom_length_file):
    """创建基因组区间(bins)"""
    bins = []
    
    # 读取染色体长度文件
    with open(chrom_length_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
                
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            
            # 跳过chrY
            if chrom == 'chrY':
                continue
                
            # 创建4Mb的区间
            for bin_start in range(start, end, BIN_SIZE):
                bin_end = min(bin_start + BIN_SIZE - 1, end)
                bins.append({
                    'chrom': chrom,
                    'start': bin_start,
                    'end': bin_end,
                    'bin_length': bin_end - bin_start + 1
                })
    
    return pd.DataFrame(bins)

def read_and_filter_bed(bed_file):
    """读取并过滤BED文件"""
    # 读取BED文件
    bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'clade', 'te_family']
    df = pd.read_csv(bed_file, sep='\t', header=None, names=bed_cols)
    
    # 过滤：只保留指定的TE家族，排除chrY
    df = df[df['te_family'].isin(TE_FAMILIES)]
    df = df[df['chrom'] != 'chrY']
    
    # 标准化clade名称（确保与指定的clade顺序匹配）
    # 注意：这里假设bed文件中的clade名称与指定名称一致
    df['clade'] = df['clade'].astype(str)
    
    return df

def assign_te_to_bins(te_df, bins_df):
    """使用IntervalIndex进行高效区间匹配"""
    results = []
    
    # 按染色体分组处理
    for chrom in te_df['chrom'].unique():
        chrom_tes = te_df[te_df['chrom'] == chrom]
        chrom_bins = bins_df[bins_df['chrom'] == chrom]
        
        if chrom_tes.empty or chrom_bins.empty:
            continue
        
        # 为当前染色体的bins创建IntervalIndex
        intervals = pd.IntervalIndex.from_arrays(
            chrom_bins['start'],
            chrom_bins['end'],
            closed='both'
        )
        
        # 为每个TE找到所在的区间
        # 使用get_indexer可以批量查找
        indices = intervals.get_indexer(chrom_tes['start'].values)
        
        # 只保留找到匹配的TE（indices != -1）
        mask = indices != -1
        matched_indices = indices[mask]
        matched_tes = chrom_tes.iloc[np.where(mask)[0]]
        
        # 构建结果
        for idx, (_, te_row) in zip(matched_indices, matched_tes.iterrows()):
            bin_row = chrom_bins.iloc[idx]
            results.append({
                'chrom': chrom,
                'bin_start': bin_row['start'],
                'bin_end': bin_row['end'],
                'clade': te_row['clade'],
                'te_family': te_row['te_family']
            })
    
    return pd.DataFrame(results)
    

def calculate_density_matrix(bins_df, te_bins_df):
    """计算密度矩阵"""
    # 创建多级索引
    bin_keys = [(row['chrom'], row['start'], row['end']) for _, row in bins_df.iterrows()]
    
    # 创建空的密度矩阵
    density_data = []
    
    # 对每个clade和TE family组合计算密度
    for clade in CLADES_ORDER:
        row_data = {}
        
        # 筛选当前clade的TE
        clade_te = te_bins_df[te_bins_df['clade'] == clade]
        
        for te_family in TE_FAMILIES:
            # 筛选当前TE family
            te_df = clade_te[clade_te['te_family'] == te_family]
            
            for _, bin_row in bins_df.iterrows():
                chrom, bin_start, bin_end = bin_row['chrom'], bin_row['start'], bin_row['end']
                
                # 生成列名
                col_name = f"{chrom}:{bin_start}-{bin_end}_{te_family}"
                
                # 计数
                count = len(te_df[
                    (te_df['chrom'] == chrom) & 
                    (te_df['bin_start'] == bin_start) & 
                    (te_df['bin_end'] == bin_end)
                ])
                
                # 计算密度
                density = count / BIN_SIZE
                row_data[col_name] = density
        
        density_data.append(row_data)
    
    # 创建DataFrame
    density_df = pd.DataFrame(density_data, index=CLADES_ORDER)
    
    return density_df

def sort_columns(density_df):
    """按照要求对列进行排序"""
    # 提取列信息
    columns_info = []
    for col in density_df.columns:
        # 解析列名：chrom:start-end_te
        chrom_part, te_part = col.split('_')
        chrom, range_part = chrom_part.split(':')
        start, end = map(int, range_part.split('-'))
        columns_info.append({
            'col_name': col,
            'chrom': chrom,
            'start': start,
            'te_family': te_part
        })
    
    # 转换为DataFrame以便排序
    cols_df = pd.DataFrame(columns_info)
    
    # 为染色体添加排序键
    chrom_order = {chrom: i for i, chrom in enumerate(CHROMOSOMES_ORDER)}
    cols_df['chrom_order'] = cols_df['chrom'].map(chrom_order)
    
    # 为TE family添加排序键
    te_order = {te: i for i, te in enumerate(TE_FAMILIES)}
    cols_df['te_order'] = cols_df['te_family'].map(te_order)
    
    # 排序：先按染色体，再按起始位置，最后按TE family
    cols_df = cols_df.sort_values(['chrom_order', 'start', 'te_order'])
    
    # 获取排序后的列名
    sorted_cols = cols_df['col_name'].tolist()
    
    # 重新排列列
    return density_df[sorted_cols]

def save_to_csv(density_df, output_file):
    """保存为CSV文件"""
    density_df.to_csv(output_file)
    print(f"密度矩阵已保存到: {output_file}")

def create_heatmap(density_df, output_file='te_density_heatmap.png'):
    """创建热图"""
    # 设置图形大小
    plt.figure(figsize=(20, 8))
    
    # 创建掩码：对于值为0的位置（即不存在数据）
    mask = density_df.values == 0
    
    # 创建热图
    ax = sns.heatmap(
        density_df.values,
        cmap='YlOrRd',  # 黄色到红色渐变
        linewidths=0.5,
        linecolor='gray',
        cbar_kws={'label': 'Density (counts per 4Mb)'},
        mask=mask,
        square=False
    )
    
    # 设置坐标轴标签
    plt.yticks(np.arange(len(CLADES_ORDER)) + 0.5, CLADES_ORDER, rotation=0)
    
    # 设置x轴标签（由于列太多，只显示部分标签）
    n_cols = len(density_df.columns)
    if n_cols > 50:  # 如果列太多，只显示部分标签
        step = max(1, n_cols // 20)
        xtick_positions = np.arange(0, n_cols, step)
        xtick_labels = [density_df.columns[i] for i in xtick_positions]
        plt.xticks(xtick_positions, xtick_labels, rotation=90, fontsize=8)
    else:
        plt.xticks(np.arange(n_cols) + 0.5, density_df.columns, rotation=90, fontsize=8)
    
    plt.title('TE Density Heatmap by Evolutionary Age and Genomic Region')
    plt.tight_layout()
    
    # 保存图像
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"热图已保存到: {output_file}")
    plt.show()
def main():
    # 文件路径
    chrom_length_file = "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/hs1.genome"
    bed_file = "/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/result2_by_newclade_family.bed"  #input files
    output_csv = "te_density_matrix.csv" #output
    output_heatmap = "te_density_heatmap.png"
    
    print("步骤1: 创建基因组区间...")
    bins_df = create_genomic_bins(chrom_length_file)
    print(f"创建了 {len(bins_df)} 个基因组区间")
    
    print("\n步骤2: 读取并过滤BED文件...")
    te_df = read_and_filter_bed(bed_file)
    print(f"读取了 {len(te_df)} 个TE元件")
    
    print("\n步骤3: 将TE分配到基因组区间...")
    te_bins_df = assign_te_to_bins(te_df, bins_df)
    print(f"分配完成，共 {len(te_bins_df)} 个分配记录")
    
    print("\n步骤4: 计算密度矩阵...")
    density_df = calculate_density_matrix(bins_df, te_bins_df)
    
    print("\n步骤5: 对列进行排序...")
    density_df = sort_columns(density_df)
    
    print("\n步骤6: 保存到CSV文件...")
    save_to_csv(density_df, output_csv)
    
    print("\n步骤7: 创建热图...")
    create_heatmap(density_df, output_heatmap)
    
    # 打印一些统计信息
    print("\n=== 统计信息 ===")
    print(f"进化年龄类别: {len(CLADES_ORDER)} 个")
    print(f"TE家族: {len(TE_FAMILIES)} 个")
    print(f"基因组区间数量: {len(bins_df)} 个")
    print(f"密度矩阵形状: {density_df.shape}")
    
    # 显示每个clade的总密度
    print("\n每个进化年龄类别的平均密度:")
    for clade in CLADES_ORDER:
        avg_density = density_df.loc[clade].mean()
        print(f"  {clade}: {avg_density:.6f}")

if __name__ == "__main__":
    main()
```

## B.2 plot the heatmap and annotation the bins
### B.2.1 prepare column annotation

#### column ==> bed3
```shell
awk -F, 'NR>1 {                                                                                         gsub(/"/, "", $1);                # 去掉首字段的引号
    split($1, a, ":");                # 分割染色体和坐标区间
    split(a[2], b, "-");              # 分割起始和结束
    print a[1], b[1], b[2]          # 输出 BED 三列：chr  start  end
}' "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/tmp/l1_data_0to1_no0_row.csv" >bins.bed

awk 'BEGIN{OFS="\t"} {                                                                                                                                        print $1, $2, $3, "bin_" NR
 }' bins.bed >cluster_assignments.named.bed #染色体和坐标区间
```

#### GC%

```shell
python gc_at_calculator.py l1_data_0to1_no0_row.csv "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/hs1.fa" annot.gc.csv
结果已保存至: annot.gc.csv
```

#### CpG fraction


```shell
awk 'BEGIN{OFS="\t"}
 $1 !~ /^(track|browser|#)/ && NF >= 3 {
     print $1,$2,$3
 }' "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/CpG_hs1.bed" |  LC_ALL=C sort -k1,1 -k2,2n >cpgIsland.clean.sorted.bed

bedtools coverage -a "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/tmp/cluster_assignments.named.bed" -b "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/tmp/cpgIsland.clean.sorted.bed" >cpgIsland_coverage.bed
```

#### Censat fraction

```shell
awk 'BEGIN{
    FS=OFS="\t"
}
NR==1{
    print $0, "family", "class"
    next
}
{
    # 提取第4列第一个 "_" 前的类别
    split($4, a, "_")
    family=a[1]

    # hor进一步区分active和inactive
    if(family=="hor"){
        if($4 ~ /L\)$/){
            family="active_HOR"
        }else{
            family="inactive_HOR"
        }
    }
    else if(family=="dhor"){
        family="divergent_HOR"
    }
    else if(family=="mon"){
        family="monomeric_alphaSat"
    }
    else if(family=="hsat1A"){
        family="HSat1A"
    }
    else if(family=="hsat1B"){
        family="HSat1B"
    }
    else if(family=="hsat2"){
        family="HSat2"
    }
    else if(family=="hsat3"){
        family="HSat3"
    }
    else if(family=="bsat"){
        family="betaSat"
    }
    else if(family=="gsat"){
        family="gammaSat"
    }

    # 根据family定义class
    if(family=="active_HOR" ||
       family=="inactive_HOR" ||
       family=="divergent_HOR" ||
       family=="monomeric_alphaSat"){
        class="alpha-satellite"
    }
    else if(family=="HSat1A" ||
            family=="HSat1B" ||
            family=="HSat2" ||
            family=="HSat3"){
        class="Classical human satellite"
    }
    else if(family=="betaSat" ||
            family=="gammaSat" ||
            family=="censat"){
        class="Other satellite"
    }
}' CenSat_hs1.bed >CenSat_hs1.named.bed
```


```R
A="/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/cluster_assignments.named.bed"
B="/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/CenSat_hs1.named.bed"

awk 'BEGIN{FS="\t"} NR>1 {print $11}' "$B" | sort -u | while read family; do     echo "Processing: $family" >&2;      awk -v fam="$family" '
    BEGIN{FS=OFS="\t"}
    NR>1 && $11==fam
    ' "$B"     | bedtools coverage -a "$A" -b -     | awk -v fam="$family" 'BEGIN{OFS="\t"} {
        print $0, fam
    }'     >> coverage_by_family.long.bed; done


awk 'BEGIN{FS="\t"} NR>1 {print $12}' "$B" | sort -u | while read class; do     awk -v cls="$class" '                   BEGIN{FS=OFS="\t"}
    NR>1 && $12==cls
    ' "$B"     | bedtools coverage -a "$A" -b -     | awk -v cls="$class" 'BEGIN{OFS="\t"} {print $0, cls}'; done > coverage_by_class.long.bed


bedtools coverage -a "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/tmp/cluster_assignments.named.bed" -b "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/tmp/CenSat_hs1.named.bed" >CenSat_coverage.bed
```

#### hs1 unique


```shell
bedtools coverage -a "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/tmp/cluster_assignments.named.bed" -b "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/tmp/hs1_unique.bed"
```

#### summarize fraction


```shell
#修改merge_bin_coverage.py中关于输入文件的文件名为我们想要整合的fraction的文件名
python merge_bin_coverage.py
```

#### organize all the annotations


```shell
annotation_col <- read.csv(
 "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/annot.gc.csv",
  header = TRUE,
  sep=",",
  #quote = "",
  check.names = FALSE,
  
)
dim(annotation_col)
head(annotation_col)
rownames(annotation_col)<-annotation_col[,1]
###2026-07-29在"cluster_assignments_annotation_final.csv"文件的基础上继续往下做；添加其它注释项目
annotation_col$chromosome <- sub(":.*", "", rownames(annotation_col))
head(annotation_col)
#                                    value chromosome
# chr1:100000001-100500000              10       chr1
# chrX:500-1000                         20        chrX
# chr2:300000-400000                    30       chr2
colnames(annotation_col)[1]<-"region_id"
annotation_col2 <- read.delim(
  "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/bins.genome_annotation.covered_fraction.tsv",
  header = TRUE,
  #row.names = 1,
  check.names = FALSE
)

head(annotation_col2)
annotation_col3<-merge(annotation_col,annotation_col2,by="region_id")
head(annotation_col3)
rownames(annotation_col3)<-annotation_col3$region_id
annotation_col4<-annotation_col3[,c(-1,-2)]#把端粒也去掉，这个范围太小了；
colnames(annotation_col4)<-c("AT%","Chromosome","CpG%","centromere%","CHM13_unique")
head(annotation_col4)
#chr_names <- paste0("chr", c(1:22, "X"))
annotation_col4$Chromosome<-factor(annotation_col4$Chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))

```

### B.2.2 plot
```r

setwd("/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/")
library(data.table)

#使用fread读取，速度比read.csv快5-10倍
data <- fread("te_density_matrix.csv", header = TRUE)

#如果需要第一列作为行名
row_names <- data[[1]]
data <- as.matrix(data[, -1])
rownames(data) <- row_names

library(data.table)
library(pheatmap)
library(dplyr)
library(tidyr)  # pivot_longer和pivot_wider在tidyr包中

#使用基础R的rownames处理
data_tibble <- as.data.frame(data)
data_tibble$TE_Family <- rownames(data_tibble)

#将行名列移到第一列
data_tibble <- data_tibble[, c("TE_Family", setdiff(colnames(data_tibble),"TE_Family"))]

#继续后续操作
long_format <- data_tibble %>%
  pivot_longer(cols = -TE_Family, 
               names_to = "Sample", 
               values_to = "Expression")
 
 #从样本名中提取类型信息
long_format <- long_format %>%
  mutate(Type = case_when(
    grepl("_L1$", Sample) ~ "L1",
    grepl("_ERVK$", Sample) ~ "ERVK",
    grepl("_Alu$", Sample) ~ "Alu",
    grepl("_ERV1$", Sample) ~ "ERV1",
    grepl("_SVA$", Sample) ~ "SVA",
    TRUE ~ "Other"
  )) %>%
  mutate(Sample_Base = gsub("_(L1|ERVK|Alu|ERV1|SVA)$", "", Sample))

#创建新的行标识符
long_format <- long_format %>%
  mutate(TE_Family_Type = paste0(TE_Family, "_", Type))

#转换回宽格式（行是TE_Family_Type，列是样本）
wide_format <- long_format %>%
  select(TE_Family_Type, Sample_Base, Expression) %>%
  pivot_wider(names_from = Sample_Base, values_from = Expression)

#转换为矩阵用于热图
final_matrix <- as.matrix(wide_format[, -1])
rownames(final_matrix) <- wide_format$TE_Family_Type

#假设 final_matrix 是你的矩阵
row_min <- apply(final_matrix, 1, min, na.rm = TRUE)
row_min

l1_data_0to1_no0_row <- final_matrix[rowSums(final_matrix) > 0,] #只删除行全为0的行
dim(l1_data_0to1_no0_row)

#scale 0-1,分subfamily scale比较对
l1_data_0to1_no0_row_final <- t(apply(l1_data_0to1_no0_row, 1, function(row) {
  if (max(row) == min(row)) {
    # 如果所有值相同，返回0.5，如果全为0呢？
    return(rep(0, length(row)))
  } else {
    # 缩放到0-1范围
    return((row - min(row)) / (max(row) - min(row)))
  }
}))

write.csv(colnames(l1_data_0to1_no0_row_final),row.names=F,"l1_data_0to1_no0_row_final.csv")

#得到这个csv文件之后，随后计算GC含量
python gc_at_calculator.py l1_data_0to1_no0_row.csv "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/hs1.fa" annot.gc.csv
#结果已保存至: annot.gc.csv

annotation_col <- read.csv(
 "/home/xxzhang/workplace/jupyter/annot.gc.t2t.1014.csv",
  header = TRUE,
  sep=",",
  #quote = "",
  check.names = FALSE,
  
)
dim(annotation_col)

head(annotation_col)
rownames(annotation_col)<-annotation_col[,1]

###2026-07-29在"cluster_assignments_annotation_final.csv"文件的基础上继续往下做；添加其它注释项目
annotation_col$chromosome <- sub(":.*", "", rownames(annotation_col))
head(annotation_col)
#                                    value chromosome
# chr1:100000001-100500000              10       chr1
# chrX:500-1000                         20        chrX
# chr2:300000-400000                    30       chr2

colnames(annotation_col)[1]<-"region_id"

###bins.genome_annotation.covered_fraction.tsv文件来源于bedtools coverage
annotation_col2 <- read.delim(
  "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/Annotation/bins.genome_annotation.covered_fraction.tsv",
  header = TRUE,
  #row.names = 1,
  check.names = FALSE
)

dim(annotation_col2)

annotation_col3<-merge(annotation_col,annotation_col2,by="region_id")
head(annotation_col3)

rownames(annotation_col3)<-annotation_col3$region_id
annotation_col4<-annotation_col3[,c(-1,-2)]#把端粒也去掉，这个范围太小了；

#colnames(annotation_col4)<-c("AT%","Chromosome","CpG%","centromere%","CHM13_unique")
colnames(annotation_col4)<-c("AT%","Chromosome","CpG%","CHM13_unique")
head(annotation_col4)

#chr_names <- paste0("chr", c(1:22, "X"))
annotation_col4$Chromosome<-factor(annotation_col4$Chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))

#####################################################################################################################

ann_colors=list(x=c("1"='#A6CEE3',"2"='#1F78B4'))
library(pheatmap)

highlight_label <- rep(
  "Other",
  ncol(l1_data_0to1_no0_row)
)

names(highlight_label) <- colnames(l1_data_0to1_no0_row)

for (group_name in names(highlight_groups)) {
  highlight_label[
    colnames(l1_data_0to1_no0_row) %in%
      highlight_groups[[group_name]]
  ] <- group_name
}

# highlight_annotation <- data.frame(
#   Highlight = factor(
#     highlight_label,
#     levels = c("Group1", "Group2", "Group3", "Other")
#   )
# )

#rownames(highlight_annotation) <-
 # colnames(l1_data_0to1_no0_row)


## =========================
## annotation_col6 merge
## =========================

## 保证原注释的顺序与矩阵列顺序完全一致
annotation_col4 <- annotation_col4[
  colnames(l1_data_0to1_no0_row),
  drop = FALSE
]

annotation_col_new <- cbind(
  annotation_col4#,
  #highlight_annotation
)


## =========================
## 6. 设置注释颜色
## =========================

## 如果此前已有 ann_colors，直接增加 Highlight
ann_colors$Highlight <- c(
  Group1 = "#00A6D6",
  Group2 = "#F39B7F",
  Group3 = "#7E57C2",
  Other  = "#FFFFFF"
) #trash

chr_colors <- c(
  "chr1"  = "#E6194B", #
  "chr2"  = "#377db8", #
  "chr3"  = "#4eae49", #
  "chr4"  = "#974ea3", #
  "chr5"  = "#f19301", #
  "chr6"  = "#f781be", #
  "chr7"  = "#bb9dcc", #
  "chr8"  = "#a75529", #
  "chr9"  = "#53afe3", #
  "chr10" = "#212e75", #
  "chr11" = "#1c9e78", #
  "chr12" = "#b1df89", #
  "chr13" = "#e2bd00", #
  "chr14" = "#fb9a97", #
  "chr15" = "#e72989", #
  "chr16" = "#910141", #
  "chr17" = "#04cdd1", #
  "chr18" = "#a6cee2", #
  "chr19" = "#ce1260", #
  "chr20" = "#5e4fa1", #
  "chr21" = "#8ba77a", #
  "chr22" = "#00441b", #
  "chrX"  = "#dedc02" #
)

ann_colors$Chromosome <- chr_colors

## =========================
## plot
## =========================

p <- pheatmap(
  l1_data_0to1_no0_row_final, #用的是这个table；

  main = "",

  color = colorRampPalette(
    c("#0a0c0d", "#f5eb0d", "#eb2627")
  )(50),

  show_rownames = TRUE,
  #show_colnames = TRUE,
  show_colnames = FALSE,

  fontsize_row = 6,
  fontsize_col = 12,

  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",

  scale = "none",

  annotation_col = annotation_col_new,
  annotation_colors = ann_colors,

  annotation_names_col = TRUE,
  annotation_legend = TRUE,

  cluster_rows = TRUE,
  cluster_cols = TRUE,

  border_color = NA
)

```




























