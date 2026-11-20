# S11

## S11.1 subregion

```r
#workplace:F:\tmp\07-31-figure2重新画图\res.gene
#file:chipseeker.R
############################################################
## 使用 ChIPseeker 注释 BED8 格式的转座子区间
## 默认参考基因组：hg38
############################################################

## 1. 安装需要的包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "ChIPseeker",
  "GenomicRanges",
  "IRanges",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

## 2. 加载包
library(ChIPseeker)
library(GenomicRanges)
library(IRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

## 3. 输入、输出文件
input_bed <- "../Hominoidea.v3.genes.centre.TE.HDG.bed"
output_file <- "Hominoidea.v3.genes.centre.HDG.TE_ChIPseeker_annotation.tsv"

## 4. 读取 BED 文件
bed_data <- read.table(
  input_bed,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)

## 检查列数
if (ncol(bed_data) < 8) {
  stop("输入文件少于8列，请检查BED文件格式。")
}

## 如果文件多于8列，只处理前8列
bed_data <- bed_data[, 1:8]

colnames(bed_data) <- c(
  "chr",
  "bed_start",
  "bed_end",
  "TE_subfamily",
  "score",
  "strand",
  "TE_family",
  "lineage"
)

## 添加唯一编号，便于检查注释前后的对应关系
bed_data$TE_ID <- paste0("TE_", seq_len(nrow(bed_data)))

## 5. 检查坐标和链方向
bed_data$bed_start <- as.integer(bed_data$bed_start)
bed_data$bed_end <- as.integer(bed_data$bed_end)

if (anyNA(bed_data$bed_start) || anyNA(bed_data$bed_end)) {
  stop("第二列或第三列存在无法转换为整数的坐标。")
}

if (any(bed_data$bed_end <= bed_data$bed_start)) {
  stop("存在终止位置小于或等于起始位置的区间。")
}

## 非标准链方向设为 *
bed_data$strand[!bed_data$strand %in% c("+", "-", "*")] <- "*"

## 6. 转换为 GRanges
##
## BED坐标：
##   start：0-based
##   end：右侧开放
##
## GRanges坐标：
##   start：1-based
##   end：闭区间
##
## 因此 BED start 需要加1，而 BED end保持不变
te_gr <- GRanges(
  seqnames = bed_data$chr,
  ranges = IRanges(
    start = bed_data$bed_start + 1,
    end = bed_data$bed_end
  ),
  strand = bed_data$strand
)

## 将原始转座子信息保存到 GRanges metadata
mcols(te_gr)$TE_ID <- bed_data$TE_ID
mcols(te_gr)$TE_subfamily <- bed_data$TE_subfamily
mcols(te_gr)$score <- bed_data$score
mcols(te_gr)$TE_family <- bed_data$TE_family
mcols(te_gr)$lineage <- bed_data$lineage
mcols(te_gr)$original_bed_start <- bed_data$bed_start
mcols(te_gr)$original_bed_end <- bed_data$bed_end

## 7. 设置 hg38 注释数据库
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## 8. 使用 ChIPseeker 注释
##
## tssRegion = c(-3000, 3000)：
## 将TSS上游3 kb至下游3 kb定义为启动子区域
##
## level = "gene"：
## 基于基因层面确定最近基因
peak_annotation <- annotatePeak(
  peak = te_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db",
  level = "gene",
  verbose = FALSE
)

## 9. 转换为数据框
annotation_df <- as.data.frame(peak_annotation)

## 10. 整理输出列
##
## ChIPseeker输出的start是GRanges的1-based坐标，
## 这里优先保留原始BED的0-based坐标
output_columns <- c(
  "seqnames",
  "original_bed_start",
  "original_bed_end",
  "TE_subfamily",
  "score",
  "strand",
  "TE_family",
  "lineage",
  "TE_ID",
  "annotation",
  "geneChr",
  "geneStart",
  "geneEnd",
  "geneLength",
  "geneStrand",
  "geneId",
  "transcriptId",
  "distanceToTSS",
  "ENTREZID",
  "ENSEMBL",
  "SYMBOL",
  "GENENAME"
)

## 不同版本的注释包输出列可能略有差异
## 只提取实际存在的列
output_columns <- output_columns[
  output_columns %in% colnames(annotation_df)
]

annotation_output <- annotation_df[, output_columns, drop = FALSE]

## 重命名染色体和原始坐标列
rename_vector <- c(
  seqnames = "chr",
  original_bed_start = "start",
  original_bed_end = "end"
)

for (old_name in names(rename_vector)) {
  if (old_name %in% colnames(annotation_output)) {
    colnames(annotation_output)[
      colnames(annotation_output) == old_name
    ] <- rename_vector[[old_name]]
  }
}

## 11. 输出注释结果
write.table(
  annotation_output,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("注释完成：", output_file, "\n")
cat("输入区间数量：", nrow(bed_data), "\n")
cat("输出注释数量：", nrow(annotation_output), "\n")

# plotAnnoBar(peakAnno)
# vennpie(peakAnno)
plotAnnoPie(peak_annotation)
# plotDistToTSS(peakAnno)

```

## S11.2 cCRE

```r
#workplace:F:\tmp\08-04-fig2-CNTN1
#cCREpieplot.R
setwd("F://08-04")
data<-unique(read.table("summary.all.cCRE.Hominoidea.v3.genes.centre.TE.txt"))#188
length(unique(data$V7))#161
library(ggplot2)

# 构建频数数据框（不含 label）
tab_df <- as.data.frame(table(data$V3))
colnames(tab_df) <- c("cCRE", "Frequency")

library(ggplot2)
#library(ggrepel)   # 安装：install.packages("ggrepel")

#library(ggplot2)
library(dplyr)

library(ggplot2)
library(ggrepel)   # 安装：install.packages("ggrepel")

# 先计算每个分组的累计比例，用于标签定位（在饼图中方便）
tab_df <- tab_df %>%
  mutate(prop = Frequency / sum(Frequency),
         cumsum = cumsum(prop),
         midpoint = cumsum - prop / 2,   # 每个扇区的中心弧度位置（0~1）
         label = paste0(cCRE, "\n", round(prop * 100, 1), "%"))
tab_df$cCRE<-factor(tab_df$cCRE,levels = c("CTCF-only,CTCF-bound","PLS","PLS,CTCF-bound",
                                           "DNase-H3K4me3","DNase-H3K4me3,CTCF-bound",
                                           "pELS","pELS,CTCF-bound",
                                           "dELS","dELS,CTCF-bound"))


p <- ggplot(tab_df, aes(x = "", y = Frequency, fill = cCRE)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Paired")+
  #scale_fill_manual(values = colors_set1)+
  geom_text_repel(
    data = tab_df,
    aes(x = 1.5, y = midpoint * sum(Frequency), label = label),
    nudge_x = 0.3,          # 向外再推移一点
    segment.color = "grey50", # 引线颜色
    segment.size = 0.3,
    direction = "y",        # 主要沿垂直方向调整
    hjust = 0.5,
    size = 3.5,
    show.legend = FALSE
  )
p
pdf("ENCODE.161loci.Hominoidea.v3.genes.centre.TE.pieplot.pdf",width =8,height=8)
p
dev.off()
```





