# S11

## S11.1 subregion

```r
#input:Hominoidea.v3.genes.centre.TE.HDG.bed
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
#这里是用时最长的；
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
#总的输入区间数：4471
cat("输出注释数量：", nrow(annotation_output), "\n")
#chipseek注释区间数：4471，也就是全部注释完成了，可以据此画一张图了。
#我们自己画，不要那张饼图；
#读取文件：Hominoidea.v3.genes.centre.HDG.TE_ChIPseeker_annotation.tsv

# plotAnnoBar(peakAnno)
# vennpie(peakAnno)
#plotAnnoPie(peak_annotation)
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

## S11.3 plotting，桑基图

```r

#base_dir <- "F://tmp//07-31-figure2重新画图//res.gene//Ploting"

#!/usr/bin/env Rscript
# ============================================================================
# FINAL PAPER FIGURE: All TE loci -> genomic annotation -> cCRE-positive subset
#
# First two layers: ALL unique TE loci in the ChIPseeker annotation file.
# Third layer: ONLY TE loci overlapping a cCRE (NO "no results" node).
# Loci lacking cCRE stop at the second layer.
#
# INPUT: two original user files (no pre-merged intermediate is required).
# OUTPUT: full-view PDF/PNG/SVG, separately magnified cCRE-positive PDF/PNG/SVG,
#         audit tables, and optional interactive HTML (networkD3 if installed).
#
# Usage:
#   Rscript TE_cCRE_final_all_loci_to_161.R /path/to/input_directory
#
# Install once if needed:
#   install.packages(c("dplyr", "stringr", "readr", "ggplot2", "svglite"))
# Optional interactive HTML:
#   install.packages(c("networkD3", "htmlwidgets"))
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else "F://tmp//07-31-figure2重新画图//res.gene//Ploting"
file_anno <- file.path(base_dir, "Hominoidea.v3.genes.centre.HDG.TE_ChIPseeker_annotation.tsv")
file_ccre <- file.path(base_dir, "summary.all.cCRE.Hominoidea.v3.genes.centre.TE.txt")
out_dir <- file.path(base_dir, "TE_cCRE_final_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file_anno) || !file.exists(file_ccre)) {
  stop("Input files not found. Supply the folder containing the two original files.")
}

keys <- c("chr", "start", "end", "TE_subfamily")

# ----- 1. Read and deduplicate the cCRE summary -----------------------------
# There are 11 whitespace-separated fields in the source summary.
# Deduplication key: columns 6-11 (the same TE locus), plus column 3 (label).
raw_ccre <- read.table(
  file_ccre, header = FALSE, sep = "", quote = "", comment.char = "",
  colClasses = "character", stringsAsFactors = FALSE
)
if (ncol(raw_ccre) != 11L) {
  stop("Expected 11 whitespace-separated columns in summary, found ",
       ncol(raw_ccre))
}
colnames(raw_ccre) <- c(
  "cCRE_chr", "cCRE_pos", "cCRE_label", "cCRE_ID", "separator",
  "chr", "start", "end", "TE_subfamily", "TE_lineage", "TE_family_summary"
)
if (!all(raw_ccre$separator == "|")) {
  stop("Unexpected separator in column 5 of the cCRE summary")
}

ccre_dedup <- raw_ccre %>%
  group_by(across(all_of(c(keys, "TE_lineage", "TE_family_summary", "cCRE_label")))) %>%
  summarise(
    cCRE_ID = paste(unique(cCRE_ID), collapse = ";"),
    cCRE_positions = paste(unique(paste(cCRE_chr, cCRE_pos, sep = ":")),
                           collapse = ";"),
    n_original_records = n(),
    .groups = "drop"
  )
write_tsv(ccre_dedup, file.path(out_dir, "01_summary_deduplicated.tsv"))

# ----- 2. Read the ChIPseeker file; split Intron (detail) ---------------------
anno <- read.delim(
  file_anno, sep = "\t", header = TRUE, check.names = FALSE,
  quote = "", colClasses = "character"
)
needed <- c(keys, "TE_family", "annotation")
if (!all(needed %in% names(anno))) {
  stop("Missing columns in ChIPseeker file: ",
       paste(setdiff(needed, names(anno)), collapse = ", "))
}

anno <- anno %>%
  mutate(
    Genomic_annotation = str_trim(str_remove(annotation, "\\s*\\(.*$")),
    Annotation_detail = str_extract(annotation, "\\(.*\\)$")
  )

# A single locus may appear repeatedly because of gene/transcript assignments.
# Check that its broad genomic annotation is consistent before collapsing rows.
conflicts <- anno %>%
  group_by(across(all_of(keys))) %>%
  summarise(n_classes = n_distinct(Genomic_annotation), .groups = "drop") %>%
  filter(n_classes > 1)
if (nrow(conflicts)) {
  write_tsv(conflicts, file.path(out_dir, "WARNING_annotation_conflicts.tsv"))
  stop("Some TE loci have conflicting annotation types; inspect WARNING file.")
}
anno_unique <- anno %>% distinct(across(all_of(keys)), .keep_all = TRUE)

# ----- 3. Merge, classify cCRE, and create per-locus weights ----------------
# Original combined labels such as 'dELS,CTCF-bound' are reduced to 'dELS'.
# Priority: PLS > pELS > dELS > DNase-H3K4me3 > CTCF-only.
collapse_ccre <- function(x) {
  if (is.na(x) || x == "" || x == "no results") return(NA_character_)
  x <- str_replace_all(x, "–", "-")
  if (str_detect(x, "\\bPLS\\b")) return("PLS")
  if (str_detect(x, "\\bpELS\\b")) return("pELS")
  if (str_detect(x, "\\bdELS\\b")) return("dELS")
  if (str_detect(x, fixed("DNase-H3K4me3"))) return("DNase-H3K4me3")
  if (str_detect(x, fixed("CTCF-only"))) return("CTCF-only")
  stop("Unknown cCRE label: ", x)
}

merged <- anno_unique %>%
  left_join(
    ccre_dedup %>%
      select(all_of(keys), cCRE_label, cCRE_ID, cCRE_positions,
             n_original_records),
    by = keys
  ) %>%
  mutate(
    locus_id = paste(chr, start, end, TE_subfamily, sep = ":"),
    cCRE_5class = vapply(cCRE_label, collapse_ccre, character(1)),
    Genomic_class = case_when(
      Genomic_annotation %in% c("5' UTR", "3' UTR") ~ "UTR",
      Genomic_annotation %in% c("Distal Intergenic", "Downstream") ~ "Intergenic",
      is.na(Genomic_annotation) | Genomic_annotation == "" ~ "Other",
      TRUE ~ Genomic_annotation
    )
  ) %>%
  group_by(locus_id) %>%
  mutate(TE_weight = 1 / n()) %>%
  ungroup()

# If a TE has 2 different cCRE labels, each contributes 0.5 TE to the cCRE layer.
# It still contributes exactly 1 TE to the full family -> annotation layer.
check_weights <- merged %>%
  group_by(locus_id) %>%
  summarise(weight = sum(TE_weight), .groups = "drop")
stopifnot(all(abs(check_weights$weight - 1) < 1e-8))

n_all <- n_distinct(merged$locus_id)
n_positive <- n_distinct(merged$locus_id[!is.na(merged$cCRE_5class)])
stopifnot(n_all >= n_positive)
write_tsv(merged, file.path(out_dir, "02_TE_full_merged_with_weights.tsv"))

# ----- 4. Build links: first link ALL, second link POSITIVE ONLY ------------
links_full <- merged %>%
  group_by(TE_family, Genomic_class) %>%
  summarise(value = sum(TE_weight), .groups = "drop") %>%
  transmute(from = TE_family, to = Genomic_class, value,
            group = paste0("FAM:", from))

links_ccre <- merged %>%
  filter(!is.na(cCRE_5class)) %>%
  group_by(Genomic_class, cCRE_5class) %>%
  summarise(value = sum(TE_weight), .groups = "drop") %>%
  transmute(from = Genomic_class, to = cCRE_5class, value,
            group = paste0("CRE:", to))

stopifnot(abs(sum(links_full$value) - n_all) < 1e-8)
stopifnot(abs(sum(links_ccre$value) - n_positive) < 1e-8)

write_tsv(links_full, file.path(out_dir, "03_all_TE_family_to_annotation.tsv"))
write_tsv(links_ccre, file.path(out_dir, "04_positive_annotation_to_cCRE.tsv"))

# Cross-check: each annotation sends no more TE to cCRE than it receives.
mid_in <- links_full %>% group_by(to) %>% summarise(n_in = sum(value), .groups = "drop")
mid_out <- links_ccre %>% group_by(from) %>% summarise(n_out = sum(value), .groups = "drop")
flow_check <- mid_in %>% left_join(mid_out, by = c("to" = "from")) %>%
  mutate(n_out = coalesce(n_out, 0), stopped_here = n_in - n_out)
stopifnot(all(flow_check$stopped_here > -1e-8))
write_tsv(flow_check, file.path(out_dir, "05_annotation_inflow_outflow.tsv"))

# ----- 5. Fixed layout order and colors -------------------------------------
family_order <- c("L1", "ERV1", "ERVK", "Alu", "SVA")
anno_order <- c("Intron", "Exon", "Promoter", "UTR", "Intergenic", "Other")
ccre_order <- c("PLS", "pELS", "dELS", "DNase-H3K4me3", "CTCF-only")

# Append any unexpected new annotation/family instead of silently losing it.
family_order <- c(intersect(family_order, unique(links_full$from)),
                  sort(setdiff(unique(links_full$from), family_order)))
anno_order <- c(intersect(anno_order, unique(links_full$to)),
                sort(setdiff(unique(links_full$to), anno_order)))
ccre_order <- c(intersect(ccre_order, unique(links_ccre$to)),
                sort(setdiff(unique(links_ccre$to), ccre_order)))

fam_colors <- c(
  "L1" = "#35B779", "ERV1" = "#31688E", "ERVK" = "#21918C",
  "Alu" = "#482878", "SVA" = "#B5DE2B"
)
extra_fam <- setdiff(family_order, names(fam_colors))
if (length(extra_fam)) fam_colors[extra_fam] <- "#888888"
ccre_colors <- c(
  "PLS" = "#CE7040", "pELS" = "#A884C4", "dELS" = "#538FBC",
  "DNase-H3K4me3" = "#C3A24E", "CTCF-only" = "#7C8894"
)
fill_colors <- c(
  setNames(unname(fam_colors), paste0("FAM:", names(fam_colors))),
  "ANN" = "#E7E9EC",
  setNames(unname(ccre_colors), paste0("CRE:", names(ccre_colors)))
)

# ----- 6. Draw a true, partially terminating Sankey as vector polygons ------
# No missing-category placeholder is added. Both link widths have the SAME
# quantitative scale; hence the third layer is much smaller than the first two.
# A separate magnified 161-locus view is also saved for legibility.
node_layout <- function(x, order, xcentre, gap, center_at = NULL) {
  tt <- data.frame(name = order, stringsAsFactors = FALSE) %>%
    left_join(x, by = "name") %>%
    mutate(value = coalesce(value, 0)) %>%
    filter(value > 0)
  n <- nrow(tt)
  if (!n) stop("Cannot draw an empty node layer")
  tt$y0 <- cumsum(c(0, head(tt$value + gap, -1)))
  tt$y1 <- tt$y0 + tt$value
  if (!is.null(center_at)) {
    delta <- center_at - (min(tt$y0) + max(tt$y1)) / 2
    tt$y0 <- tt$y0 + delta
    tt$y1 <- tt$y1 + delta
  }
  tt$x <- xcentre
  tt
}

# Build nonoverlapping source/target intervals for each link.
position_links <- function(links, left_nodes, right_nodes,
                           left_order, right_order, prefix) {
  e <- links %>% mutate(
    edge_id = paste0(prefix, seq_len(n())),
    from_rank = match(from, left_order),
    to_rank = match(to, right_order)
  )
  e <- e %>% arrange(from_rank, to_rank) %>% group_by(from) %>%
    mutate(left_offset = cumsum(value) - value) %>% ungroup()
  e <- e %>% arrange(to_rank, from_rank) %>% group_by(to) %>%
    mutate(right_offset = cumsum(value) - value) %>% ungroup()
  
  e <- e %>%
    left_join(left_nodes %>% select(name, y0) %>%
                rename(from = name, left_node_y0 = y0), by = "from") %>%
    left_join(right_nodes %>% select(name, y0) %>%
                rename(to = name, right_node_y0 = y0), by = "to") %>%
    mutate(
      y_left = left_node_y0 + left_offset,
      y_right = right_node_y0 + right_offset
    )
  stopifnot(!anyNA(e$y_left), !anyNA(e$y_right))
  e
}

ribbon_polygons <- function(edges, x_left, x_right, nsteps = 45) {
  t <- seq(0, 1, length.out = nsteps)
  s <- 3 * t^2 - 2 * t^3
  xv <- x_left + (x_right - x_left) * t
  bind_rows(lapply(seq_len(nrow(edges)), function(i) {
    r <- edges[i, ]
    low <- r$y_left + (r$y_right - r$y_left) * s
    high <- low + r$value
    data.frame(
      x = c(xv, rev(xv)),
      y = c(low, rev(high)),
      edge_id = r$edge_id,
      fill_group = r$group,
      stringsAsFactors = FALSE
    )
  }))
}

pretty_count <- function(x) {
  ifelse(abs(x - round(x)) < 1e-7,
         format(round(x), big.mark = ",", trim = TRUE, scientific = FALSE),
         format(round(x, 1), big.mark = ",", trim = TRUE, scientific = FALSE))
}

plot_partial_sankey <- function(edge12, edge23, title, subtitle, n_total,
                                n_pos, zoom = FALSE) {
  # Node counts are explicitly based on all edge12 links for the first 2 axes.
  left_totals <- edge12 %>% group_by(name = from) %>%
    summarise(value = sum(value), .groups = "drop")
  mid_totals <- edge12 %>% group_by(name = to) %>%
    summarise(value = sum(value), .groups = "drop")
  right_totals <- edge23 %>% group_by(name = to) %>%
    summarise(value = sum(value), .groups = "drop")
  
  gap_left <- max(7, n_total * 0.026)
  gap_right <- max(9, n_total * if (zoom) 0.045 else 0.030)
  n1 <- node_layout(left_totals, family_order, 0, gap_left)
  n2 <- node_layout(mid_totals, anno_order, 1, gap_left)
  H <- max(n1$y1, n2$y1)
  n3 <- node_layout(right_totals, ccre_order, 2, gap_right,
                    center_at = H / 2)
  
  e12 <- position_links(edge12, n1, n2, family_order, anno_order, "A")
  e23 <- position_links(edge23, n2, n3, anno_order, ccre_order, "B")
  r12 <- ribbon_polygons(e12, 0.105, 0.895)
  r23 <- ribbon_polygons(e23, 1.105, 1.895)
  
  bars <- bind_rows(
    n1 %>% mutate(fill_group = paste0("FAM:", name), stage = "Family"),
    n2 %>% mutate(fill_group = "ANN", stage = "Annotation"),
    n3 %>% mutate(fill_group = paste0("CRE:", name), stage = "cCRE")
  ) %>% mutate(y_mid = (y0 + y1)/2)
  
  label1 <- n1 %>% mutate(label = paste0(name, "  (", pretty_count(value), ")"),
                          x_text = -0.15)
  label2 <- n2 %>% mutate(label = paste0(name, "  (", pretty_count(value), ")"),
                          x_text = 1)
  label3 <- n3 %>% mutate(
    label = paste0(ifelse(name == "DNase-H3K4me3", "DNase–H3K4me3", name),
                   "  (", pretty_count(value), ")"),
    x_text = 2.15
  )
  
  ggplot() +
    geom_polygon(
      data = r12,
      aes(x, y, group = edge_id, fill = fill_group),
      alpha = 0.64, color = NA
    ) +
    geom_polygon(
      data = r23,
      aes(x, y, group = edge_id, fill = fill_group),
      alpha = 0.68, color = NA
    ) +
    geom_rect(
      data = bars,
      aes(xmin = x - 0.105, xmax = x + 0.105,
          ymin = y0, ymax = y1, fill = fill_group),
      color = "#50565C", linewidth = 0.27
    ) +
    geom_text(
      data = label1,
      aes(x = x_text, y = (y0 + y1) / 2, label = label),
      hjust = 1, size = 3.2, color = "#252B32"
    ) +
    geom_label(
      data = label2,
      aes(x = x_text, y = (y0 + y1) / 2, label = label),
      size = 3, label.size = 0, label.padding = grid::unit(0.12, "lines"),
      fill = "white", color = "#252B32"
    ) +
    geom_text(
      data = label3,
      aes(x = x_text, y = (y0 + y1) / 2, label = label),
      hjust = 0, size = 3, color = "#252B32"
    ) +
    annotate("text", x = c(0, 1, 2), y = -0.08 * H,
             label = c("TE family", "Genomic annotation", "cCRE"),
             fontface = "bold", size = 4, color = "#252B32") +
    scale_fill_manual(values = fill_colors, guide = "none") +
    scale_y_reverse() +
    coord_cartesian(xlim = c(-0.55, 2.9),
                    ylim = c(-0.135 * H, H + 0.035 * H), clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      caption = paste0(
        sprintf("%s / %s TE loci (%.1f%%) overlap a cCRE. ",
                pretty_count(n_pos), pretty_count(n_total), 100 * n_pos/n_total),
        "Only cCRE-positive loci extend into the third layer; ",
        "flow widths are proportional to TE-locus weights."
      )
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15, color = "#252B32"),
      plot.subtitle = element_text(size = 11, color = "#50565C"),
      plot.caption = element_text(size = 9, color = "#626971", hjust = 0),
      plot.margin = margin(15, 110, 22, 125)
    )
}

p_full <- plot_partial_sankey(
  links_full, links_ccre,
  title = "TE family → Genomic annotation → cCRE",
  subtitle = "All TE loci in the first two layers; cCRE-positive loci only in the third",
  n_total = n_all, n_pos = n_positive, zoom = FALSE
)

# In the magnified companion panel, the first two layers are deliberately
# restricted to the positive subset. This panel has its OWN quantitative scale
# and must not be interpreted as the global TE distribution.
links_pos_12 <- merged %>%
  filter(!is.na(cCRE_5class)) %>%
  group_by(TE_family, Genomic_class) %>%
  summarise(value = sum(TE_weight), .groups = "drop") %>%
  transmute(from = TE_family, to = Genomic_class, value,
            group = paste0("FAM:", from))

p_zoom <- plot_partial_sankey(
  links_pos_12, links_ccre,
  title = "cCRE-positive TE loci (magnified view)",
  subtitle = "Magnified companion panel; all three layers represent the same positive subset",
  n_total = n_positive, n_pos = n_positive, zoom = TRUE
)

# Vector outputs for manuscript layout; PNG files are 400 dpi.
ggsave(file.path(out_dir, "Fig_A_all_TE_to_cCRE_subset.pdf"),
       p_full, width = 13.5, height = 8.5)
ggsave(file.path(out_dir, "Fig_A_all_TE_to_cCRE_subset.png"),
       p_full, width = 13.5, height = 8.5, dpi = 400)
ggsave(file.path(out_dir, "Fig_B_cCRE_positive_magnified.pdf"),
       p_zoom, width = 13.5, height = 8.5)
ggsave(file.path(out_dir, "Fig_B_cCRE_positive_magnified.png"),
       p_zoom, width = 13.5, height = 8.5, dpi = 400)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(out_dir, "Fig_A_all_TE_to_cCRE_subset.svg"),
         p_full, width = 13.5, height = 8.5, device = "svg")
  ggsave(file.path(out_dir, "Fig_B_cCRE_positive_magnified.svg"),
         p_zoom, width = 13.5, height = 8.5, device = "svg")
}


# ----- 8. Console QC --------------------------------------------------------
cat("\n===============================================\n")
cat("Unique TE loci (layers 1 and 2):", n_all, "\n")
cat("cCRE-positive loci (layer 3): ", n_positive, "\n")
cat("Positive fraction: ", sprintf("%.2f%%", 100*n_positive/n_all), "\n")
cat("Sum of family -> annotation links:", sum(links_full$value), "\n")
cat("Sum of annotation -> cCRE links: ", sum(links_ccre$value), "\n")
cat("Outputs:", out_dir, "\n")
cat("===============================================\n")



```








