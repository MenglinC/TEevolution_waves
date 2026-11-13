# A

## A.0 consensus sequence
### 共识序列是怎么得到的？
```
以6kb为全长的阈值，找到hg38人类基因组上全长的位点，在考虑链的信息下bedtools getfasta提取fasta序列，随后用mafft --auto --thread 8 L1PA6_6kb.fasta >align/align_L1PA6_6kb.fasta进行多序列比对；随后用consensusFretch.py找到每个家族的共识序列，最后去除gap区域（-），即为最终的共识序列；
对于小鼠L1MdTf-I-5end的5-UTR共识序列，则直接从Dfam中下载，序列号为DF000001864；
```

```python

import argparse
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

# IUPAC 混合碱基符号（支持 GAP）
IUPAC_CODES = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['-']): '-',  # GAP
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N',
}

# # 从比重生成 IUPAC 符号
# def get_iupac_from_frequencies(freqs):
    # # 按比重降序排列碱基
    # sorted_bases = sorted(freqs.items(), key=lambda x: x[1], reverse=True)
    # max_freq = sorted_bases[0][1]
    # selected_bases = {base for base, freq in sorted_bases if freq == max_freq}

    # # 如果 gap 存在但不是唯一最高频率，排除 gap
    # if '-' in selected_bases and len(selected_bases) > 1:
        # selected_bases.remove('-')

    # # 打印调试信息：查看当前位点的选中碱基和最大比重
    # print(f"DEBUG: Selected bases: {selected_bases}, Frequencies: {freqs}")
    
    # # 返回对应的 IUPAC 符号
    # iupac = IUPAC_CODES.get(frozenset(selected_bases), 'N')  # 默认返回 N 表示不确定
    # print(f"DEBUG: Returning IUPAC: {iupac}")
    # return iupac

# # 计算共识序列和每个位点比重
# def get_frequencies_and_consensus(alignment):
    # frequencies = []
    # consensus_sequence = []
    # for i in range(alignment.get_alignment_length()):
        # column = alignment[:, i]
        # counts = Counter(column)
        # total = sum(counts.values())
        # freqs = {base: round(count / total, 3) for base, count in counts.items()}
        # frequencies.append(freqs)

        # consensus_base = get_iupac_from_frequencies(freqs)
        # consensus_sequence.append(consensus_base)

        # # 调试输出每个位点的共识碱基
        # print(f"DEBUG: Position {i+1}, Frequencies: {freqs}, Consensus Base: {consensus_base}")
    
    # consensus_str = "".join(consensus_sequence)
    # print(f"DEBUG: Final consensus sequence: {consensus_str}")
    # return frequencies, consensus_str

# # 从比重生成 IUPAC 符号
# def get_iupac_from_frequencies(freqs):
    # sorted_bases = sorted(freqs.items(), key=lambda x: x[1], reverse=True)
    # max_freq = sorted_bases[0][1]
    # selected_bases = {base for base, freq in sorted_bases if freq == max_freq}

    # # 如果 GAP 存在但不是唯一最高频率，排除 GAP
    # if '-' in selected_bases and len(selected_bases) > 1:
        # selected_bases.remove('-')

    # # 如果只有一个碱基具有最高频率，直接返回
    # if len(selected_bases) == 1:
        # return next(iter(selected_bases))

    # # 返回对应的 IUPAC 符号
    # iupac = IUPAC_CODES.get(frozenset(selected_bases), 'N')
    # return iupac
# 从比重生成 IUPAC 符号
def get_iupac_from_frequencies(freqs):
    sorted_bases = sorted(freqs.items(), key=lambda x: x[1], reverse=True)
    max_freq = sorted_bases[0][1]
    selected_bases = {base.upper() for base, freq in sorted_bases if freq == max_freq}

    # 如果 GAP 存在但不是唯一最高频率，排除 GAP
    if '-' in selected_bases and len(selected_bases) > 1:
        selected_bases.remove('-')

    # 如果只有一个碱基具有最高频率，直接返回
    if len(selected_bases) == 1:
        return next(iter(selected_bases))

    # 返回对应的 IUPAC 符号
    iupac = IUPAC_CODES.get(frozenset(selected_bases), 'N')
    return iupac.upper()


# 核心函数：计算比重和共识序列
def get_frequencies_and_consensus(alignment):
    frequencies = []
    consensus_sequence = []
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        counts = Counter(column)
        total = sum(counts.values())
        freqs = {base: round(count / total, 3) for base, count in counts.items()}
        frequencies.append(freqs)

        consensus_base = get_iupac_from_frequencies(freqs)
        consensus_sequence.append(consensus_base)
    return frequencies, "".join(consensus_sequence)


# 将共识序列写入 FASTA 文件
def write_consensus_to_fasta(consensus_sequence, output_file):
    with open(output_file, 'w') as f:
        f.write(f">Consensus_Sequence\n{consensus_sequence}\n")
        print(f"DEBUG: Writing consensus sequence to FASTA file: {consensus_sequence}")

# 清理 FASTA 文件中的换行符，整理成连续序列
def clean_fasta_sequences(input_file, cleaned_file):
    records = []
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        current_id = None
        current_seq = []
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                # 使用正则去除标题行中的非标准字符（如双引号）
                line = line.replace('"', '').replace('>', '')
                if current_id:
                    records.append(SeqRecord(Seq("".join(current_seq)), id=current_id))
                current_id = line  # 更新基因ID
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:  # 添加最后一个序列
            records.append(SeqRecord(Seq("".join(current_seq)), id=current_id))

    # 写入清理后的文件
    with open(cleaned_file, 'w') as outfile:
        for record in records:
            outfile.write(f">{record.id}\n{record.seq}\n")

    return cleaned_file

# 将频率写入文本文件
def write_frequencies_to_txt(frequencies, output_file):
    with open(output_file, 'w') as f:
        for i, freq in enumerate(frequencies):
            f.write(f"Position {i+1}: {freq}\n")
            print(f"DEBUG: Writing frequencies for position {i+1}: {freq}")

# 主函数
def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="计算多序列比对的共识序列和每个位点的比重，并输出到文件")
    parser.add_argument(
        "input_file", 
        type=str, 
        help="输入的多序列比对文件（FASTA 格式）"
    )
    parser.add_argument(
        "output_prefix",
        type=str,
        help="输出文件的前缀"
    )
    args = parser.parse_args()

    # 清理输入文件
    cleaned_file = f"{args.output_prefix}_cleaned.fasta"
    cleaned_file = clean_fasta_sequences(args.input_file, cleaned_file)

    # 读取比对文件
    try:
        alignment = AlignIO.read(cleaned_file, "fasta")
    except FileNotFoundError:
        print(f"错误: 文件 '{args.input_file}' 不存在。")
        return
    except Exception as e:
        print(f"错误: 无法读取文件 '{args.input_file}'。请确保文件格式正确。")
        print(f"详细错误: {e}")
        return
    
    # 计算比重和共识序列
    frequencies, consensus_sequence = get_frequencies_and_consensus(alignment)
    
    # 输出文件名
    fasta_file = f"{args.output_prefix}_consensus.fasta"
    txt_file = f"{args.output_prefix}_frequencies.txt"
    
    # 写入文件
    write_consensus_to_fasta(consensus_sequence, fasta_file)
    write_frequencies_to_txt(frequencies, txt_file)
    
    print(f"共识序列已保存到 {fasta_file}")
    print(f"位点比重已保存到 {txt_file}")

# 调用主函数
if __name__ == "__main__":
    main()
```




## A.1 alignment

#input:L1PA1-17conseq.261001.fastq

```shell
#L-INS-i
mafft-linsi L1PA1-17conseq.261001.fastq > L1PA.LINSI.261001.fasta


# G-INS-i
mafft --globalpair --maxiterate 1000 \
    L1PA1-17conseq.261001.fastq \
    > L1PA.GINSI.261001.fasta

# E-INS-i
mafft --genafpair --maxiterate 1000 \
    L1PA1-17conseq.261001.fastq \
    > L1PA.EINSI.261001.fasta
```

## A.2 analyse the alignment

#file:figure4a.R
#input:L1PA.LINSI.261001.fasta

```r

library(Biostrings)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(zoo)


#每一步用到的数据太多了，然后一点点的来看来理解；
setwd("F://10-01")
aln <- readDNAStringSet("L1PA.LINSI.261001.fasta") #读取局部比对的结果

#length(aln)          # 应该是17
#width(aln)           # 每条alignment后的长度应该一致
#names(aln)

# 转换成 matrix：
# row = sequence
# column = alignment position
aln_mat <- do.call(
  rbind,
  strsplit(as.character(aln), "")
)

dim(aln_mat) #做格式转换
# ==============================
# 2. Column-level statistics
# ==============================

calc_column_stats <- function(x) {
  
  x <- toupper(x)
  
  # 1. gap比例：只统计 "-"
  gap_density <- mean(x == "-")
  
  # 2. 去掉gap，剩余字符都参与碱基组成复杂度计算
  bases <- x[x != "-"]
  
  # 如果这一列全是gap
  if (length(bases) == 0) {
    return(c(
      gap_density = gap_density,
      entropy = NA
    ))
  }
  
  # 各字符频率
  freq <- table(bases) / length(bases)
  
  # Shannon entropy
  entropy <- -sum(freq * log2(freq))
  
  c(
    gap_density = gap_density,
    entropy = entropy
  )
}
#函数calc_column_stats，计算熵和gap的比例


column_stats <- t(
  apply(aln_mat, 2, calc_column_stats)
) %>%
  as.data.frame()



column_stats$position <- seq_len(nrow(column_stats)) #生成的数据，丢掉了位置信息；

#head(column_stats)



############################################
# 计算比对稳定性
###########################################
library(Biostrings)
L <- readDNAStringSet("L1PA.LINSI.261001.fasta")
G <- readDNAStringSet("L1PA.GINSI.261001.fasta")
E <- readDNAStringSet("L1PA.EINSI.261001.fasta")

# 确保sequence names完全一致
stopifnot(
  setequal(names(L), names(G)),
  setequal(names(L), names(E))
)

# 按L-INS-i的顺序排列
G <- G[names(L)]
E <- E[names(L)]

make_residue_map <- function(aln) {
  
  mat <- do.call(
    rbind,
    strsplit(as.character(aln), "")
  )
  
  rownames(mat) <- names(aln)
  
  maps <- vector("list", nrow(mat))
  names(maps) <- rownames(mat)
  
  for (i in seq_len(nrow(mat))) {
    
    x <- mat[i, ]
    
    # 每一个alignment position对应原序列第几个base
    pos <- cumsum(x != "-")
    pos[x == "-"] <- NA
    
    maps[[i]] <- pos
  }
  
  list(
    matrix = mat,
    map = maps
  )
} #将比对后的位置与序列原来的位置进行映射


map_L <- make_residue_map(L)
map_G <- make_residue_map(G)
map_E <- make_residue_map(E)

get_aligned_pairs <- function(aln_map) {
  
  mat <- aln_map$matrix
  maps <- aln_map$map
  
  pair_set <- new.env(hash = TRUE)
  
  nseq <- nrow(mat)
  ncol_aln <- ncol(mat)
  
  for (j in seq_len(ncol_aln)) {
    
    present <- which(mat[, j] != "-")
    
    if (length(present) < 2)
      next
    
    cmb <- combn(present, 2)
    
    for (k in seq_len(ncol(cmb))) {
      
      i1 <- cmb[1, k]
      i2 <- cmb[2, k]
      
      seq1 <- rownames(mat)[i1]
      seq2 <- rownames(mat)[i2]
      
      pos1 <- maps[[seq1]][j]
      pos2 <- maps[[seq2]][j]
      
      key <- paste(
        seq1, pos1,
        seq2, pos2,
        sep = "|"
      )
      
      pair_set[[key]] <- TRUE
    }
  }
  
  pair_set
}

pairs_G <- get_aligned_pairs(map_G)
pairs_E <- get_aligned_pairs(map_E)

calc_stability <- function(ref_map, alt_pair_sets) {
  
  mat <- ref_map$matrix
  maps <- ref_map$map
  
  stability <- rep(NA_real_, ncol(mat))
  
  for (j in seq_len(ncol(mat))) {
    
    present <- which(mat[, j] != "-")
    
    if (length(present) < 2)
      next
    
    cmb <- combn(present, 2)
    
    scores <- c()
    
    for (k in seq_len(ncol(cmb))) {
      
      i1 <- cmb[1, k]
      i2 <- cmb[2, k]
      
      seq1 <- rownames(mat)[i1]
      seq2 <- rownames(mat)[i2]
      
      pos1 <- maps[[seq1]][j]
      pos2 <- maps[[seq2]][j]
      
      key <- paste(
        seq1, pos1,
        seq2, pos2,
        sep = "|"
      )
      
      for (pair_set in alt_pair_sets) {
        
        scores <- c(
          scores,
          exists(
            key,
            envir = pair_set,
            inherits = FALSE
          )
        )
      }
    }
    
    stability[j] <- mean(scores)
  }
  
  stability
}

column_stats$stability <- calc_stability(
  map_L,
  list(pairs_G, pairs_E)
)

# ==========================================
window_size=10
column_stats <- column_stats %>%
  mutate(
    gap_density_50bp = zoo::rollapply(
      gap_density,
      width = window_size,
      FUN = mean,
      partial = TRUE,
      align = "center"
    ),
    
    entropy_50bp = zoo::rollapply(
      entropy,
      width = window_size,
      FUN = mean,
      partial = TRUE,
      align = "center"
    ),
    stability_50bp = zoo::rollapply(
      stability,
      width = window_size,
      FUN = mean,
      partial = TRUE,
      align = "center"
    )
  )



# ==========================================
# 上面就完成了计算，下面就是画图环节
# ==========================================

theme_track <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 8)
    ),
    
    axis.text = element_text(size = 9),
    
    plot.margin = margin(
      t = 3,
      r = 8,
      b = 3,
      l = 8
    )
  )

plot_dat<-column_stats
# ==========================================
# A. gap density
# ==========================================
p_gap <- ggplot(
  plot_dat,
  aes(x = position)
) +
  
  # raw per-site values
  geom_line(
    aes(y = gap_density),
    linewidth = 0.25,
    alpha = 0.25
  ) +
  
  # 50-bp rolling mean
  geom_line(
    aes(y = gap_density_50bp),
    linewidth = 0.2
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = "Gap density"
  ) +
  
  theme_track +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_gap
# ==========================================
# B. Alignment stability
# ==========================================

p_stability <- ggplot(
  plot_dat,
  aes(x = position)
) +
  
  geom_line(
    aes(y = stability),
    linewidth = 0.25,
    alpha = 0.25
  ) +
  
  geom_line(
    aes(y = stability_50bp),
    linewidth = 0.2
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = "Alignment\nstability"
  ) +
  
  theme_track +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_stability
# ==========================================
# C. Sequence variability
# ==========================================

p_var <- ggplot(
  plot_dat,
  aes(x = position)
) +
  
  geom_line(
    aes(y = entropy),
    linewidth = 0.25,
    alpha = 0.25
  ) +
  
  geom_line(
    aes(y = entropy_50bp),
    linewidth = 0.2
  ) +
  
  scale_y_continuous(
    limits = c(0, 0.8),
    breaks = c(0, 0.4, 0.8),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Alignment position (bp)",
    y = "Sequence\nvariability"
  ) +
  
  theme_track

p_var

####################################################
#绘制L1的位置模型图
####################################################

ref_name <- grep(
  "L1HS",
  rownames(aln_mat),
  value = TRUE
)[1]

ref_name

# 提取L1PA1在alignment中的序列
ref <- toupper(aln_mat[ref_name, ])

# 哪些alignment columns在L1PA1中是真实碱基
#ref_is_base <- ref %in% c("A", "C", "G", "T")
ref_is_base <- ref != "-" #除了gap之外就是有效碱基

# 从5'到3'累计编号
ref_coordinate <- cumsum(ref_is_base)

# L1PA1本身是gap的位置没有reference coordinate
ref_coordinate[!ref_is_base] <- NA

# 检查
length(ref_coordinate)
nrow(column_stats)

# 两者必须相同
stopifnot(length(ref_coordinate) == nrow(column_stats))

# 加入column_stats
column_stats$ref_coordinate <- ref_coordinate



ref_length <- max(
  column_stats$ref_coordinate,
  na.rm = TRUE
)

annotation_ref <- data.frame(
  feature = c(
    "5'UTR",
    "ORF1",
    "Inter-ORF",
    "ORF2",
    "3'UTR"
  ),
  # start = c(
  #   1,
  #   2326,
  #   3493,
  #   3556,
  #   7396
  # ),
  # end = c(
  #   2325,
  #   3492,
  #   3555,
  #   min(7395, ref_length),
  #   ref_length
  # )
  # start = c(
  #   1,
  #   2360,
  #   3537,
  #   3600,
  #   7436
  # ),
  # end = c(
  #   2359,
  #   3536,
  #   3599,
  #   min(7435, ref_length),
  #   ref_length
  # )
  start = c(
    1,
    905,
    1923,
    1986,
    5814
  ),
  end = c(
    906,
    1922,
    1985,
    min(5813, ref_length),
    ref_length
  )
)

annotation_ref <- annotation_ref %>%
  filter(start <= end)

annotation_ref

get_alignment_range <- function(start, end, column_stats) {
  
  idx <- which(
    !is.na(column_stats$ref_coordinate) &
      column_stats$ref_coordinate >= start &
      column_stats$ref_coordinate <= end
  )
  
  if (length(idx) == 0) {
    return(
      c(
        xmin = NA,
        xmax = NA
      )
    )
  }
  
  c(
    xmin = min(column_stats$position[idx]),
    xmax = max(column_stats$position[idx])
  )
}

ranges <- lapply(
  seq_len(nrow(annotation_ref)),
  function(i) {
    
    get_alignment_range(
      annotation_ref$start[i],
      annotation_ref$end[i],
      column_stats
    )
  }
)

ranges <- do.call(rbind, ranges)

orf_dat <- cbind(
  annotation_ref,
  ranges
)

orf_dat

p_structure <- ggplot() +
  
  geom_segment(
    aes(
      x = min(column_stats$position),
      xend = max(column_stats$position),
      y = 0.5,
      yend = 0.5
    ),
    linewidth = 0.5
  ) +
  
  geom_rect(
    data = orf_dat %>%
      filter(feature %in% c("5'UTR", "3'UTR")),
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = 0.38,
      ymax = 0.62
    ),
    fill = "grey90",
    color = "black",
    linewidth = 0.4
  ) +
  
  geom_rect(
    data = orf_dat %>%
      filter(feature %in% c("ORF1", "ORF2")),
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = 0.25,
      ymax = 0.75
    ),
    fill = "grey65",
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_text(
    data = orf_dat %>%
      filter(feature != "Inter-ORF"),
    aes(
      x = (xmin + xmax) / 2,
      y = 0.5,
      label = feature
    ),
    size = 3.3
  ) +
  
  scale_x_continuous(
    limits = range(column_stats$position),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = NULL
  ) +
  
  theme_classic() +
  
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

p_structure

#################################################
#同一坐标
#################################################

x_scale <- scale_x_continuous(
  limits = range(column_stats$position),
  expand = c(0, 0)
)

p_structure <- p_structure + x_scale
p_gap       <- p_gap + x_scale
p_stability <- p_stability + x_scale
p_var       <- p_var + x_scale

################################################

library(patchwork)
p_final <- (
  p_structure /
    p_gap /
    p_stability /
    p_var
) +
  plot_layout(
    heights = c(
      0.4,
      1,
      1,
      1
    )
  )

p_final

```
