# A

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
