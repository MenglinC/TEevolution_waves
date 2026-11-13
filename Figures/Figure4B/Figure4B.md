# B

## B.1

#input:L1PA.LINSI.261001.5utr.fasta

```shell
conda install -c bioconda emboss

water \
  -asequence L1_consensus.LINSI-5utr.fasta\
  -bsequence L1_consensus.LINSI-5utr.fasta\
  -gapopen 10 \
  -gapextend 0.5 \
  -outfile L1_5UTR_all_vs_all.water
```


## B.2

#file:figure4b.R

```r
library(dplyr)
library(stringr)
library(readr)

parse_emboss <- function(file, method = "global") {
  
  x <- readLines(file, warn = FALSE)
  
  res <- list()
  cur <- list()
  k <- 0
  
  for (line in x) {
    
    # sequence 1
    if (str_detect(line, "^# 1:")) {
      cur$seq1 <- str_trim(str_remove(line, "^# 1:"))
    }
    
    # sequence 2
    if (str_detect(line, "^# 2:")) {
      cur$seq2 <- str_trim(str_remove(line, "^# 2:"))
    }
    
    # alignment length
    if (str_detect(line, "^# Length:")) {
      cur$length <- as.numeric(
        str_match(line, "^# Length:\\s+(\\d+)")[,2]
      )
    }
    
    # Identity
    if (str_detect(line, "^# Identity:")) {
      
      m <- str_match(
        line,
        "^# Identity:\\s+(\\d+)/(\\d+)\\s+\\(\\s*([0-9.]+)%\\)"
      )
      
      cur$identity_n   <- as.numeric(m[2])
      cur$identity_den <- as.numeric(m[3])
      cur$identity_pct <- as.numeric(m[4])
    }
    
    # Similarity
    if (str_detect(line, "^# Similarity:")) {
      
      m <- str_match(
        line,
        "^# Similarity:\\s+(\\d+)/(\\d+)\\s+\\(\\s*([0-9.]+)%\\)"
      )
      
      cur$similarity_n   <- as.numeric(m[2])
      cur$similarity_den <- as.numeric(m[3])
      cur$similarity_pct <- as.numeric(m[4])
    }
    
    # Gaps
    if (str_detect(line, "^# Gaps:")) {
      
      m <- str_match(
        line,
        "^# Gaps:\\s+(\\d+)/(\\d+)\\s+\\(\\s*([0-9.]+)%\\)"
      )
      
      cur$gaps_n   <- as.numeric(m[2])
      cur$gaps_den <- as.numeric(m[3])
      cur$gaps_pct <- as.numeric(m[4])
    }
    
    # Score = 一个alignment结束
    if (str_detect(line, "^# Score:")) {
      
      cur$score <- as.numeric(
        str_match(line, "^# Score:\\s+([-0-9.]+)")[,2]
      )
      
      cur$method <- method
      
      k <- k + 1
      res[[k]] <- as.data.frame(cur)
      
      cur <- list()
    }
  }
  
  bind_rows(res)
}


needle <- parse_emboss(
  "L1_5UTR_all_vs_all.261001.pair",
  method = "global"
)



L1_order <- c(
  "L1HS",
  "L1PA2",
  "L1PA3",
  "L1PA4",
  "L1PA5",
  "L1PA6",
  "L1PA7",
  "L1PA8",
  "L1PA8A",
  "L1PA10",
  "L1PA11",
  "L1PA12",
  "L1PA13",
  "L1PA14",
  "L1PA15",
  "L1PA16",
  "L1PA17"
)

needle <- needle %>%
  mutate(
    seq1 = factor(seq1, levels = rev(L1_order)),
    seq2 = factor(seq2, levels = L1_order)
  )
library(ggplot2)

ggplot(
  needle,
  aes(x = seq2, y = seq1, fill = identity_pct)
) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Identity (%)"
  ) +
  coord_equal() +
  labs(
    x = "",
    y = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  )

################
# group
################
library(dplyr)
library(tidyr)

# 先保证 seq1/seq2 是字符
needle2 <- needle %>%
  mutate(
    seq1 = as.character(seq1),
    seq2 = as.character(seq2)
  )

# 构建 identity matrix
id_mat <- matrix(
  100,
  nrow = length(L1_order),
  ncol = length(L1_order),
  dimnames = list(L1_order, L1_order)
)

for (i in seq_len(nrow(needle2))) {
  s1 <- needle2$seq1[i]
  s2 <- needle2$seq2[i]
  id <- needle2$identity_pct[i]
  
  id_mat[s1, s2] <- id
  id_mat[s2, s1] <- id
}

diag(id_mat) <- 100

dist_mat <- as.dist(100 - id_mat)

hc <- hclust(
  dist_mat,
  method = "complete"
)

groups <- cutree(
  hc,
  h = 30 #调整阈值
)

group_df <- tibble(
  L1 = names(groups),
  group = groups
)

group_df
#########################



```
