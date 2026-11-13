# D
## D.1.1

```shell
blastn \
  -query L1MdA_I.fa \
  -subject L1HS.fa \
  -word_size 7 \
  -dust no \
  -evalue 1000 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  > L1MdA_I_vs_L1HS.blast.tsv

```


## D.1.2
#input:L1MdA_I_vs_L1HS.blast.tsv

```r
library(ggplot2)
library(dplyr)

blast <- read.table(
  "L1MdA_I_vs_L1HS.blast.tsv",
  sep = "\t",
  header = FALSE
)

colnames(blast) <- c(
  "query", "subject", "identity", "length",
  "mismatch", "gapopen",
  "qstart", "qend",
  "sstart", "send",
  "evalue", "bitscore"
)
summary(blast$length)
summary(blast$identity)

blast2 <- blast %>%
  filter(
    length >= 8,
    identity >= 65
  ) %>%
  arrange(qstart)

library(ggplot2)
library(dplyr)

# 两条 UTR 的长度
len_mouse <- 1328
len_human <- 905

# 每个 alignment block 构造成 polygon
ribbon_df <- bind_rows(
  lapply(seq_len(nrow(blast2)), function(i) {
    
    x <- blast2[i, ]
    
    data.frame(
      block = i,
      x = c(
        x$qstart,
        x$qend,
        x$send,
        x$sstart
      ),
      y = c(
        1,
        1,
        0,
        0
      ),
      identity = x$identity
    )
  })
)

p <- ggplot() +
  
  # conserved blocks
  geom_polygon(
    data = ribbon_df,
    aes(
      x = x,
      y = y,
      group = block,
      fill = identity
    ),
    alpha = 0.65,
    color = NA
  ) +
  
  # mouse 5'UTR
  geom_segment(
    aes(
      x = 1,
      xend = len_human,
      y = 1,
      yend = 1
    ),
    linewidth = 4,
    color = "grey30"
  ) +
  
  # human 5'UTR
  geom_segment(
    aes(
      x = 1,
      xend = len_mouse,
      y = 0,
      yend = 0
    ),
    linewidth = 4,
    color = "grey30"
  ) +
  
  scale_fill_gradient(
    low = "#D9EAF7",
    high = "#2166AC",
    name = "Identity (%)"
  ) +
  
  scale_y_continuous(
    breaks = c(0, 1),
    labels = c(
      "L1MdA_I nMouse",
      "L1HS Human"
    ),
    limits = c(-0.15, 1.15)
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  
  labs(
    x = "Position in 5′UTR (bp)",
    y = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "right"
  )

p

pdf("5utr_align.pdf",width =5,height =1.5)
p
dev.off()
```
