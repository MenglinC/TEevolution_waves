# C
## C.1.1

```shell
conda activate MAFFT
mafft-linsi L1_group5_consensus.fasta > L1_group5_linsi.fasta

```

## C.1.2
#inner group identical region
```r
#workplace:/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/motif_evolution_fig4_0810/L1_subfamily_conseq/261001_4.1/
#workplace:/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/motif_evolution_fig4_0810/L1_subfamily_conseq/261001_4.1/group/flmo/

library(Biostrings)
library(dplyr)

# ==============================
# Parameters
# ==============================
setwd("F://10-01//group")
input_fasta <- "L1_group5_linsi.fasta"
output_fasta <- "L1_group5_identical_regions.fasta"
output_table <- "L1_group5_identical_regions.tsv"

min_length <- 10


# ==============================
# Read alignment
# ==============================

aln <- readDNAStringSet(input_fasta)

aln_mat <- do.call(
  rbind,
  strsplit(toupper(as.character(aln)), "")
)

rownames(aln_mat) <- names(aln)


# ==============================
# Identify completely identical columns
# ==============================

identical_col <- apply(
  aln_mat,
  2,
  function(x) {
    
    length(unique(x)) == 1 &&
      x[1] %in% c("A", "T", "C", "G")
    
  }
)


# ==============================
# Merge consecutive columns
# ==============================

r <- rle(identical_col)

end_pos <- cumsum(r$lengths)
start_pos <- end_pos - r$lengths + 1

regions <- tibble(
  start = start_pos,
  end = end_pos,
  identical = r$values
) %>%
  filter(identical) %>%
  mutate(
    length = end - start + 1
  ) %>%
  filter(length >= min_length)


# ==============================
# Extract sequence
# ==============================

regions <- regions %>%
  rowwise() %>%
  mutate(
    sequence = paste0(
      aln_mat[1, start:end],
      collapse = ""
    )
  ) %>%
  ungroup() %>%
  mutate(
    region_id = paste0(
      "identical_region_",
      row_number(),
      "_aln",
      start,
      "-",
      end,
      "_len",
      length
    )
  )


# ==============================
# Write FASTA
# ==============================

out_fasta <- DNAStringSet(regions$sequence)

names(out_fasta) <- regions$region_id

writeXStringSet(
  out_fasta,
  filepath = output_fasta
)


# ==============================
# Write region table
# ==============================

write.table(
  regions,
  file = output_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


```

## C.1.3

```shell

fimo -o fimo_group9  /home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/motif_evolution_fig4_0810/batch/motif_db/jaspar_9606.meme L1PA17_consensus.fasta
```


## C.1.4

```r
# 目录
setwd("F://10-01//summary_fimo_group")

# 找到所有 txt 文件
files <- list.files(
  pattern = "\\.txt$",
  full.names = TRUE
)

length(files)   # 应该是 9

data_list <- lapply(files, function(f) {
  
  dat <- read.table(
    f,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
  )
  
  unique(dat[[1]])
})

# 用文件名作为 list 名称
#names(data_list) <- basename(files)
names(data_list) <- c("group1","group2","group3","group5","group4","group6","group7","group8","group9")



intersection <- Reduce(
  intersect,
  data_list
)

length(intersection)
intersection

all_elements <- unique(unlist(data_list))

presence <- sapply(data_list, function(x) {
  all_elements %in% x
})

presence_df <- data.frame(
  element = all_elements,
  presence,
  n_files = rowSums(presence),
  check.names = FALSE
)

presence_df <- presence_df[
  order(presence_df$n_files, decreasing = TRUE),
]

head(presence_df)
write.csv(presence_df,"TF_summary.group.csv")

##################################
#画图
##################################
library(tidyverse)

# 1. 9个group的顺序
#group_order <- names(data_list)
group_order <- c("group1","group2","group3","group4","group5",
                 "group6","group7","group8","group9")


# 2. 所有TF
all_TF <- unique(unlist(data_list))

# 3. 构建长表
presence_long <- expand_grid(
  TF = all_TF,
  group = group_order
) %>%
  mutate(
    present = map2_lgl(TF, group, ~ .x %in% data_list[[.y]]),
    present = as.integer(present),
    group_idx = match(group, group_order)
  )
presence_long <- presence_long %>%
  mutate(
    TF = sub("_", "|", TF)
  )

tf_stat <- presence_long %>%
  group_by(TF) %>%
  summarise(
    first_group = min(group_idx[present == 1]),
    last_group  = max(group_idx[present == 1]),
    n_group     = sum(present),
    .groups = "drop"
  ) %>%
  filter(n_group >= 6) %>%   # 如果想改成“恰好5个group”，这里改成 == 5
  arrange(
    first_group,
    desc(n_group),
    last_group,
    TF
  )

tf_order <- tf_stat$TF

plot_dat <- presence_long %>%
  filter(TF %in% tf_order) %>%
  mutate(
    group = factor(group, levels = group_order),
    TF = factor(TF, levels = rev(tf_order))
  )

library(ggplot2)

p_tf <- ggplot(
  plot_dat,
  aes(x = group, y = TF)
) +
  geom_point(
    aes(
      fill = factor(present),
      color = factor(present)
    ),
    shape = 21,
    size = 3.2,
    stroke = 0.4
  ) +
  scale_fill_manual(
    values = c(
      "0" = "#dfdfdf",
      "1" = "#559a61"
    )
  ) +
  scale_color_manual(
    values = c(
      "0" = "#dfdfdf",
      "1" = "#559a61"
    )
  ) +
  labs(
    x = "",
    y = ""
  ) +
  guides(
    fill = "none",
    color = "none"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    )
  )

p_tf

pdf(paste("TF_summary_group.pdf"),width =3,height =6)
print(p_tf)
dev.off()
```
