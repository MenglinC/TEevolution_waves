# D.gene cluster&TE density

## D.1 prepare input files
```shell
#workplace:/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/gene_family_fig1_1108/
#SD files download from UCSC tablebrowser
#CNV files download from ClinGen_recurrent_CNV_V2.0 hg38
#gene family download from HGNC
#gene download from gencode
```


## D.2 do the enrichment

```shell

pip install  statsmodels

python window_500kb_enrichment.py \
    --te TE_500kb.tsv \
    --genes genes.bed \
    --families gene_family.tsv \
    --sd SD.bed \
    --cnv CNV.bed \
    --top 0.05 \
    --sd-cutoff 0.10 \
    --cnv-cutoff 0.10 \
    --prefix Fig1

```

## D.3 plot

```r
# Fig1 TE hotspot enrichment bubble plot — white background (no grey tiles)
# Input: Fig1_enrichment.tsv
# Required columns: feature, feature_type, TE_group, FDR, log2_OR_plot
# Packages: install.packages(c('dplyr','ggplot2','scales'))

library(dplyr)
library(ggplot2)
library(scales)

input_file <- 'Fig1_enrichment.tsv'
#output_file <- 'Fig1_TE_enrichment_OR_FDR_white.png'
df <- read.delim(input_file, sep = '\t', check.names = FALSE,
                 stringsAsFactors = FALSE)
stopifnot(all(c('feature', 'feature_type', 'TE_group',
                'FDR', 'log2_OR_plot') %in% names(df)))

# 1) 与上一版示意图一致：gene_family_multi_2plus 中取最常显著的10项，
#    再补充两类 olfactory receptors；SD/CNV 单列于底部。
#    如果想人工选择家族，直接替换 selected_families 即可。
top_families <- df %>%
  filter(feature_type == 'gene_family_multi_2plus') %>%
  group_by(feature) %>%
  summarise(n_sig = sum(FDR < 0.05, na.rm = TRUE),
            min_FDR = min(FDR, na.rm = TRUE), .groups = 'drop') %>%
  filter(n_sig >= 2) %>%
  arrange(desc(n_sig), min_FDR) %>%
  slice_head(n = 10) %>%
  pull(feature)

selected_families <- unique(c(
  top_families,
  'Olfactory receptors, family 10',
  'Olfactory receptors, family 4'
))

region_features <- c('SD_coverage', 'recurrent_CNV_region')
region_labels <- c('Segmental duplication', 'Recurrent CNV')

te_order <- c(
  'non-Primates_L1', 'Primates_L1', 'Simiiformes_L1',
  'Catarrhini_L1', 'Hominoidea_L1',
  'Primates_Alu', 'Simiiformes_Alu', 'Catarrhini_Alu', 'Hominoidea_Alu',
  'non-Primates_ERV1', 'Primates_ERV1', 'Simiiformes_ERV1',
  'Catarrhini_ERV1', 'Hominoidea_ERV1',
  'Simiiformes_ERVK', 'Catarrhini_ERVK', 'Hominoidea_ERVK',
  'Hominoidea_SVA'
)

# feature 自下向上依次为：CNV、SD、其他家族（图中顶部为ZNF）
feature_order <- rev(c(selected_families, region_labels))

plot_df <- df %>%
  filter((feature_type == 'gene_family_multi_2plus' &
            feature %in% selected_families) |
           (feature_type == 'genomic_region' &
              feature %in% region_features)) %>%
  mutate(
    feature = recode(feature,
                     SD_coverage = 'Segmental duplication',
                     recurrent_CNV_region = 'Recurrent CNV'),
    feature = factor(feature, levels = feature_order),
    TE_group = factor(TE_group, levels = te_order),
    fdr_score = pmin(-log10(pmax(FDR, .Machine$double.xmin)), 20)
  ) %>%
  # 只显示FDR<0.05的点；其他格子完全留白
  filter(!is.na(FDR), FDR < 0.05)

stopifnot(!anyDuplicated(plot_df[, c('TE_group', 'feature')]))

p <- ggplot(plot_df, aes(x = TE_group, y = feature)) +
  # 可选的家族分隔线（没有灰色方格背景）
  geom_vline(xintercept = c(5.5, 9.5, 14.5, 17.5),
             colour = 'grey80', linewidth = 0.4) +
  geom_hline(yintercept = 2.5, colour = 'grey65', linewidth = 0.5) +
  geom_point(aes(fill = log2_OR_plot, size = fdr_score),
             shape = 21, stroke = 0) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis_c(
    option = "D",
    direction = 1,
    name = expression(log[2](OR))
  )+
  # scale_fill_gradientn(
  #   colours = c('#fff7bc', '#fec44f', '#fe5f2c', '#e41a1c', '#7f002b'),
  #   limits = c(0, 3.5), oob = scales::squish,
  #   breaks = c(0, 1, 2, 3, 3.5),
  #   name = expression(log[2](OR)~'(corrected)')
  # ) +
  scale_size_continuous(
    limits = c(0, 20), #oob = scales::squish,
    range = c(1, 8), breaks = c(2, 5, 10, 20),
    name = expression(-log[10](FDR))
  ) +
  labs(
   # title = 'TE hotspot enrichment across gene families and genomic regions',
    #subtitle = 'Gene family: >=2 genes per 500-kb window · dots indicate FDR < 0.05',
    x = NULL, y = NULL
  ) +
  guides(fill = guide_colorbar(order = 1),
         size = guide_legend(order = 2,
                             override.aes = list(fill = 'black'))) +
  theme_classic(base_size = 12) +
  theme(
    # 关键：无灰底、无tile、无网格
    plot.background = element_rect(fill = 'white', colour = NA),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 16, face = 'bold', margin = margin(b = 7)),
    plot.subtitle = element_text(size = 11, colour = 'grey35',
                                 margin = margin(b = 15)),
    legend.position = 'right',
    plot.margin = margin(12, 14, 12, 12)
  )

print(p)

pdf("dotplot.pdf",width =10,height = 6)
print(p)
dev.off()

```

