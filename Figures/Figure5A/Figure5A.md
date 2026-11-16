# A

```r
setwd("F://11-15-画图//5A")
library(tidyverse)
# 目录
setwd("F://10-01//summary_fimo_group")

files <- list.files(
  path    = "F:/tmp/10-01/summary_fimo_group",
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

#从tf_stat这里开始筛选的
tf_stat <- presence_long %>%
  group_by(TF) %>%
  summarise(
    first_group = min(group_idx[present == 1]),
    last_group  = max(group_idx[present == 1]),
    n_group     = sum(present),
    .groups = "drop"
  ) %>%
  filter(n_group <= 1) %>%   # 如果想改成“恰好5个group”，这里改成 == 5 #筛选特异的
  arrange(
    first_group,
    desc(n_group),
    last_group,
    TF
  )

# tf_order <- tf_stat$TF
# 
# plot_dat <- presence_long %>%
#   filter(TF %in% tf_order) %>%
#   mutate(
#     group = factor(group, levels = group_order),
#     TF = factor(TF, levels = rev(tf_order))
#   )

tf_pick <- c(
  "MA0047.4|FOXA2",
  "MA0684.1|RUNX3",
#  "MA0007.2|AR",
  "MA0486.1|HSF1",
 # "MA1520.1|MAF",
  "MA1523.2|MSANTD3", #3
  "MA0595.1|SREBF1",#4
  "MA0105.2|NFKB1",#5
  "MA1116.1|RBPJ", #6
  "MA1634.1|BATF",#7
  "MA1541.1|NR6A1",#8
  "MA0660.1|MEF2B"#9
)

# 只保留这些 TF，并按 tf_pick 的顺序设因子水平
plot_dat <- presence_long %>%
  filter(TF %in% tf_pick) %>%
  mutate(
    group = factor(group, levels = group_order),
    TF    = factor(TF, levels = rev(tf_pick))   # rev 使第一个在图的顶部/底部按需
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

pdf(paste("TF_summary_group-specific.pdf"),width =3,height =4)
print(p_tf)
dev.off()



```
