# S3
## S3.1

```shell
#workplace:/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/chr_enrichment/

cp "/home/xxzhang/Resource/Reference/Genome/homo_sapiens/T2T/hs1.genome" ./
awk 'BEGIN{OFS="\t"}{print $1,0,$2}' hs1.genome > hs1.genome.bed
bedtools intersect -a hs1.genome.bed -b "/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/clade_family_bp/non-Primates_L1.bed" -c > chr_counts_non-Primates.bed
bedtools intersect -a hs1.genome.bed -b "/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/clade_family_bp/Primates_L1.bed" -c > chr_counts_Primates.bed
bedtools intersect -a hs1.genome.bed -b "/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/clade_family_bp/Catarrhini_L1.bed" -c > chr_counts_Catarrhini.bed
bedtools intersect -a hs1.genome.bed -b "/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/clade_family_bp/Simiiformes_L1.bed" -c >chr_counts_Simiiformes.bed
bedtools intersect -a hs1.genome.bed -b "/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/clade_family_bp/Hominoidea_L1.bed" -c >chr_counts_Hominoidea.bed
```


## S3.2

```r
#file:http://localhost:8888/notebooks/chr_enrichment_hs1.ipynb

library(dplyr)
# 设置目录
base_dir <- "/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/chr_enrichment/"

# 获取所有以 chr_count 开头的文件
files <- list.files(
  path = base_dir,
  pattern = "^chr_count",
  full.names = TRUE
)

# 批量处理
lapply(files, function(file) {

  # 读取数据
  chr_data <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "start", "end", "observed")
  )

  # 计算染色体富集
  chr_result <- chr_data %>%
    mutate(
      observed = as.numeric(observed),
      start = as.numeric(start),
      end = as.numeric(end),
      chr_length = end - start
    ) %>%
    mutate(
      genome_fraction = chr_length / sum(chr_length, na.rm = TRUE),
      expected = sum(observed, na.rm = TRUE) * genome_fraction,
      density_per_Mb = observed / chr_length * 1e6,
      enrichment = observed / expected
    )

  # 生成输出文件名
  output_file <- file.path(
    base_dir,
    sub("^chr_count", "chr_enrich", basename(file))
  )

  # 输出结果
  write.table(
    chr_result,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("Finished:", basename(file), "\n")
})

cat("All files processed!\n")

```

## S3.3

```r
setwd("/home/xxzhang/Resource/Reference/GTF/TE_gtf/T2T/chr_enrichment/")

dat1<-read.table("chr_enrichs_non-Primates_L1.tsv",header=T)[,c(1,9)]
dat2<-read.table("chr_enrichs_Primates_L1.tsv",header=T)[,c(1,9)]
dat3<-read.table("chr_enrichs_Simiiformes_L1.tsv",header=T)[,c(1,9)]
dat4<-read.table("chr_enrichs_Catarrhini_L1.tsv",header=T)[,c(1,9)]
dat5<-read.table("chr_enrichs_Hominoidea_L1.tsv",header=T)[,c(1,9)]

dat12<-merge(dat1,dat2,by="chr")
dat123<-merge(dat12,dat3,by="chr")
dat1234<-merge(dat123,dat4,by="chr")
dat12345<-merge(dat1234,dat5,by="chr")

colnames(dat12345)<-c("chr","non-Primates","Primates","Simiiformes","Catarrhini","Hominoidea")
head(dat12345)

library(tidyr)
dat_long <- gather(dat12345, key = "stage", value = "observed_expected",
                    -`chr`)

dat_long$chr<-factor(dat_long$chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                           "chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                                           "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"))

dat_long$stage<-factor(dat_long$stage,levels=c("non-Primates",
    "Primates",
    "Simiiformes",
    "Catarrhini",
    "Hominoidea"))

options(repr.plot.width =12 , repr.plot.height = 4)
library(ggplot2)
p1<-ggplot(dat_long, aes(x = chr, y = observed_expected, fill = stage)) +
  geom_col(position = "dodge") +
scale_fill_manual(values = c(
    "non-Primates" = "#812e29",
    "Primates" = "#549357",
    "Simiiformes" = "#e7d968",
    "Catarrhini" = "#8b6a93",
    "Hominoidea" = "#92cbe1"
)) +
 geom_hline(aes(yintercept=1),, colour="#990000", linetype="dashed")+
  #scale_fill_viridis_d() +
  theme_bw(base_size = 14) +
  labs(x = "Chromosome", y = "Observed/expected") #+
  #theme_classic()

pdf("hs1_L1.chr.enrich.pdf",width =15 ,height =4 )
p1
dev.off()

```

