# A.
## A.1

## A.2
```r

```R
data<-read.table("/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/HM_global_0122/all.giggle.result",sep = "\t",fill=TRUE,na.strings = "",quote="",header=T)

library(tidyr)
library(dplyr)
data1<- data %>%
  separate(group, 
           into = c("clade", "TE"), 
           sep = "_(?=[^_]*$)",  # 正则表达式：匹配最后一个下划线
           remove = FALSE,      # 保留原列
           fill = "right")      # 如果拆分不出两部分，右边填充NA
data2<- data1 %>%
  separate(file, 
           into = c("tissue", "state"), 
           sep = "_(?=[^_]*$)",  # 正则表达式：匹配最后一个下划线
           remove = FALSE,      # 保留原列
           fill = "right")      # 如果拆分不出两部分，右边填充NA
df_summary <- data2 %>%
  group_by(state, group) %>%
  summarise(
    mean_num = mean(overlaps, na.rm = TRUE), #修改此處
    .groups = "drop"
  )
df_summary1<- df_summary %>%
  separate(group, 
           into = c("clade", "TE"), 
           sep = "_(?=[^_]*$)",  # 正则表达式：匹配最后一个下划线
           remove = FALSE,      # 保留原列
           fill = "right")      # 如果拆分不出两部分，右边填充NA
df_summary2<-df_summary1[df_summary1$TE%in%c("Alu","ERVK","ERV1","SVA","L1"),]
library(dplyr)
df_summary2 <- df_summary2 %>%
  group_by(group) %>%
  mutate(prop_state = mean_num / sum(mean_num, na.rm = TRUE)) %>%
  ungroup()
df_summary2$TE<-factor(df_summary2$TE,levels=c("L1","ERV1","Alu","ERVK","SVA"))
df_summary2$clade<-factor(df_summary2$clade,levels=c("non-Primates","Primates","Simiiformes","Catarrhini","Hominoidea"))
#df_summary2$state<-factor(df_summary2$state,levels=c("TssA","TssAFlnk","TxFlnk","Tx","TxWk",
 #                                                   "EnhG","Enh","ZNF/Rpts","Het","TssBiv",
  #                                                  "BivFlnk","EnhBiv","ReprPC","ReprPCWk","Quies"))
df_summary2$state<-factor(df_summary2$state,levels=c("Active.TSS","Flanking.Active.TSS","Transcr.at.gene.5.and.3","Strong.transcription","Weak.transcription",
                                                    "Genic.enhancers","Enhancers","ZNF.genes.and.repeats","Heterochromatin","Bivalent.Poised.TSS",
                                                    "Flanking.Bivalent.TSS.Enh","Bivalent.Enhancer","Repressed.PolyComb","Weak.Repressed.PolyComb","Quiescent.Low"))
library(ggplot2)
options(repr.plot.width =12 , repr.plot.height = 4)
p<-ggplot(df_summary2, aes(fill=state, y=mean_num, x=TE)) + 
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(values = c("Active.TSS" = "#d72929", 
                                 "Flanking.Active.TSS" = "#a4e298", 
                                 "Transcr.at.gene.5.and.3" = "#f0e442",
                                "Strong.transcription" = "#2079b4", 
                                 "Weak.transcription" = "#aec6e8", 
                                 "Genic.enhancers" = "#afd063",
                                "Enhancers" = "#fcb23d", 
                                 "ZNF.genes.and.repeats" = "#626262", 
                                 "Heterochromatin" = "#747474",
                                "Bivalent.Poised.TSS" = "#858585", 
                                 "Flanking.Bivalent.TSS.Enh" = "#979797", 
                                 "Bivalent.Enhancer" = "#a8a8a8",
                                "Repressed.PolyComb" = "#b9b9b9", 
                                 "Weak.Repressed.PolyComb" = "#cbcbcb", 
                                 "Quiescent.Low" = "#dcdcdc"))+
    labs(title = "",x = "Family",y = "Percentage") +
    theme_classic()+
    facet_grid(. ~ clade)+
    theme(legend.position="bottom")

pdf("family_clade_percentage_mean.pdf",width =10 ,height =4 )
p
dev.off()

```
