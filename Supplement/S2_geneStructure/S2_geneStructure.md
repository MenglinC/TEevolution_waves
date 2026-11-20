# S2 

## S2.1
```shell
#workplace:/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/clade_family/giggle_geneRegion/
#giggle index:/mnt/backup/Zhu_lab/xxzhang/TE_InDel/4-geneRegion/combine_index_v412/combine_chr1-X_sort_b/
giggle search -i /mnt/backup/Zhu_lab/xxzhang/TE_InDel/4-geneRegion/combine_index_v412/combine_chr1-X_sort_b/ -q "/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/clade_family/Simiiformes_L1.bed.gz" -s > Simiiformes_L1.giggle.result"

```

## S2.2
```r
#file:arrangeData_OR.R
setwd("/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/clade/giggle/")
library(vroom)
filelist <- list.files("./")
n <-length(filelist) 
test<- vroom(filelist[1],id="./",col_names= T,delim="\t")
test<-test[test$fishers_two_tail<0.05,c(2,5)]

for (i in 2:n)
{
  txt_data<-vroom(filelist[i],id="./",col_names= T,show_col_types=F,delim = "\t")
  txt_data<-txt_data[txt_data$fishers_two_tail<0.05,c(2,5)]
  test<- merge(test,txt_data,by="#file",all.x=T,all.y=T)
}
library(stringr)
superfamily<-as.data.frame(str_split(filelist,"\\.",n=3,simplify=TRUE))$V1
name<-c("TF",as.vector(superfamily))
colnames(test)<-name
write.csv(test,"../family_geneRegion_all_OR.csv",row.names = F)

```

## S2.3
```r
#file:geneRegionplot.R
setwd("F://06-16")
data<-read.csv("family_geneRegion_all_OR.csv")
##
#此处在画图的时候遇到了一些问题，non-primates和primates在画图的时候存在离群点，导致在画图的时候，离群点特别大，而其它比离群点小的数反而看不太出来了
#请教了DS，解决方法是做一个变换（原理我现在还不太清楚）
#del就是在删除离群点之后得到的结果（但是在思考这样做的合理性）
##

clade<-c("non.Primates","Primates","Simiiformes","Catarrhini","Hominoidea")
#seldata<-data[,-grep("hAT",colnames(data))]
seldata<-data
if(TRUE){
i=3
seldata2<-seldata[,grep(clade[i],colnames(seldata))]
seldata2<-asinh(seldata2 - 1) #关键在这一行，做了变换
}else{
  #Primates是特例，因为会同时带出来non-Primates的选项，应该会有解决办法，但是我懒得想了
seldata2<-seldata[,grep("^Primates",colnames(seldata))]
seldata2<-asinh(seldata2 - 1)
}

seldata2<-as.data.frame(lapply(seldata2,as.numeric))
rownames(seldata2)<-data$TF
#install.packages("pheatmap")
#seq1<-seq(min(seldata2,na.rm = T),asinh(0),by=0.1)
#n<-length(seq1)
#seq2<-seq(asinh(0)+0.1,max(seldata2,na.rm = T),by=0.1)
bk <- c(seq(min(seldata2,na.rm = T),asinh(0),length.out =10),seq(asinh(0)+0.1,max(seldata2,na.rm = T),length.out =10))
library(pheatmap)
p<-pheatmap(seldata2,cluster_rows = F,cluster_cols = F,
         color = c(colorRampPalette(colors = c("white","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         #legend_breaks=seq(-8,8,2),
         breaks=bk)

pdf(paste("pheatmap_sinh",clade[i],"pdf",sep = "."),width =6 ,height = 5)
#pdf(paste("pheatmap_sinh","Primates","pdf",sep = "."),width =4 ,height = 4)
print(p)
dev.off()


##############################################################################################
##重新画图##
##############################################################################################
#短数据转换为长数据
setwd("F://06-16")
data<-read.csv("family_geneRegion_all_OR.csv")
data[sapply(data, is.numeric)] <- lapply(data[sapply(data, is.numeric)], function(x) {
  x[x < 1] <- NA
  x
})

library(tidyr)
n<-dim(data)[2]
asinhdata<-asinh(data[,2:n] - 1)
asinhdata$region<-data$TF
data_long <- gather(asinhdata,key="clade",value="OR",-`region`)
data_long <- separate(data_long,clade,into = c("clade", "family"), sep = "_")
write.csv(data_long,"data_long.csv",row.names = F)
##############################################################################################
##开始画图
##############################################################################################
# 加载必要的库
library(ggplot2)
library(dplyr)
library(tidyr)

# 读取数据
datatmp <- read.csv("data_long.csv", stringsAsFactors = FALSE)

data_long<-datatmp[datatmp$family%in%c("L1","Alu","SVA","ERV1","ERVK"),]

family_order <-c("L1","ERV1","Alu","ERVK","SVA")

# 定义family的顺序
# family_order <- c(
#   "5S.Deu.L2", "CR1", "Dong.R4", "Gypsy", "I-Jockey", "Kolobok", "L1-Tx1", "L2",
#   "Long.Terminal.Repeat.Element", "MIR", "PIF.Harbinger", "RTE.BovB", "RTE.X",
#   "TcMar", "TcMar.Pogo", "TcMar.Tc1", "TcMar.Tc2", "Transposase", "hAT", "hAT.Ac",
#   "hAT.Blackjack", "hAT.Tag1", "hAT.Tip100", "tRNA", "tRNA.Deu", "tRNA.RTE",
#   "hAT.Charlie", "ERVL", "ERVL.MaLR", "Merlin", "TcMar.Tigger", "L1", "TcMar.Mariner",
#   "MULE.MuDR", "ERV1", "PiggyBac", "Alu", "ERVK", "SVA"
# )

# 定义region的顺序
region_order <- c(
  "intergenic", "lncRNA.downstream", "lncRNA.exon", "lncRNA.intron", 
  "lncRNA.promoter", "miRNA.2k", "protein_coding.3_UTR", "protein_coding.5_UTR",
  "protein_coding.downstream", "protein_coding.exon", "protein_coding.intron",
  "protein_coding.promoter"
)

# 定义clade的顺序
clade_order <- c("non.Primates", "Primates", "Simiiformes", "Catarrhini", "Hominoidea")

# 数据预处理：
# 1. 筛选在family_order中存在的family
# 2. 创建region_clade组合变量
# 3. 设置因子水平
plot_data <- data_long %>%
  filter(family %in% family_order) %>%
  mutate(
    family = factor(family, levels = family_order),  # rev使排序从上到下与列表一致
    region = factor(region, levels = rev(region_order)),
    clade = factor(clade, levels = rev(clade_order)),
    region_clade = paste(clade,region,sep = ":"),
    # 创建region_clade的排序变量（先按region顺序，再按clade顺序）
    region_clade = factor(
      region_clade,
      levels = apply(
        expand.grid(clade_order,region_order), 
        1, 
        paste, collapse = ":"
      )
    )
  ) %>%
  filter(!is.na(OR))  # 移除NA值
# common_levels <- c(
#   "5S.Deu.L2", "CR1", "Dong.R4", "Gypsy", "I.Jockey", "Kolobok", "L1.Tx1", "L2",
#   "Long.Terminal.Repeat.Element", "MIR", "PIF.Harbinger", "RTE.BovB", "RTE.X",
#   "TcMar", "TcMar.Pogo", "TcMar.Tc1", "TcMar.Tc2", "Transposase", "hAT", "hAT.Ac",
#   "hAT.Blackjack", "hAT.Tag1", "hAT.Tip100", "tRNA", "tRNA.Deu", "tRNA.RTE",
#   "hAT.Charlie", "ERVL", "ERVL.MaLR", "Merlin", "TcMar.Tigger", "L1", "TcMar.Mariner",
#   "MULE.MuDR", "ERV1", "PiggyBac", "Alu", "ERVK", "SVA"
# )

common_levels <-c("L1","ERV1","Alu","ERVK","SVA")

shared_scale <- function() {
  list(
    scale_x_discrete(limits = common_levels),
    theme_minimal() # 可添加统一主题
  )
}

#plot_data2<-plot_data[plot_data$region%in%c("protein_coding.promoter","miRNA.2k","lncRNA.promoter","lncRNA.exon","intergenic"),]

plot_data2<-plot_data
p2 <- ggplot(plot_data2, aes(x = family, y = clade)) +
               facet_grid(region ~ .)+
               #geom_point(aes(size = abs(OR), color = "black", shape = 21) +  # 点的大小表示|OR|
               geom_point(aes(size = abs(OR), fill = OR),shape=21,stroke = 0.2)+  # 颜色和透明度表示OR值
               scale_size_continuous(
                 range = c(0.5, 2),  # 点的大小范围
                 name = "|OR|"
               ) +
               scale_fill_gradient2(
                 low = "blue",      # 负值用蓝色
                 mid = "white",     # 中性值用白色
                 high = "#891a1b",      # 正值用红色
                 midpoint = 0,      # 中点为0
                 name = "OR Value"  # 图例名称
               ) +
               scale_alpha_continuous(guide = "none") +  # 隐藏alpha图例
               scale_size_continuous(
                  range = c(0.5, 2)   # 最小和最大尺寸（推荐2-8之间）
                  #breaks = seq(50, 350, 100)  # 可选：设置图例断点
                ) +
               labs(
                 x = "Region / Clade",
                 y = "TE Family",
                 title = "Odds Ratio (OR) of TE Family Enrichment"
               ) +
               theme_classic(base_size = 6) +
               shared_scale()+
               theme(
                 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # x轴标签旋转90度
                 axis.text.y = element_text(size = 10),
                 #panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
                 legend.position = "right",
                 legend.box = "vertical",
                 plot.title = element_text(hjust = 0.5, face = "bold")
               )
p2
pdf(paste("combine3_260129_v4.pdf"),width =4 ,height = 6)
#pdf(paste("pheatmap_sinh","Primates","pdf",sep = "."),width =4 ,height = 4)
print(p2)
dev.off()

```

