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
