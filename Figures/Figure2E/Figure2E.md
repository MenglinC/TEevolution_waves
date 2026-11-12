# E

## E.1:prepare the .ini files
> 工作路径：F:\tmp\09-12-MBE示例文章

```
[x-axis]
#optional
#fontsize=20
# default is bottom meaning below the axis line
# where=top

[spacer]
# height of space in cm (optional)
height = 0.5

[genecode gene]
# height of space in cm (optional)
file=/home/xxzhang/Resource/Reference/GTF/Gencode_v42/gencode.v42.annotation.protein_coding.gtf
height = 0.5
title =
file_type = gtf
 
[spacer]
height = 0.5

[SRR2056014_sorted]
file=/mnt/backup02/xxzhang/ESC_dataset/mix-seq_Het/results/bowtie2/SRR3826166_sorted.bw

[spacer]
height = 0.5
[SRR2056014_sorted]
file=/mnt/backup02/xxzhang/ESC_dataset/mix-seq_Het/results/bowtie2/SRR3826163_sorted.bw

[spacer]
height = 0.5
[SRR2056014_sorted]
file=/mnt/backup02/xxzhang/ESC_dataset/mix-seq_Het/results/bowtie2/H3K9me3_FE.bw

[spacer]
height = 0.5
[SRR2056014_sorted]
file=/mnt/backup02/xxzhang/ESC_dataset/mix-seq_Het/results/bowtie2/TRIM28_FE.bw

[spacer]
height = 0.5
[vhighlight_track]
file = /home/xxzhang/workplace/software/pyGenomeTracks/TE/BMPR1B/Hominoidea.sorted.bed
type = vhighlight
title = Highlight
color = yellow
alpha = 0.3
zorder = -10  
line_style = solid

[Hominoidea_TE]
file=/home/xxzhang/workplace/software/pyGenomeTracks/TE/BMPR1B/Hominoidea.sorted.bed
title = Hominoidea_TE
color = #92cbe1
height = 0.5
labels = off
file_type = bed
zorder = 10

[spacer]
height = 0.5
[Hominoidea_TE]
file=/home/xxzhang/workplace/software/pyGenomeTracks/TE/NRXN3/highlight.bed
title = het_loci
color = #92cbe1
height = 0.5
labels = off
file_type = bed
zorder = 10

[spacer]
height = 0.5

[spacer]
height = 0.5
```
## E.2 plot the track

```shell
[xxzhang@mu02 TE]$ conda activate pygenometracks
(pygenometracks) [xxzhang@mu02 TE]$ mkdir NRXN3
(pygenometracks) [xxzhang@mu02 TE]$ cd NRXN3/
(pygenometracks) [xxzhang@mu02 NRXN3]$ pyGenomeTracks --tracks NRXN3.ini --region chr14:78170373-79868291  -o NRXN3.v1.pdf
```

