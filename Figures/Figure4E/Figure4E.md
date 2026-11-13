# E

## E.1.0 download data

```shell
prefetch SRR3112448/SRR3112447

```

## E.1.1

#workplace:/mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_zfp281/
```shell

#genome=/home/xxzhang/Resource/Reference/Genome/homo_sapiens/hg38/bowtie2/hg38_dm6 #bowtie2的索引目录
genome=/mnt/backup02/xxzhang/Reference/Mus_musculus/bowtie2/GRCm39
ID=/mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_zfp281/meta/SRR_Acc_List.txt #样本的ID信息
fastqdir=/mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_zfp281/raw_data/fastq/
bowdir=/mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_zfp281/results/bowtie2/ #bowtie2比对结果文件的位置
outdir=/mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_zfp281/results/deepTools/

#Run bowtie2，进行比对，single
cat ${ID} | parallel \
        bowtie2 -p 6 -q --local -x ${genome} -U ${fastqdir}/{}.fastq.gz -S ${bowdir}/{}_unsorted.sam
#paired
#cat ${ID} | parallel \
#        bowtie2 -p 6 -q --local -x ${genome} -1 ${fastqdir}/{}_1.fastq.gz -2 ${fastqdir}/{}_2.fastq.gz -S ${bowdir}/{}_unsorted.sam

# Create BAM from SAM，进行格式转换
cat ${ID} | parallel \
        samtools view -h -S -b -@ 6 -o ${bowdir}/{}_unsorted.bam ${bowdir}/{}_unsorted.sam

# Sort BAM file by genomic coordinates，将bam文件排序
cat ${ID} | parallel \
        sambamba sort -t 6 -o ${bowdir}/{}_sorted.bam ${bowdir}/{}_unsorted.bam

# Filter out duplicates，过滤多重比对序列，只保留唯一比对序列
#cat ${ID} | parallel \
#        sambamba view -h -t 6 -f bam -F "[XS] == null and not unmapped " ${bowdir}/{}_sorted.bam > ${bowdir}/{}_aln.bam

# Create indices for all the bam files for visualization and QC，建立bam文件的索引
#cat ${ID} | parallel \
#        samtools index ${bowdir}/{}_aln.bam

```



## E.1.2

```shell

conda activate macs2_python2.7

macs2 callpeak \
-t SRR3112448.sorted.bam \
-c SRR3112447.sorted.bam \ #input
-g hs \
-n ZNF281 \
-B \
--SPMR

macs2 bdgcmp \
-t ZNF281_treat_pileup.bdg \
-c ZNF281_control_lambda.bdg \
-m FE \
-o ZNF281_FE.bdg

LC_COLLATE=C sort -k1,1 -k2,2n \
ZNF281_FE.bdg \
-o ZNF281_FE.sorted.bdg

在转化为bw文件：
bedGraphToBigWig \
ZNF281_FE.sorted.bdg \
hg38.chrom.sizes \
ZNF281_FE.bw
```

## E.1.3

```shell
parallel computeMatrix reference-point --regionBodyLength 1000 -p 32 --startLabel "5'" --endLabel "3'" --missingDataAsZero -R "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdTf_I.bed"  "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdTf_II.bed"  "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdTf_III.bed" "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdA_I.bed"  "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdA_II.bed"  "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdA_III.bed"  "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdA_IV.bed"  "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdA_V.bed"  "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdA_VI.bed"  "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdA_VII.bed"   "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdGf_I.bed"  "/mnt/backup02/xxzhang/Reference/Mus_musculus/TEtranscript/L1/bed/L1MdGf_II.bed" --sortUsingSamples 1 --binSize 100 --skipZeros -S ZNF281_FE.bw -b 100 -a 100 --samplesLabel ZNF281 -o ZNF281.L1.matrix.mat.gz
parallel plotHeatmap -m ZNF281.L1.matrix.mat.gz --colorMap summer --regionsLabel "L1MdTf_I" "L1MdTf_II" "L1MdTf_III" "L1MdA_I" "L1MdA_II" "L1MdA_III" "L1MdA_IV" "L1MdA_V" "L1MdA_VI" "L1MdA_VII" "L1MdGf_I" "L1MdGf_II" -o ZNF281.L1.pdf
```

