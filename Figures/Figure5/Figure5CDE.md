# C

```shell
computeMatrix reference-point  -p 32 --missingDataAsZero -R "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1HS_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA2_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA3_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA4_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA5_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA6_6kb.bed"  --sortUsingSamples 1 --skipZeros -S "/mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_L1_TF/ENCFF008ABX.bigWig" "/mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_L1_TF/ENCFF348LLH.bigWig" "/mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_L1_TF/ENCFF737TUE.bigWig"   -b 2000 -a 2000  -o ./FOXA.6kb.L1-6.matrix.mat.gz
plotProfile -m FOXA.6kb.L1-6.matrix.mat.gz  --plotFileFormat  pdf --plotWidth 5 --plotHeight 5  -out  FOXA.6kb.L1-6.v2.pdf

```

# D

```shell
computeMatrix reference-point  -p 32 --missingDataAsZero -R "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1HS_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA2_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA3_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA4_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA5_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA6_6kb.bed"  --sortUsingSamples 1 --skipZeros -S "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/deepTools/ATAC-seq_BW/SRR7126315_sorted.bw"  "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/deepTools/ATAC-seq_BW/SRR7126316_sorted.bw"   -b 2000 -a 2000  -o ./FOXA.ATAC-seq.6kb.L1-6.matrix.mat.gz
plotProfile -m "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/deepTools/FOXA.ATAC-seq.6kb.L1-6.matrix.mat.gz"  --plotFileFormat  pdf --plotWidth 5 --plotHeight 5  -out  FOXA2.ATAC-seq.1516.6kb.L1-6.v2.pdf

```

# E
## Prepare input

```shell
#workplace:/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/macs/idr_FG/
bedtools intersect -a L1_all.bed  -b FG_IDR_0.05.bed  -wa -wb >FG_IDR_L1.bed
bedtools intersect -a FG_IDR_L1.bed  -b FG_ATAC_KO_down.sorted.bed  -wa -wb >FG_IDR_L1_ATAC.bed
bedtools intersect -a FG_IDR_L1_ATAC.bed  -b L1_1-3.bed  -wa -wb |awk -v FS="\t" -v OFS="\t" '{print $18,$19,$20,$21,$22,$23}' >FG_IDR_L1_ATAC.test.bed
```

## Plot
```shell
computeMatrix reference-point  --missingDataAsZero -R "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/macs/idr_FG/FG_IDR_L1_ATAC.bed" --sortUsingSamples 1  --skipZeros -S "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/macs/tf_2rep/WT_FG_SRR7525114_SRR7525115/bigwig/TF_rep1_FE.bw" "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/macs/tf_2rep/SRR7525160_sorted.bam_SRR7525120_SRR7525121/bigwig/TF_rep2_FE.bw"   -b 2000 -a 2000  -o ./H3K27ac.FG_IDR_L1_ATAC.matrix.mat.gz
plotProfile -m "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/deepTools/H3K27ac.FG_IDR_L1_ATAC.matrix.mat.gz"   --plotFileFormat   pdf  --plotWidth 5 --plotHeight 5  -out  H3K27ac.FG_IDR_L1_ATAC.v2.pdf
```
