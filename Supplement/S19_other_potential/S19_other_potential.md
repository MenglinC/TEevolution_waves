# S19

## S19.1 H3K4me1

```shell
computeMatrix reference-point  --missingDataAsZero -R "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/macs/idr_FG/FG_IDR_L1_ATAC.bed" --sortUsingSamples 1  --skipZeros -S "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/macs/tf_2rep/WT_FG_SRR7525102_SRR7525103/bigwig/TF_rep1_FE.bw"  "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/macs/tf_2rep/SRR7525160_sorted.bam_SRR7525108_SRR7525109/bigwig/TF_rep2_FE.bw"   -b 2000 -a 2000  -o ./H3K4me1.FG_IDR_L1_ATAC.matrix.mat.gz
plotProfile -m "/mnt/backup02/xxzhang/ESC_dataset/mix-seq_FOXA2/results/deepTools/H3K4me1.FG_IDR_L1_ATAC.matrix.mat.gz"   --plotFileFormat   pdf  --plotWidth 5 --plotHeight 5  -out  H3K4me1.FG_IDR_L1_ATAC.v2.pdf

```

## S19.2 RNA-seq

## S19.3 nearby genes

