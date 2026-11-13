# F

#workplace:/mnt/backup02/xxzhang/ESC_dataset/ChIP-seq_ZNF281/ENCODE/
#file:ENCFF855CVL.bigWig

```shell
computeMatrix scale-regions -S ENCFF855CVL.bigWig -R "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1HS_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA2_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA3_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA4_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA5_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA6_6kb.bed" "/home/xxzhang/workplace/project/CRISPRa/basic/bed/6kb_fasta/L1PA7_6kb.bed"  -b  2000  -a  2000  --regionBodyLength  6000  --binSize 50  --averageTypeBins  mean  -p  8  -o  ZNF281.matrix.gz
plotProfile  -m  ZNF281.matrix.gz   -out ZNF281.profile.pdf –regionsLabel  "L1HS" "L1PA2" "L1PA3" "L1PA4" "L1PA5" "L1PA6" "L1PA7"  --samplesLabel  "Fold enrichment"  --plotType se --averageType mean –startLabel  "5-prime end" --endLabel  "3-prime end"  --yAxisLabel "Fold enrichment over control"

```
