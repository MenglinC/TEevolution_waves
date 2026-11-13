#!/usr/bin/env bash
set -euo pipefail

input=${1:-Hominoidea.v3.genes.txt}
output=${2:-Hominoidea.v3.genes.summary.tsv}

awk -F'\t' 'BEGIN {
    OFS="\t"
    print "gene_name", "gene_length", "genomic_location", "row_count", \
          "col10_types", "col10_type_counts", "col11_types", "col11_type_counts"
}
{
    # 兼容 Windows CRLF 行尾
    sub(/\r$/, "", $NF)

    gene=$4
    len=$3-$2
    loc=$1 ":" $2 "-" $3
    type10=$10
    type11=$11

    # 记录基因首次出现顺序
    if (!(gene in gene_seen)) {
        gene_seen[gene]=1
        ngene++
        gene_order[ngene]=gene
    }

    row_count[gene]++

    # 同一基因名可能对应多个坐标或长度，均去重后保留
    k=gene SUBSEP len
    if (!(k in len_seen)) {
        len_seen[k]=1
        gene_length[gene]=(gene_length[gene] == "" ? len : gene_length[gene] "," len)
    }

    k=gene SUBSEP loc
    if (!(k in loc_seen)) {
        loc_seen[k]=1
        genomic_location[gene]=(genomic_location[gene] == "" ? loc : genomic_location[gene] ";" loc)
    }

    # 第10列：记录类型顺序及每种类型出现次数
    k=gene SUBSEP type10
    if (!(k in type10_seen)) {
        type10_seen[k]=1
        type10_n[gene]++
        idx=type10_n[gene]
        type10_value[gene SUBSEP idx]=type10
    }
    type10_count[k]++

    # 第11列：记录类型顺序及每种类型出现次数
    k=gene SUBSEP type11
    if (!(k in type11_seen)) {
        type11_seen[k]=1
        type11_n[gene]++
        idx=type11_n[gene]
        type11_value[gene SUBSEP idx]=type11
    }
    type11_count[k]++
}
END {
    for (i=1; i<=ngene; i++) {
        gene=gene_order[i]

        types10=""
        counts10=""
        for (j=1; j<=type10_n[gene]; j++) {
            value=type10_value[gene SUBSEP j]
            types10=(types10 == "" ? value : types10 "," value)
            counts10=(counts10 == "" ? value "(" type10_count[gene SUBSEP value] ")" : \
                                      counts10 "," value "(" type10_count[gene SUBSEP value] ")")
        }

        types11=""
        counts11=""
        for (j=1; j<=type11_n[gene]; j++) {
            value=type11_value[gene SUBSEP j]
            types11=(types11 == "" ? value : types11 "," value)
            counts11=(counts11 == "" ? value "(" type11_count[gene SUBSEP value] ")" : \
                                      counts11 "," value "(" type11_count[gene SUBSEP value] ")")
        }

        print gene, gene_length[gene], genomic_location[gene], row_count[gene], \
              types10, counts10, types11, counts11
    }
}' "$input" > "$output"

echo "Output written to: $output" >&2
