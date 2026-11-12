#!/bin/bash

# 检查参数数量
if [ $# -ne 2 ]; then
    echo "Usage: $0 <newclade.annot.csv> <input.gtf>"
    exit 1
fi

csv_file="$1"
gtf_file="$2"

# 步骤1: 从CSV文件中提取subfamily, family, newclade映射
awk -F, '
NR > 1 {
    # 移除字段两侧的双引号
    gsub(/^"|"$/, "", $1);
    gsub(/^"|"$/, "", $4);
    gsub(/^"|"$/, "", $6);
    print $1 "\t" $4 "\t" $6
}' "$csv_file" > mapping.tmp

# 步骤2: 处理GTF文件并生成中间BED文件
awk -F'\t' '
BEGIN {
    OFS = "\t";
    # 加载映射关系
    while (getline < "mapping.tmp") {
        subfamily = $1;
        family = $2;
        newclade = $3;
        map[subfamily] = 1;
        family_map[subfamily] = family;
        newclade_map[subfamily] = newclade;
    }
    close("mapping.tmp");
}
{
    # 跳过注释行
    if ($0 ~ /^#/) next;

    # 在属性字段中查找gene_id
    if ($9 ~ /gene_id *"[^"]+";/) {
        # 提取gene_id值
        split($9, attr, ";");
        for (i in attr) {
            if (attr[i] ~ /gene_id *"/) {
                split(attr[i], parts, "\"");
                gene_id = parts[2];
                if (gene_id in map) {
                    # 转换坐标: GTF(1-based) -> BED(0-based)
                    chrom = $1;
                    start = $4 - 1;
                    end = $5;
                    strand = $7;
                    
                    # 获取映射信息
                    family = family_map[gene_id];
                    newclade = newclade_map[gene_id];
                    
                    # 构建BED行
                    bed_line = chrom OFS start OFS end OFS gene_id OFS "0" OFS strand;
                    
                    # 输出到中间文件
                    print bed_line, newclade > "result1.tmp";
                    print bed_line, newclade, family > "result2.tmp";
                }
            }
        }
    }
}' "$gtf_file"

# 步骤3: 排序并生成最终结果
sort -k7,7 "result1.tmp" -o "result1_by_newclade.bed"
sort -k7,7 -k8,8 "result2.tmp" -o "result2_by_newclade_family.bed"

# 清理临时文件
rm -f mapping.tmp result1.tmp result2.tmp

echo "Processing completed. Results saved to:"
echo "1. result1_by_newclade.bed (grouped by newclade)"
echo "2. result2_by_newclade_family.bed (grouped by newclade and family)"
