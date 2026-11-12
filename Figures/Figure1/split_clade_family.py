from collections import defaultdict

# 读取文件并分组
groups = defaultdict(list)
with open('result2_by_newclade_family.bed') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 8: 
            continue  # 跳过无效行
        key = f"{parts[6]}_{parts[7]}"  # 组合第7列和第8列
        groups[key].append(line)

# 写入分组文件
for key, lines in groups.items():
    with open(f"./clade_family/{key}.bed", 'w') as f_out:
        f_out.writelines(lines)
        