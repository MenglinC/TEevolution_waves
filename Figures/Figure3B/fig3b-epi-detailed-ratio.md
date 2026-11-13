```python
import csv
from collections import defaultdict, Counter
import re

def process_file(input_file, output_csv):
    """
    处理包含##标记区间和条目的文件，统计各区间下条目的类型分布
    """
    results = []
    current_interval = None
    current_entries = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # 处理每一行
    for line in lines:
        line = line.strip()
        
        # 跳过空行
        if not line:
            continue
            
        # 如果是新的区间标记行
        if line.startswith('##'):
            # 保存前一个区间的统计结果
            if current_interval is not None and current_entries:
                stats = calculate_statistics(current_interval, current_entries)
                results.append(stats)
            
            # 开始新的区间
            current_interval = line[2:].strip()  # 去掉##
            current_entries = []
        else:
            # 这是条目的行
            current_entries.append(line)
    
    # 处理最后一个区间
    if current_interval is not None and current_entries:
        stats = calculate_statistics(current_interval, current_entries)
        results.append(stats)
    
    # 输出汇总的CSV文件
    write_csv(results, output_csv)
    
    return results

def calculate_statistics(interval, entries):
    """
    计算单个区间的统计信息
    """
    # 解析区间信息
    interval_parts = interval.split('\t')
    if len(interval_parts) >= 4:
        interval_id = f"{interval_parts[0]}:{interval_parts[1]}-{interval_parts[2]} ({interval_parts[3]})"
    else:
        interval_id = interval
    
    # 统计条目
    total_count = len(entries)
    type_counter = Counter()
    
    for entry in entries:
        parts = entry.split('\t')
        if len(parts) >= 3:
            # 提取第6列的类型（列索引从0开始）
            entry_type = parts[3].strip()
            # 提取类型编号（如9_Het中的9）
            type_match = re.match(r'(\d+)_', entry_type)
            if type_match:
                type_num = type_match.group(1)
                type_counter[type_num] += 1
            else:
                type_counter[entry_type] += 1
    
    # 计算比例
    type_stats = {}
    for type_name, count in type_counter.items():
        percentage = (count / total_count * 100) if total_count > 0 else 0
        type_stats[type_name] = {
            'count': count,
            'percentage': round(percentage, 2)
        }
    
    return {
        'interval': interval_id,
        'total_count': total_count,
        'type_stats': type_stats
    }

def write_csv(results, output_file):
    """
    将统计结果写入CSV文件
    """
    # 收集所有出现过的类型
    all_types = set()
    for result in results:
        all_types.update(result['type_stats'].keys())
    
    # 排序类型
    all_types = sorted(all_types, key=lambda x: (int(x) if x.isdigit() else float('inf'), x))
    
    # 写入CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # 写入标题行
        headers = ['Interval', 'Total_Count']
        for type_name in all_types:
            headers.extend([f'Type_{type_name}_Count', f'Type_{type_name}_Percentage'])
        writer.writerow(headers)
        
        # 写入数据行
        for result in results:
            row = [result['interval'], result['total_count']]
            for type_name in all_types:
                if type_name in result['type_stats']:
                    stats = result['type_stats'][type_name]
                    row.extend([stats['count'], f"{stats['percentage']}%"])
                else:
                    row.extend([0, "0%"])
            writer.writerow(row)
    
    print(f"统计结果已保存到: {output_file}")
```


```python
import os
dir1="/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/HM_global_0122/"
os.chdir(dir1)
```


```python
clade_TE=["Hominoidea_L1","Hominoidea_ERV1","Hominoidea_Alu","Hominoidea_ERVK","Hominoidea_SVA",
           "non-Primates_L1","non-Primates_ERV1",
           "Primates_L1","Primates_ERV1","Primates_Alu",
           "Simiiformes_L1","Simiiformes_ERV1","Simiiformes_Alu","Simiiformes_ERVK",
           "Catarrhini_L1","Catarrhini_ERV1","Catarrhini_Alu","Catarrhini_ERVK"]
```


```python
if __name__ == "__main__":
    # 输入文件和输出文件
    for i in range(0,len(clade_TE)):
        input_file = "/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/HM_global_0122/giggle.detailed."+ clade_TE[i] +".txt"
        output_csv = clade_TE[i] +"_interval_statistics.csv"
        
        # 处理文件
        results = process_file(input_file, output_csv)
    
    # 打印摘要
    #print_summary(results)
```

    统计结果已保存到: Hominoidea_L1_interval_statistics.csv
    统计结果已保存到: Hominoidea_ERV1_interval_statistics.csv
    统计结果已保存到: Hominoidea_Alu_interval_statistics.csv
    统计结果已保存到: Hominoidea_ERVK_interval_statistics.csv
    统计结果已保存到: Hominoidea_SVA_interval_statistics.csv
    统计结果已保存到: non-Primates_L1_interval_statistics.csv
    统计结果已保存到: non-Primates_ERV1_interval_statistics.csv
    统计结果已保存到: Primates_L1_interval_statistics.csv
    统计结果已保存到: Primates_ERV1_interval_statistics.csv
    统计结果已保存到: Primates_Alu_interval_statistics.csv
    统计结果已保存到: Simiiformes_L1_interval_statistics.csv
    统计结果已保存到: Simiiformes_ERV1_interval_statistics.csv
    统计结果已保存到: Simiiformes_Alu_interval_statistics.csv
    统计结果已保存到: Simiiformes_ERVK_interval_statistics.csv
    统计结果已保存到: Catarrhini_L1_interval_statistics.csv
    统计结果已保存到: Catarrhini_ERV1_interval_statistics.csv
    统计结果已保存到: Catarrhini_Alu_interval_statistics.csv
    统计结果已保存到: Catarrhini_ERVK_interval_statistics.csv



```python
if __name__ == "__main__":
    # 输入文件和输出文件
    #for i in range(0,len(clade_TE)):
        input_file = "giggle.detailed.Hominoidea.v3.genes.centre.TE.txt"
        output_csv = "Hominoidea.v3.genes.centre.TE_interval_statistics.csv"
        
        # 处理文件
        results = process_file(input_file, output_csv)
    
    # 打印摘要
    #print_summary(results)
```

    统计结果已保存到: Hominoidea.v3.genes.centre.TE_interval_statistics.csv



```python

```
