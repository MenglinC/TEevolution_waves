import argparse
from collections import defaultdict
from typing import List, Tuple, Dict

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='查找TE中间点是否在基因区域内')
    parser.add_argument("-a",
                       help="基因的BED文件路径",
                       dest="file1",
                       type=str,
                       default="proteinGene.bed")
    parser.add_argument("-b",
                       help="TE的BED文件路径",
                       dest="file2",
                       type=str,
                       default="Hominoidea.bed")
    return parser.parse_args()

def load_genes(filename: str) -> Dict[str, List[Tuple[int, int, str, str]]]:
    """
    加载基因数据并建立索引
    返回字典: {chromosome: [(start, end, name, strand), ...]}
    """
    gene_dict = defaultdict(list)
    
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():  # 跳过空行
                fields = line.strip().split()
                if len(fields) >=3:  # 确保有足够的列
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    name = fields[3]
                    strand = fields[4]
                    gene_dict[chrom].append((start, end, name, strand))
    
    # 按起始位置排序，便于后续二分查找
    for chrom in gene_dict:
        gene_dict[chrom].sort(key=lambda x: x[0])
    
    return gene_dict

def find_overlapping_genes(gene_list: List[Tuple[int, int, str, str]], 
                        midpoint: int) -> List[Tuple[int, int, str, str]]:
    """
    在有序的基因列表中查找包含指定点的基因
    使用线性搜索，因为基因区域可能重叠
    """
    overlapping = []
    for gene in gene_list:
        start_gene, end_gene, name, strand = gene
        # 如果当前TE的起始位置已经超过中点，后面的TE也都会超过（因为已排序）
        if start_gene > midpoint:
            break
        if start_gene <= midpoint < end_gene:
            overlapping.append(gene)
    return overlapping

def process_tes(te_file: str, gene_dict: Dict[str, List[Tuple[int, int, str, str]]]):
    """处理基因文件，查找与TE重叠的基因"""
    with open(te_file, 'r') as f:
        for line in f:
            if not line.strip():  # 跳过空行
                continue
                
            fields = line.strip().split()
            if len(fields) < 3:  # 确保有足够的列
                continue
                
            chrom_te = fields[0]
            start_te = int(fields[1])
            end_te = int(fields[2])
            sclass = fields[3]
            subfamily = fields[7]
            family = fields[6]
            # 计算基因中间点（整数运算）
            midpoint = start_te + (end_te - start_te) // 2
            
            # 只在对应染色体的TE中查找
            if chrom_te in gene_dict:
                overlapping_genes = find_overlapping_genes(gene_dict[chrom_te], midpoint)
                
                for start_gene, end_gene,name,strand in overlapping_genes:
                    # 输出结果，保持与原格式一致
                    print(f"{chrom_te}\t{start_gene}\t{end_gene}\t{name}\t{strand}\t|\t"
                          f"{chrom_te}\t{start_te}\t{end_te}\t"
                          f"{sclass}\t{subfamily}\t{family}")

def main():
    """主函数"""
    args = parse_args()
    
    # 预先加载并索引TE数据
    print(f"正在加载基因数据: {args.file1}", file=sys.stderr)
    gene_dict = load_genes(args.file1)
    
    # 处理基因数据
    print(f"正在处理TE数据: {args.file2}", file=sys.stderr)
    process_tes(args.file2, gene_dict)

if __name__ == "__main__":
    import sys
    main()
