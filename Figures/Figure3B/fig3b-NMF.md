```R
setwd("/home/xxzhang/Resource/Reference/GTF/TE_gtf/hg38_rmsk_TE_20240807/HM_global_0122/")
```


```R
#data<-read.csv("Hominoidea_ERV1_interval_statistics.csv")
```


```R
#head(data)
```


```R
# 加载必要的包（如果没有安装，请先运行 install.packages("dplyr")）
library(dplyr)
# 获取所有匹配的文件路径
file_list <- list.files(pattern = "*_interval_statistics.csv$", full.names = TRUE)
#file_list <- list.files(pattern = "Hominoidea_ERV1_interval_statistics.csv", full.names = TRUE)
# 读取所有文件到一个列表
data_list <- lapply(file_list, read.csv, stringsAsFactors = FALSE)
# 合并所有数据框
combined_data <- bind_rows(data_list)
# 查看结果
head(combined_data)
```


<table class="dataframe">
<caption>A data.frame: 6 × 32</caption>
<thead>
	<tr><th></th><th scope=col>Interval</th><th scope=col>Total_Count</th><th scope=col>Type_1_Count</th><th scope=col>Type_1_Percentage</th><th scope=col>Type_2_Count</th><th scope=col>Type_2_Percentage</th><th scope=col>Type_3_Count</th><th scope=col>Type_3_Percentage</th><th scope=col>Type_4_Count</th><th scope=col>Type_4_Percentage</th><th scope=col>⋯</th><th scope=col>Type_11_Count</th><th scope=col>Type_11_Percentage</th><th scope=col>Type_12_Count</th><th scope=col>Type_12_Percentage</th><th scope=col>Type_13_Count</th><th scope=col>Type_13_Percentage</th><th scope=col>Type_14_Count</th><th scope=col>Type_14_Percentage</th><th scope=col>Type_15_Count</th><th scope=col>Type_15_Percentage</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr10:10160728-10162335 (LTR12C)  </td><td>159</td><td>0</td><td>0%   </td><td>0</td><td>0%   </td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>⋯</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td> 0</td><td>0%   </td><td> 1</td><td>0.63% </td><td> 98</td><td>61.64%</td></tr>
	<tr><th scope=row>2</th><td>chr10:102199775-102201425 (LTR12C)</td><td>160</td><td>0</td><td>0%   </td><td>0</td><td>0%   </td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>⋯</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>10</td><td>6.25%</td><td>39</td><td>24.38%</td><td> 82</td><td>51.25%</td></tr>
	<tr><th scope=row>3</th><td>chr10:103205440-103205872 (LTR7Y) </td><td>150</td><td>2</td><td>1.33%</td><td>8</td><td>5.33%</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>⋯</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td> 0</td><td>0%   </td><td>11</td><td>7.33% </td><td> 99</td><td>66.0% </td></tr>
	<tr><th scope=row>4</th><td>chr10:103213230-103213657 (LTR7Y) </td><td>131</td><td>0</td><td>0%   </td><td>0</td><td>0%   </td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>⋯</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td> 0</td><td>0%   </td><td> 2</td><td>1.53% </td><td>113</td><td>86.26%</td></tr>
	<tr><th scope=row>5</th><td>chr10:103521833-103522264 (LTR7Y) </td><td>131</td><td>0</td><td>0%   </td><td>1</td><td>0.76%</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>⋯</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td> 3</td><td>2.29%</td><td>24</td><td>18.32%</td><td> 96</td><td>73.28%</td></tr>
	<tr><th scope=row>6</th><td>chr10:104244986-104246358 (LTR12C)</td><td>156</td><td>0</td><td>0%   </td><td>0</td><td>0%   </td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>⋯</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td> 0</td><td>0%   </td><td>46</td><td>29.49%</td><td> 75</td><td>48.08%</td></tr>
</tbody>
</table>




```R
library(dplyr)

# 获取所有匹配的文件路径
file_list <- list.files(pattern = "*_interval_statistics.csv$", full.names = TRUE)

# 读取所有文件，并为每个文件添加一个包含文件名关键部分的列
data_list <- lapply(file_list, function(file) {
  # 读取CSV文件
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # 提取想要的字符串：去掉路径和末尾的 "_interval_statistics.csv"
  # 例如 "path/Hominoidea_ERV1_interval_statistics.csv" -> "Hominoidea_ERV1"
  name_part <- sub("_interval_statistics\\.csv$", "", basename(file))
  
  # 将提取的字符串作为新列添加到数据框
  df$group <- name_part   # 你可以将列名改为更合适的名称，如 "dataset" 或 "source"
  
  df
})

# 合并所有数据框
combined_data <- bind_rows(data_list)

# 查看结果
head(combined_data)
```


<table class="dataframe">
<caption>A data.frame: 6 × 33</caption>
<thead>
	<tr><th></th><th scope=col>Interval</th><th scope=col>Total_Count</th><th scope=col>Type_1_Count</th><th scope=col>Type_1_Percentage</th><th scope=col>Type_2_Count</th><th scope=col>Type_2_Percentage</th><th scope=col>Type_3_Count</th><th scope=col>Type_3_Percentage</th><th scope=col>Type_4_Count</th><th scope=col>Type_4_Percentage</th><th scope=col>⋯</th><th scope=col>Type_11_Percentage</th><th scope=col>Type_12_Count</th><th scope=col>Type_12_Percentage</th><th scope=col>Type_13_Count</th><th scope=col>Type_13_Percentage</th><th scope=col>Type_14_Count</th><th scope=col>Type_14_Percentage</th><th scope=col>Type_15_Count</th><th scope=col>Type_15_Percentage</th><th scope=col>group</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr10:100012056-100012304 (AluY)</td><td>130</td><td>1</td><td>0.77%</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>  0</td><td>0%    </td><td>⋯</td><td>0%</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>19</td><td>14.62%</td><td> 91</td><td>70.0% </td><td>Catarrhini_Alu</td></tr>
	<tr><th scope=row>2</th><td>chr10:100021880-100022162 (AluY)</td><td>127</td><td>0</td><td>0%   </td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>  0</td><td>0%    </td><td>⋯</td><td>0%</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>18</td><td>14.17%</td><td> 86</td><td>67.72%</td><td>Catarrhini_Alu</td></tr>
	<tr><th scope=row>3</th><td>chr10:1000235-1000536 (AluY)    </td><td>137</td><td>0</td><td>0%   </td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>114</td><td>83.21%</td><td>⋯</td><td>0%</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td> 0</td><td>0%    </td><td>  1</td><td>0.73% </td><td>Catarrhini_Alu</td></tr>
	<tr><th scope=row>4</th><td>chr10:100023937-100024247 (AluY)</td><td>127</td><td>0</td><td>0%   </td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>  0</td><td>0%    </td><td>⋯</td><td>0%</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>13</td><td>10.24%</td><td> 94</td><td>74.02%</td><td>Catarrhini_Alu</td></tr>
	<tr><th scope=row>5</th><td>chr10:100031501-100031793 (AluY)</td><td>127</td><td>0</td><td>0%   </td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>  0</td><td>0%    </td><td>⋯</td><td>0%</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td> 6</td><td>4.72% </td><td>114</td><td>89.76%</td><td>Catarrhini_Alu</td></tr>
	<tr><th scope=row>6</th><td>chr10:100063112-100063420 (AluY)</td><td>127</td><td>0</td><td>0%   </td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>  1</td><td>0.79% </td><td>⋯</td><td>0%</td><td>0</td><td>0%</td><td>0</td><td>0%</td><td>56</td><td>44.09%</td><td> 66</td><td>51.97%</td><td>Catarrhini_Alu</td></tr>
</tbody>
</table>




```R
combined_data3<-combined_data[,c(1,33)]
```


```R
library(tidyr)
combined_data3<-separate(combined_data[,c(1,33)],group, into = c("lineage","family"), sep = "[_]")
head(combined_data3)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>Interval</th><th scope=col>lineage</th><th scope=col>family</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr10:100012056-100012304 (AluY)</td><td>Catarrhini</td><td>Alu</td></tr>
	<tr><th scope=row>2</th><td>chr10:100021880-100022162 (AluY)</td><td>Catarrhini</td><td>Alu</td></tr>
	<tr><th scope=row>3</th><td>chr10:1000235-1000536 (AluY)    </td><td>Catarrhini</td><td>Alu</td></tr>
	<tr><th scope=row>4</th><td>chr10:100023937-100024247 (AluY)</td><td>Catarrhini</td><td>Alu</td></tr>
	<tr><th scope=row>5</th><td>chr10:100031501-100031793 (AluY)</td><td>Catarrhini</td><td>Alu</td></tr>
	<tr><th scope=row>6</th><td>chr10:100063112-100063420 (AluY)</td><td>Catarrhini</td><td>Alu</td></tr>
</tbody>
</table>




```R
write.csv(combined_data3, "merge_intervals_group.csv", row.names = TRUE)
```


```R
dim(combined_data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>3330</li><li>32</li></ol>




```R
combined_data2<-combined_data[,c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)]/combined_data[,2]
```


```R
rownames(combined_data2)<-combined_data[,1]
```


```R
head(combined_data2)
```


<table class="dataframe">
<caption>A data.frame: 6 × 15</caption>
<thead>
	<tr><th></th><th scope=col>Type_1_Count</th><th scope=col>Type_2_Count</th><th scope=col>Type_3_Count</th><th scope=col>Type_4_Count</th><th scope=col>Type_5_Count</th><th scope=col>Type_6_Count</th><th scope=col>Type_7_Count</th><th scope=col>Type_8_Count</th><th scope=col>Type_9_Count</th><th scope=col>Type_10_Count</th><th scope=col>Type_11_Count</th><th scope=col>Type_12_Count</th><th scope=col>Type_13_Count</th><th scope=col>Type_14_Count</th><th scope=col>Type_15_Count</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>chr10:10160728-10162335 (LTR12C)</th><td>0.00000000</td><td>0.000000000</td><td>0</td><td>0</td><td>0.006289308</td><td>0</td><td>0.00000000</td><td>0.025157233</td><td>0.345911950</td><td>0</td><td>0</td><td>0</td><td>0.00000000</td><td>0.006289308</td><td>0.6163522</td></tr>
	<tr><th scope=row>chr10:102199775-102201425 (LTR12C)</th><td>0.00000000</td><td>0.000000000</td><td>0</td><td>0</td><td>0.006250000</td><td>0</td><td>0.01875000</td><td>0.012500000</td><td>0.143750000</td><td>0</td><td>0</td><td>0</td><td>0.06250000</td><td>0.243750000</td><td>0.5125000</td></tr>
	<tr><th scope=row>chr10:103205440-103205872 (LTR7Y)</th><td>0.01333333</td><td>0.053333333</td><td>0</td><td>0</td><td>0.046666667</td><td>0</td><td>0.14666667</td><td>0.000000000</td><td>0.006666667</td><td>0</td><td>0</td><td>0</td><td>0.00000000</td><td>0.073333333</td><td>0.6600000</td></tr>
	<tr><th scope=row>chr10:103213230-103213657 (LTR7Y)</th><td>0.00000000</td><td>0.000000000</td><td>0</td><td>0</td><td>0.038167939</td><td>0</td><td>0.01526718</td><td>0.015267176</td><td>0.053435115</td><td>0</td><td>0</td><td>0</td><td>0.00000000</td><td>0.015267176</td><td>0.8625954</td></tr>
	<tr><th scope=row>chr10:103521833-103522264 (LTR7Y)</th><td>0.00000000</td><td>0.007633588</td><td>0</td><td>0</td><td>0.030534351</td><td>0</td><td>0.01526718</td><td>0.000000000</td><td>0.007633588</td><td>0</td><td>0</td><td>0</td><td>0.02290076</td><td>0.183206107</td><td>0.7328244</td></tr>
	<tr><th scope=row>chr10:104244986-104246358 (LTR12C)</th><td>0.00000000</td><td>0.000000000</td><td>0</td><td>0</td><td>0.051282051</td><td>0</td><td>0.06410256</td><td>0.006410256</td><td>0.102564103</td><td>0</td><td>0</td><td>0</td><td>0.00000000</td><td>0.294871795</td><td>0.4807692</td></tr>
</tbody>
</table>




```R
write.csv(combined_data2, "Hominoidea_ERV1_intervals_ratio.csv", row.names = TRUE)
```


```R
write.csv(combined_data2, "merged_intervals.csv", row.names = TRUE)
```


```R
dim(combined_data2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>2284835</li><li>15</li></ol>




```R
BiocManager::install("devtools")
```

    'getOption("repos")' replaces Bioconductor standard repositories, see
    'help("repositories", package = "BiocManager")' for details.
    Replacement repositories:
        CRAN: https://cran.r-project.org
    
    Bioconductor version 3.20 (BiocManager 1.30.26), R 4.4.3 (2025-02-28)
    
    Installing package(s) 'devtools'
    
    还安装依赖关系‘systemfonts’, ‘textshaping’, ‘gert’, ‘ragg’, ‘usethis’, ‘pak’, ‘pkgbuild’, ‘pkgdown’, ‘pkgload’, ‘rcmdcheck’, ‘roxygen2’, ‘testthat’
    
    
    Warning message in install.packages(...):
    “安装程序包‘systemfonts’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘gert’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘pak’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘pkgbuild’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘textshaping’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘usethis’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘pkgload’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘rcmdcheck’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘ragg’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘roxygen2’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘testthat’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘pkgdown’时退出狀態的值不是0”
    Warning message in install.packages(...):
    “安装程序包‘devtools’时退出狀態的值不是0”
    更新'.Library'里的HTML程序包列表
    
    Making 'packages.html' ...
     做完了。
    
    Old packages: 'base64enc', 'BH', 'BiocManager', 'BiocParallel', 'Biostrings',
      'blob', 'broom', 'bslib', 'car', 'carData', 'checkmate', 'cluster', 'cpp11',
      'data.table', 'DBI', 'dbplyr', 'digest', 'doBy', 'downlit', 'dplyr',
      'dtplyr', 'emmeans', 'factoextra', 'fansi', 'fitdistrplus', 'flashClust',
      'foreign', 'futile.logger', 'future', 'future.apply', 'gargle',
      'GenomeInfoDb', 'ggplot2', 'ggpubr', 'ggsci', 'globals', 'gplots', 'highr',
      'Hmisc', 'htmltools', 'httr', 'httr2', 'igraph', 'IRanges', 'irlba',
      'isoband', 'knitr', 'later', 'lattice', 'leidenbase', 'listenv', 'lme4',
      'lsa', 'lubridate', 'MatrixGenerics', 'mgcv', 'mvtnorm', 'openssl',
      'parallelly', 'plotly', 'png', 'progressr', 'promises', 'purrr', 'rappdirs',
      'rbibutils', 'Rcpp', 'RcppAnnoy', 'RcppArmadillo', 'Rdpack', 'reformulas',
      'reshape2', 'reticulate', 'ROCR', 'RSQLite', 'rstatix', 'rstudioapi', 'S7',
      'sctransform', 'Seurat', 'SeuratObject', 'shiny', 'sp', 'spam',
      'SparseArray', 'spatstat.data', 'spatstat.explore', 'spatstat.geom',
      'spatstat.random', 'spatstat.univar', 'spatstat.utils', 'stringr',
      'survival', 'tibble', 'tidyr', 'tidytree', 'timechange', 'tinytex', 'uuid',
      'vctrs', 'vegan', 'viridisLite', 'vroom', 'WGCNA', 'xfun', 'XML', 'xtable',
      'yaml', 'yulab.utils', 'zoo'
    



```R
install.packages("devtools")
```

    还安装依赖关系‘systemfonts’, ‘textshaping’, ‘gert’, ‘ragg’, ‘usethis’, ‘fs’, ‘lifecycle’, ‘pak’, ‘pkgbuild’, ‘pkgdown’, ‘pkgload’, ‘rcmdcheck’, ‘rlang’, ‘roxygen2’, ‘testthat’
    
    
    Warning message in install.packages("devtools"):
    “安装程序包‘gert’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘pak’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘pkgbuild’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘rcmdcheck’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘systemfonts’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘usethis’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘pkgload’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘textshaping’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘roxygen2’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘testthat’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘ragg’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘pkgdown’时退出狀態的值不是0”
    Warning message in install.packages("devtools"):
    “安装程序包‘devtools’时退出狀態的值不是0”
    更新'.Library'里的HTML程序包列表
    
    Making 'packages.html' ...
     做完了。
    



```R
library(devtools)
```


    Error in library(devtools): 不存在叫‘devtools’这个名称的程序包
    Traceback:


    1. stop(packageNotFoundError(package, lib.loc, sys.call()))



```R
# 安装 RcppML（如果未安装）
if (!requireNamespace("RcppML", quietly = TRUE))
    install.packages("RcppML")
library(RcppML)
```

    Warning message in download.file(url, destfile, method, mode = "wb", ...):
    “URL 'https://cran.r-project.org/src/contrib/RcppML_0.3.7.1.tar.gz': Timeout of 60 seconds was reached”


    Error in download.file(url, destfile, method, mode = "wb", ...) : 
      无法打开URL'https://cran.r-project.org/src/contrib/RcppML_0.3.7.1.tar.gz'


    Warning message in download.packages(pkgs, destdir = tmpd, available = available, :
    “下载程序包‘RcppML’时出了问题”



    Error in library(RcppML): 不存在叫‘RcppML’这个名称的程序包
    Traceback:


    1. stop(packageNotFoundError(package, lib.loc, sys.call()))



```R

```
