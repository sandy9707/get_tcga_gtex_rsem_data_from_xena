---
title: TCGA 联合 GTEX 分析流程
tags:
-
categories:
- 生物信息学
- TCGA
date: Sep 23, 2023 at 16:23:01
author: yeyezi
cover: /img/叶子的神奇小软件.jpg
---

## 引言

在我们使用 TCGA 数据进行差异分析的时候，可能会遇到肿瘤和正常数据之间的不平衡的问题。为了得到更科学的结果，结合 GTEx 项目提供的正常样本数据是一种较为常见的解决方法，因此本篇文章整理和运行了一种可行的数据下载和处理流程。

>高通量 RNA 测序（RNA-Seq）已成为转录组分析的强大方法，广泛用于了解基因功能和生物模式，找到候选药物靶点，并识别疾病分类和诊断的生物标志物。近年来，癌症基因组图谱（TCGA）和基因型组织表达（GTEx）项目为数万个癌症和非癌症样本提供了 RNA-Seq 数据，为包括癌症生物学在内的许多相关领域提供了前所未有的机会。到目前为止，TCGA 已经为 33 种癌症类型的 9736 个肿瘤样本提供了 RNA-Seq 数据，此外还有 726 个相邻正常组织的数据。肿瘤和正常数据之间的不平衡可能导致各种差异分析的效率低下。幸运的是，GTEx 项目为 8000 多个正常样本提供了 RNA-Seq 数据，尽管这些样本来自不相关的捐赠者。由于数据处理管道和基因模型等方面的许多差异，此类数据无法直接组合进行综合分析。为了使来自不同来源的数据更加兼容，UCSC Xena 项目（<http://xena.ucsc.edu/>）基于标准管道重新计算了所有表达式原始数据，以尽量减少与不同来源的差异，从而允许形成最新的最全面的表达式数据。[1]

## 效果展示

所有代码和需要的文件已提交至 github 仓库，仓库地址：
<https://github.com/sandy9707/get_tcga_gtex_rsem_data_from_xena>

```txt
├── data
│   ├── GTEX
│   │   ├── Adipose Tissue
│   │   │   ├── Adipose Tissue_expr.csv
│   │   │   ├── Adipose Tissue_expr_ph.Rdata
│   │   │   └── Adipose Tissue_ph.csv
│   │   ├── ...
│   ├── TCGA
│   │   ├── Acute Myeloid Leukemia
│   │   │   ├── Acute Myeloid Leukemia_expr.csv
│   │   │   ├── Acute Myeloid Leukemia_expr_ph.Rdata
│   │   │   └── Acute Myeloid Leukemia_ph.csv
│   │   ├── ...
│   ├── tcga_gtex_conbined
│   │   ├── ACC
│   │   │   ├── Adrenocortical Cancer_expr.csv
│   │   │   ├── Adrenocortical Cancer_expr_ph.Rdata
│   │   │   └── Adrenocortical Cancer_ph.csv
│   │   ├── ...
├── reference
│   ├── detailed_category.csv
│   ├── TCGA_catagory.csv
│   └── tcga_gtex_catagory.csv
├── results
│   ├── DEG.Rdata
│   ├── DESeq2_DEG.csv
│   ├── edgeR_DEG.csv
│   ├── limma_voom_DEG.csv
│   ├── tj.csv
│   └── venn.pdf
```

## 过程

### 获取 xena ucsc 数据

数据可以直接从<https://xenabrowser.net/datapages/>网站获得。

获得合并后数据或者分别获得 TCGA 和 GTEx 数据均可，为确保可以比较，本文选择直接下载合并数据[TCGA TARGET GTEx (13 datasets)](https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). 网页提供了 gene expression RNAseq 的 RSEM expected_count (DESeq2 standardized), 可以直接用于差异分析。

### 处理

具体观看<https://github.com/sandy9707/get_tcga_gtex_rsem_data_from_xena>的`README.md`文件。

大致步骤分为：

1. 分离为 GTEX 和 TCGA 数据集
2. 对于两个数据集中的每个主要部位，保存相应的表达和表型数据到单独的文件中。
3. 根据预设的分类合并来自 TCGA 和 GTEX 数据集的相关样本，并将其表达数据和临床数据保存以供进一步分析。

实际上就是利用 tcga 项目和患病部位的对应关系添加上相应部位的 gtex 数据，需要注意的就是在 gtex 中实际部位和标注部位不一定一致，需要自行配置对应文件，本篇中使用的配置文件均放入`reference`文件夹。

### 篇外

#### 插曲

对于`TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz`文件，无法解压和读入，原因未知。解压命令为`gzip -d`.

无论使用`readr`的`read.delim`函数，还是`data.table`的`fread`函数，都无法读取成功。读入命令和显示信息如下：

```R
## read.table
r$> data <- read.table(file.path(data_folder, "TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz"), sep = "\t")
# Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#   line 11069 did not have 19040 elements

## read.delim
data <- readr::read.delim(file.path(data_folder, "TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz"), header = TRUE, as.is = TRUE)
data <- data.table::fread(file.path(data_folder, "TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz"), data.table = FALSE)
```

### Brain 组织数量与 GEPIA 不一致的处理方法

在经过 UCSC 处理 GTEx 中，Brain 分为 13 个组织 1148 个样本，但 GEPIA2 将他们分成了由 207+945=1152 两组，那么哪些是可以作为 LGG 或 GBM 的对照的呢？

从 TCGA 中下载 LGG 的临床数据，可得都来源于额叶、顶叶、颞叶和枕叶四个脑叶，也就是说都来源于大脑 (Cerebrum), 因此我们要找的也是 Cerebrum.

GTEx Brain 的 13 个组织：

```txt
Brain - Amygdala
脑杏仁核（Amygdala）69 正常组织
Brain - Anterior Cingulate Cortex (Ba24)
大脑前扣带皮质（Ba24）83 正常组织
Brain - Caudate (Basal Ganglia)
脑尾（基底神经节）（Caudate）108 正常组织
Brain - Cerebellar Hemisphere
小脑半球（Cerebellar Hemisphere）97 正常组织
Brain - Cerebellum
小脑（Cerebellum）117 正常组织
Brain - Cortex
大脑皮质（Cortex）105 正常组织
Brain - Frontal Cortex (Ba9)
大脑额叶皮质（Ba9）101 正常组织
Brain - Hippocampus
大脑海马（Hippocampus）84 正常组织
Brain - Hypothalamus
大脑下丘脑（Hypothalamus）82 正常组织
Brain - Nucleus Accumbens (Basal Ganglia)
脑积聚核（基底神经节）（Nucleus Accumbens）104 正常组织
Brain - Putamen (Basal Ganglia)
脑壳 Putamen（基底神经节）81 正常组织
Brain - Spinal Cord (Cervical C-1)
脊髓（颈 C-1）（Spinal Cord）60 正常组织
Brain - Substantia Nigra
大脑黑质（Substantia Nigra）57 正常组织
```

首先，小脑 (Cerebellum) 和小脑半球 (Cerebellar Hemisphere) 显然属于小脑半球，所以 214 例属于 945 组。

而后，根据 [基底核 - 维基百科，自由的百科全书](https://zh.wikipedia.org/wiki/基底核) 的说明，
基底核（basal ganglia）包括尾/壳/黑质，共 393 例，因此属于 945 组。

最后，海马，下丘脑，脊髓，杏仁核，前扣带都属于单独的结构，所以推测只有皮质（Cortex）和额叶皮质（Ba9）属于大脑 (Cerebrum), 样本量为 206, 在误差范围内，可以被认为是 LGG 或 GBM 的对照组。

基底核（basal ganglia）:

- 前侧
  - 纹状体（Striatum）包括
    - 尾状核（Caudate nucleus）
    - 壳（Putamen）
    - 伏隔核（Nucleus accumbens）
    - 外苍白球（External segment of globus pallidus，GPe）
    - 内苍白球（Internal segment of globus pallidus，GPi）
- 后侧，以下这些结构在大脑中更靠下，靠后。
  - 丘脑下核（Subthalamic nucleus, STN）
  - 黑质（Substantia nigra, SN），根据内部结构可分为
    - 黑质致密部（Substantia nigra pars compacta，SNc）
    - 黑质网状部（Substantia nigra pars reticulata，SNr）
    - 黑质侧部（Substantia nigra pars lateralis，SNl）

## 结论

## 引用

1. [GEPIA：用于癌症和正常基因表达分析和交互式分析的网络服务器](https://academic.oup.com/nar/article/45/W1/W98/3605636?login=false)
2. [3 大数据库超 2 万 RNA-seq 数据重新统一处理——关于 TCGA-GTEx 是否需要标准化 – 王进的个人网站](https://www.jingege.wang/2023/05/24/3大数据库超2万rna-seq数据重新统一处理-关于tcga-gtex是否/)
3. [GTEx 联合 TCGA 数据库差异分析（更新） – 王进的个人网站](https://www.jingege.wang/2022/03/16/gtex联合tcga数据库差异分析/)
4. [TCGA 和 GTEx 的数据联合分析实战 - 简书](https://www.jianshu.com/p/46b048220b88)
5. [GitHub - xjsun1221/RSEM_with_limma_edgeR_Deseq2](https://github.com/xjsun1221/RSEM_with_limma_edgeR_Deseq2)
6. [GEPIA 2 - Dataset Sources](http://gepia2.cancer-pku.cn/#dataset)
