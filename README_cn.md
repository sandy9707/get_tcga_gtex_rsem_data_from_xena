# TCGA + GTEX 数据处理流程

该仓库包含用于处理和合并 TCGA（癌症基因组图谱）和 GTEX（基因型 - 组织表达）数据集的脚本集。目标是为下游分析准备这些数据集，包括基因表达分析和临床数据整合。

## 功能

- **数据加载与合并**：处理并合并大型 TCGA 和 GTEX 基因表达数据集。
- **自定义数据分割**：自动按研究类型和详细分类分离并组织数据。
- **灵活的目录管理**：自动创建和管理数据及结果存储的目录结构。

## 先决条件

确保你已经安装了以下 R 包：

- `dplyr`
- `stringr`
- `data.table`
- `librarian`（用于管理 R 包依赖）

你可以使用以下 R 命令来安装这些包：

```R
install.packages(c("dplyr", "stringr", "data.table"))
```

你还需要以下文件（可以从 TCGA 和 GTEX 数据库获取）：

- `TcgaTargetGtex_gene_expected_count.gz`：包含基因表达数据的文件。
- `TcgaTargetGTEX_phenotype.txt.gz`：包含表型和临床数据的文件。
- `detailed_category.csv`：将主要部位映射到详细分类的参考文件。
- `TCGA_catagory.csv` 和 `tcga_gtex_catagory.csv`：这些文件按站点和研究对 TCGA 和 TCGA-GTEX 合并数据集进行分类。

## 使用方法

### 1. 初始设置

```R
# 清空环境
rm(list = ls())

# 设置主要工作目录
root_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
data_dir <- shelf_dir("data", root_dir)
results_dir <- shelf_dir("results", root_dir)
ref_dir <- shelf_dir("reference", root_dir)
```

### 2. 数据下载与预处理

脚本会检查是否存在已处理的数据文件（`gtex_tcga_expr_ph.Rdata`）。如果这些文件不存在，它将加载并预处理原始表达和表型数据。

```R
# 如果尚未处理数据，则加载并预处理
if (!("gtex_tcga_expr_ph.Rdata" %in% list.files())) {
    TcgaTargetGtex_gene_expected_count_file <- "TcgaTargetGtex_gene_expected_count.gz"
    counts <- data.table::fread(TcgaTargetGtex_gene_expected_count_file, data.table = FALSE)
    ph <- data.table::fread("TcgaTargetGTEX_phenotype.txt.gz", data.table = FALSE)

    # 在此处执行额外的过滤和处理步骤
    save(expr, ph, file = "gtex_tcga_expr_ph.Rdata")
}
```

### 3. 按数据来源分离数据

预处理后，脚本将数据分离为 GTEX 和 TCGA 数据集，并进一步按主要部位对它们进行分类。

```R
# 分离GTEX数据
ph_gtex <- ph[ph$`_study` == "GTEX", ]
# 分离TCGA数据
ph_TCGA <- ph[ph$`_study` == "TCGA", ]
```

### 4. 进一步数据分割

对于数据集中的每个主要部位，脚本将保存相应的表达和表型数据到单独的文件中。其中，GTEX 数据的主要部位和实际部位标识的对照关系来自于`detailed_category.csv`, TCGA 数据的项目和实际部位标识的对照关系来自于`TCGA_catagory`。

```R
# 按主要部位保存GTEX和TCGA数据
for (primary_site in primary_sites) {
    # 保存GTEX数据
    primary_site_dir <- shelf_dir(primary_site, gtex_dir)
    save(selected_ph, selected_expr, file = file.path(primary_site_dir, paste0(primary_site, "_expr_ph.Rdata")))

    # 保存TCGA数据
    catagory_dir <- shelf_dir(catagory_name, tcga_dir)
    save(selected_ph, selected_expr, file = file.path(catagory_dir, paste0(catagory_name, "_expr_ph.Rdata")))
}
```

### 5. 合并 TCGA 和 GTEX 数据

此外，脚本还根据预设的分类合并来自 TCGA 和 GTEX 数据集的相关样本，并将其保存以供进一步分析。

GTEX 和 TCGA 数据的对照关系来自于`tcga_gtex_catagory.csv`. 另外，还包含了癌症样本、癌旁样本、正常样本的实际数量。

## 文件结构

- **data/**：原始和处理后的数据文件。
- **results/**：数据处理结果。
- **reference/**：参考文件，如 TCGA 和 GTEX 数据集的分类映射。

## 许可证

此项目采用 MIT 许可证。
