# ! 准备----
## 清除当前环境中的所有对象
rm(list = ls())
## 设置主文件夹路径, 并设置工作目录
(root_dir <- dirname(rstudioapi::getSourceEditorContext()$path))
# root_dir <- "/dev/sda1/home/tenney/RStudio/"
shelf_dir <- function(dir_name, parent_dir = ".") {
    target_dir <- file.path(parent_dir, dir_name)
    if (!file.exists(target_dir)) {
        dir.create(target_dir, recursive = TRUE)
    }
    setwd(target_dir)
    return(target_dir)
}
data_dir <- shelf_dir("data", root_dir)
results_dir <- shelf_dir("results", root_dir)
tcga_gtex_dir <- shelf_dir("tcga_gtex_conbined", data_dir)

## main
library(librarian)
shelf(dplyr, stringr, quiet = TRUE)

# 以KIRC为例
catagory_name_brief <- "KIRC"
catagory_dir <- shelf_dir(catagory_name_brief, tcga_gtex_dir)
expr <- read.csv(file = file.path(catagory_dir, list.files(catagory_dir, pattern = "expr.csv")), row.names = 1, check.names = FALSE)
dim(expr)
ph <- read.csv(file = file.path(catagory_dir, list.files(catagory_dir, pattern = "ph.csv")), row.names = 1, check.names = FALSE)
## log逆转
expr <- 2^expr - 1
expr[1:4, 1:4]
# tumor和normal统计
library(stringr)
k1 <- str_starts(rownames(ph), "TCGA")
k2 <- as.numeric(str_sub(rownames(ph), 14, 15)) < 10
k_norm_tcga <- as.numeric(str_sub(rownames(ph), 14, 15)) > 10
k_norm_gtex <- str_starts(rownames(ph), "GTEX")
table(k1 & k2)
table(k1 & k_norm_tcga)
table(k_norm_gtex)
group_list <- ifelse(k1 & k2, "tumor", "normal")
group_list <- factor(group_list, levels = c("normal", "tumor"))
table(group_list)
expr[1:4, 1:4]
k3 <- apply(expr, 1, function(x) {
    sum(x > 1)
}) > 150
table(k3)
expr <- expr[k3, ]
### 2.DESeq2---------
expr <- floor(expr)
expr[1:4, 1:4]
library(DESeq2)
colData <- data.frame(
    row.names = colnames(expr),
    condition = group_list
)
dds <- DESeqDataSetFromMatrix(
    countData = expr,
    colData = colData,
    design = ~condition
)
# 参考因子应该是对照组 dds$condition <- relevel(dds$condition, ref = "untrt")

dds <- DESeq(dds)
# 两两比较
res <- results(dds, contrast = c("condition", rev(levels(group_list))))
resOrdered <- res[order(res$pvalue), ] # 按照P值排序
DEG <- as.data.frame(resOrdered)
head(DEG)
# 去除NA值
DEG <- na.omit(DEG)

# 添加change列标记基因上调下调
logFC_cutoff <- with(DEG, mean(abs(log2FoldChange)) + 2 * sd(abs(log2FoldChange)))
# logFC_cutoff <- 1
DEG$change <- as.factor(
    ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
        ifelse(DEG$log2FoldChange > logFC_cutoff, "UP", "DOWN"), "NOT"
    )
)
head(DEG)
table(DEG$change)
DESeq2_DEG <- DEG
### 3.edgeR---------
expr <- expr
library(edgeR)

dge <- DGEList(counts = expr, group = group_list)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + group_list)
rownames(design) <- colnames(dge)
colnames(design) <- levels(group_list)

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast = c(-1, 1))

DEG <- topTags(fit2, n = nrow(exp))
DEG <- as.data.frame(DEG)
logFC_cutoff <- with(DEG, mean(abs(logFC)) + 2 * sd(abs(logFC)))
# logFC_cutoff <- 2
DEG$change <- as.factor(
    ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,
        ifelse(DEG$logFC > logFC_cutoff, "UP", "DOWN"), "NOT"
    )
)
head(DEG)
table(DEG$change)
edgeR_DEG <- DEG
### limma----
library(limma)

design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp)

dge <- DGEList(counts = expr)
dge <- calcNormFactors(dge)

v <- voom(dge, design, normalize = "quantile")
fit <- lmFit(v, design)

constrasts <- paste(rev(levels(group_list)), collapse = "-")
cont.matrix <- makeContrasts(contrasts = constrasts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

DEG <- topTable(fit2, coef = constrasts, n = Inf)
DEG <- na.omit(DEG)
logFC_cutoff <- with(DEG, mean(abs(logFC)) + 2 * sd(abs(logFC)))
# logFC_cutoff <- 2
DEG$change <- as.factor(
    ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
        ifelse(DEG$logFC > logFC_cutoff, "UP", "DOWN"), "NOT"
    )
)
head(DEG)
limma_voom_DEG <- DEG
tj <- data.frame(
    deseq2 = as.integer(table(DESeq2_DEG$change)),
    edgeR = as.integer(table(edgeR_DEG$change)),
    limma_voom = as.integer(table(limma_voom_DEG$change)),
    row.names = c("down", "not", "up")
)
tj
results_dir <- shelf_dir("results", root_dir)
save(DESeq2_DEG, edgeR_DEG, limma_voom_DEG, group_list, tj, file = "DEG.Rdata")
write.csv(DESeq2_DEG, file = file.path(results_dir, "DESeq2_DEG.csv"))
write.csv(edgeR_DEG, file = file.path(results_dir, "edgeR_DEG.csv"))
write.csv(limma_voom_DEG, file = file.path(results_dir, "limma_voom_DEG.csv"))
write.csv(tj, file = file.path(results_dir, "tj.csv"))

load("DEG.Rdata")
list(
    row.names(DESeq2_DEG)[DESeq2_DEG$change == "DOWN"],
    edgeR_DEG$change == "DOWN",
    limma_voom_DEG$change == "DOWN",
    DESeq2_DEG$change == "UP",
    edgeR_DEG$change == "UP",
    limma_voom_DEG$change == "UP"
)
venn_list <- list(
    row.names(DESeq2_DEG)[DESeq2_DEG$change == "DOWN"],
    row.names(edgeR_DEG)[edgeR_DEG$change == "DOWN"],
    row.names(limma_voom_DEG)[limma_voom_DEG$change == "DOWN"],
    row.names(DESeq2_DEG)[DESeq2_DEG$change == "UP"],
    row.names(edgeR_DEG)[edgeR_DEG$change == "UP"],
    row.names(limma_voom_DEG)[limma_voom_DEG$change == "UP"]
)
venn_list <- setNames(venn_list, c("deseq2Down", "edgeRDown", "limma_voomDown", "deseq2Up", "edgeRUp", "limma_voomUp"))
names(venn_list)
shelf(venn, VennDiagram)
results_dir <- shelf_dir("results", root_dir)
pdf("venn.pdf", width = 12, height = 12)
venn(venn_list[1:3],
    zcolor = "style", # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
    opacity = 0.3, # 调整颜色透明度
    box = FALSE, # 是否添加边框
    ilcs = 4, # 数字大小
    sncs = 1, # 组名字体大小
    ilabels = "counts", # 是否显示数字
)
venn(venn_list[4:6],
    zcolor = "style", # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
    opacity = 0.3, # 调整颜色透明度
    box = FALSE, # 是否添加边框
    ilcs = 4, # 数字大小
    sncs = 1, # 组名字体大小
    ilabels = "counts", # 是否显示数字
)
venn(venn_list,
    zcolor = "style", # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
    opacity = 0.3, # 调整颜色透明度
    box = FALSE, # 是否添加边框
    ilcs = 3, # 数字大小
    sncs = 1, # 组名字体大小
    ilabels = "counts", # 是否显示数字
)
dev.off()
