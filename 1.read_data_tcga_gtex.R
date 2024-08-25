# ! 准备----
## 清除当前环境中的所有对象
rm(list = ls())
## 设置主文件夹路径, 并设置工作目录
(root_dir <- dirname(rstudioapi::getSourceEditorContext()$path))
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
ref_dir <- shelf_dir("reference", root_dir)

# !main----
library(librarian)
shelf(dplyr, stringr, quiet = TRUE)
data_dir <- shelf_dir("data", root_dir)
## 第一次运行需要下载文件并仅输出tcga和gtex的表达量和临床信息

if (!("gtex_tcga_expr_ph.Rdata" %in% list.files())) {
    # 输入"TcgaTargetGtex_gene_expected_count.gz"的真实位置
    TcgaTargetGtex_gene_expected_count_file <- "TcgaTargetGtex_gene_expected_count.gz"
    counts <- data.table::fread(TcgaTargetGtex_gene_expected_count_file, data.table = FALSE)
    ph <- data.table::fread("TcgaTargetGTEX_phenotype.txt.gz", data.table = FALSE)
    counts[1:4, 1:4]
    rownames(counts) <- counts$sample
    counts <- counts[, -1]
    rownames(ph) <- ph$sample
    ph <- ph[, -1]
    expr <- as.matrix(counts)
    dim(expr)
    dim(ph)
    # 匹配含临床信息样本名, 结果显示22例样本没有匹配到表达量信息, 有19109例成功匹配
    table(rownames(ph) %in% colnames(expr))
    print(ph[!rownames(ph) %in% colnames(expr), ])
    # 匹配表达量信息, 结果19109例全部成功匹配
    table(colnames(expr) %in% rownames(ph))
    ph <- ph[match(colnames(expr), rownames(ph)), ]
    identical(colnames(expr), rownames(ph))
    # 挑选出TCGA和gtex样本
    table(ph$`_study`)
    keep <- ph$`_study` != "TARGET"
    ph <- ph[keep, ]
    expr <- expr[, keep]
    save(expr, ph, file = "gtex_tcga_expr_ph.Rdata")
    write.csv(ph, file = "gtex_tcga_ph.csv")
}

## 挑出不同数据库数据
load("gtex_tcga_expr_ph.Rdata")
head(ph)
colnames(ph)

## GTEX
ph_gtex <- ph[ph$`_study` == "GTEX", ]
gtex_dir <- shelf_dir("GTEX", data_dir)
# The detailed_category is organized by the lev.
detailed_category <- read.csv(file.path(ref_dir, "detailed_category.csv"))
colnames(detailed_category)
primary_sites <- unique(detailed_category$`X_primary_site`)
# primary_sites <- primary_sites[5:8]
primary_sites
for (primary_site in primary_sites) {
    # primary_site <- primary_sites[1]
    detailed_category_site <- detailed_category[detailed_category$`X_primary_site` == primary_site, "detailed_category"]
    keep2 <- ph$`detailed_category` %in% detailed_category_site
    selected_expr <- expr[, keep2]
    selected_ph <- ph[keep2, ]
    # nrow(selected_ph)
    primary_site_dir <- shelf_dir(primary_site, gtex_dir)
    save(selected_ph, selected_expr, file = file.path(primary_site_dir, paste0(primary_site, "_expr_ph.Rdata")))
    write.csv(selected_expr, file = file.path(primary_site_dir, paste0(primary_site, "_expr.csv")))
    write.csv(selected_ph, file = file.path(primary_site_dir, paste0(primary_site, "_ph.csv")))
}

## TCGA
tcga_dir <- shelf_dir("TCGA", data_dir)
colnames(ph)
ph_TCGA <- ph[ph$`_study` == "TCGA", ]
table(ph_TCGA$`detailed_category`) %>% length()
### get TCGA catagory
# 第一次运行填TRUE
if (FALSE) {
    TCGA_catagory <- table(ph_TCGA$`detailed_category`) %>% data.frame()
    i <- 1
    # detailed_category为空是Control Analyte, 可去除
    # ph_TCGA[ph_TCGA$`detailed_category` == TCGA_catagory[i, 1], ]
    TCGA_catagory <- TCGA_catagory[-1, ]
    TCGA_catagory_sites <- list()
    for (i in 1:nrow(TCGA_catagory)) {
        TCGA_catagory_sites[[i]] <- ph_TCGA[ph_TCGA$`detailed_category` == TCGA_catagory[i, 1], "_primary_site"] %>% unique()
    }
    TCGA_catagory$site <- TCGA_catagory_sites %>% unlist()
    # 保存并**手动**修改------
    colnames(TCGA_catagory) <- c("category_ph", "category_num", "site")
    write.csv(TCGA_catagory, file = file.path(ref_dir, "TCGA_catagory.csv"), row.names = FALSE)
}
TCGA_catagory <- read.csv(file = file.path(ref_dir, "TCGA_catagory.csv"))
head(TCGA_catagory)
dim(TCGA_catagory)
head(ph_TCGA)
for (i in 1:nrow(TCGA_catagory)) {
    catagory_name <- TCGA_catagory[i, 1]
    keep3 <- ph$`detailed_category` %in% TCGA_catagory[i, "category_ph"]
    # sum(keep3)
    selected_expr <- expr[, keep3]
    selected_ph <- ph[keep3, ]
    # nrow(selected_ph)
    catagory_dir <- shelf_dir(catagory_name, tcga_dir)
    save(selected_ph, selected_expr, file = file.path(catagory_dir, paste0(catagory_name, "_expr_ph.Rdata")))
    write.csv(selected_expr, file = file.path(catagory_dir, paste0(catagory_name, "_expr.csv")))
    write.csv(selected_ph, file = file.path(catagory_dir, paste0(catagory_name, "_ph.csv")))
}

## tcga_gtex
tcga_gtex_dir <- shelf_dir("tcga_gtex_conbined", data_dir)
ph_tcga <- ph[ph$`_study` %in% c("TCGA"), ]
head(ph_tcga)
dim(ph_tcga)
ph_gtex <- ph[ph$`_study` %in% c("GTEX"), ]
head(ph_gtex)
dim(ph_gtex)
tcga_gtex_catagory <- read.csv(file = file.path(ref_dir, "tcga_gtex_catagory.csv"))
detailed_category <- read.csv(file.path(ref_dir, "detailed_category.csv"))
head(tcga_gtex_catagory)
for (i in 1:nrow(tcga_gtex_catagory)) {
    # i <- 4
    catagory_name <- tcga_gtex_catagory[i, "category_ph"]
    selected_tcga_ph <- ph_tcga[ph_tcga$`detailed_category` %in% catagory_name, ]
    # print(nrow(selected_tcga_ph))
    primary_site <- tcga_gtex_catagory[i, "site_custom"]
    # print(primary_site)
    detailed_category_site <- detailed_category[detailed_category$`X_primary_site` == primary_site, "detailed_category"]
    keep_gtex_index <- ph$`detailed_category` %in% detailed_category_site
    selected_gtex_ph <- ph[keep_gtex_index, ]
    selected_ph <- rbind(selected_tcga_ph, selected_gtex_ph)
    # print(nrow(selected_ph))
    selected_expr <- expr[, rownames(selected_ph)]
    catagory_name_brief <- tcga_gtex_catagory[i, "TCGA"]
    catagory_dir <- shelf_dir(catagory_name_brief, tcga_gtex_dir)
    save(selected_ph, selected_expr, file = file.path(catagory_dir, paste0(catagory_name, "_expr_ph.Rdata")))
    write.csv(selected_expr, file = file.path(catagory_dir, paste0(catagory_name, "_expr.csv")))
    write.csv(selected_ph, file = file.path(catagory_dir, paste0(catagory_name, "_ph.csv")))
}
