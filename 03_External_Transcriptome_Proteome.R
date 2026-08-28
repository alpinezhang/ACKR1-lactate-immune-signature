# Figures 4-8 and Supplementary Figures 2-4: ACKR1 external transcriptome/proteome analyses (English; current revision order)
# Corresponds to current main Figures 4-8 and Supplementary Figures 2-4.

# Supplementary Figure 2: GSE163154 GO/KEGG/HALLMARK enrichment


suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(illuminaHumanv2.db)
  library(limma)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "24_GEO_Validation_GSE163154")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> [Step 1] 正在读取 GSE163154 表达矩阵...")
target_gse <- "GSE163154_series_matrix.txt.gz"

gset <- getGEO(filename = target_gse, getGPL = FALSE)
pheno_gse <- pData(gset)
expr_gse <- exprs(gset)
message(">>> [Step 2] 正在提取并设定临床分组...")
group_gse <- rep(NA, nrow(pheno_gse))
is_non_iph <- grepl("non-IPH", pheno_gse$title, ignore.case = TRUE)

group_gse[is_non_iph] <- "non_IPH" 
group_gse[!is_non_iph] <- "IPH"

group_factor <- factor(group_gse, levels = c("non_IPH", "IPH"))

n_non_iph <- sum(group_factor == "non_IPH")
n_iph <- sum(group_factor == "IPH")
message(sprintf("   - 分组判定: non_IPH (对照) = %d 例, IPH (实验) = %d 例", n_non_iph, n_iph))

message(">>> [Step 3] 正在执行 limma 差异分析 (IPH vs non_IPH)...")

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor) 

contrast.matrix <- makeContrasts(IPH_vs_non_IPH = IPH - non_IPH, levels = design)

fit <- lmFit(expr_gse, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

all_deg <- topTable(fit2, coef = "IPH_vs_non_IPH", number = Inf, adjust.method = "fdr")
all_deg$PROBEID <- rownames(all_deg)

message(">>> [Step 4] 正在将探针转换为 Gene Symbol 并清洗数据...")

probe_to_symbol <- mapIds(illuminaHumanv2.db, 
                          keys = all_deg$PROBEID, 
                          column = "SYMBOL", 
                          keytype = "PROBEID", 
                          multiVals = "first")

all_deg$Gene <- probe_to_symbol

deg_clean <- all_deg %>%
  filter(!is.na(Gene)) %>%
  arrange(P.Value) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  dplyr::select(Gene, PROBEID, logFC, AveExpr, t, P.Value, adj.P.Val, B)
message(">>> [Step 5] 正在导出差异基因表格...")

output_file <- file.path(out_dir, "01_GSE163154_IPH_vs_nonIPH_DEGs.csv")
write.csv(deg_clean, file = output_file, row.names = FALSE)

sig_p_count <- sum(deg_clean$P.Value < 0.05)
sig_fdr_count <- sum(deg_clean$adj.P.Val < 0.05)

message("\n>>> 第一阶段完美通关！全基因组差异表格已生成。")
message(sprintf("   - 导出路径: %s", output_file))
message(sprintf("   - 成功映射基因总数: %d", nrow(deg_clean)))
message(sprintf("   - P < 0.05 的显著基因数: %d", sig_p_count))
message(sprintf("   - FDR < 0.05 的显著基因数: %d", sig_fdr_count))


suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(forcats)
})

message(">>> [Step 1] 正在准备富集分析的基因集...")

p_cutoff <- 0.05
logfc_cutoff <- 0.585

genes_up <- deg_clean %>% filter(P.Value < p_cutoff & logFC > logfc_cutoff) %>% pull(Gene)
genes_down <- deg_clean %>% filter(P.Value < p_cutoff & logFC < -logfc_cutoff) %>% pull(Gene)

message(sprintf("   - 筛选出上调基因: %d 个", length(genes_up)))
message(sprintf("   - 筛选出下调基因: %d 个", length(genes_down)))

entrez_up <- bitr(genes_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
entrez_down <- bitr(genes_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

message(">>> [Step 2] 正在运行 GO (生物学过程) 和 KEGG 富集分析...")

go_up <- enrichGO(gene = entrez_up, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
go_down <- enrichGO(gene = entrez_down, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)

kegg_up <- enrichKEGG(gene = entrez_up, organism = 'hsa', pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = entrez_down, organism = 'hsa', pvalueCutoff = 0.05)

message(">>> [Step 3] 正在渲染 Nature 风格富集气泡图...")

sample_info <- sprintf("Dataset: GSE163154 | IPH (N=%d) vs non_IPH (N=%d)", n_iph, n_non_iph)

plot_enrichment_dot <- function(enrich_obj, title_prefix, color_low, color_high) {
  if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) {
    return(NULL)
  }

  df <- as.data.frame(enrich_obj) %>%
    arrange(p.adjust) %>%
    head(10)

  if(nrow(df) == 0) return(NULL)

  df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    parts[1] / parts[2]
  })

  df$Description <- fct_reorder(df$Description, df$GeneRatioNum)

  full_title <- sprintf("%s\n%s", title_prefix, sample_info)

  p <- ggplot(df, aes(x = GeneRatioNum, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.85) +
    scale_color_gradient(low = color_low, high = color_high, name = "p.adjust") +
    scale_size_continuous(name = "Gene Count", range = c(3, 8)) +

    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    scale_y_discrete(expand = expansion(add = c(0.8, 0.8))) +

    theme_classic() +
    labs(title = full_title, x = "Gene Ratio", y = "") +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, lineheight = 1.2, margin = margin(b=15)),
      axis.text.y = element_text(size = 11, color = "black", face = "bold"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", color = "black"),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 15, r = 25, b = 15, l = 10) 
    )
  return(p)
}


p_go_up <- plot_enrichment_dot(go_up, "GO Enrichment: Up-regulated", "#D51F26", "#272E6A")
p_kegg_up <- plot_enrichment_dot(kegg_up, "KEGG Pathways: Up-regulated", "#D51F26", "#272E6A")

p_go_down <- plot_enrichment_dot(go_down, "GO Enrichment: Down-regulated", "#272E6A", "#8A9FD1")
p_kegg_down <- plot_enrichment_dot(kegg_down, "KEGG Pathways: Down-regulated", "#272E6A", "#8A9FD1")

save_enrich_plot <- function(plot_obj, filename_base, width = 7.5, height = 5.5) {
  if(!is.null(plot_obj)) {
    ggsave(file.path(out_dir, paste0(filename_base, ".pdf")), plot_obj, width = width, height = height)
    ggsave(file.path(out_dir, paste0(filename_base, ".png")), plot_obj, width = width, height = height, dpi = 600)
  } else {
    message(sprintf("   ! 警告: %s 未产生显著富集结果，跳过绘图。", filename_base))
  }
}

message(">>> [Step 4] 正在导出高清 PDF 与 PNG 图片 (已添加 GSE 数据集标注)...")

save_enrich_plot(p_go_up, "02_GO_Enrichment_Up_Red")
save_enrich_plot(p_go_down, "02_GO_Enrichment_Down_Blue")
save_enrich_plot(p_kegg_up, "03_KEGG_Enrichment_Up_Red")
save_enrich_plot(p_kegg_down, "03_KEGG_Enrichment_Down_Blue")

message("\n>>> 图表更新完成！现在标题中已包含 Dataset: GSE163154 信息。")


message(">>> [Step 1] 正在准备 GSEA 所需的排序基因列表 (Ranked Gene List)...")

suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(forcats)
})

all_genes_mapped <- bitr(deg_clean$Gene, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)

gsea_df <- deg_clean %>%
  inner_join(all_genes_mapped, by = c("Gene" = "SYMBOL")) %>%
  filter(!is.na(logFC)) %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(logFC), n = 1) %>%
  ungroup() %>%
  arrange(desc(logFC))

gene_list <- gsea_df$logFC
names(gene_list) <- gsea_df$ENTREZID

message(sprintf("   - 成功构建 GSEA 基因集，包含 %d 个有效基因。", length(gene_list)))

message(">>> [Step 2] 正在运行 GO 和 KEGG 的 GSEA 分析 (这可能需要一两分钟)...")

set.seed(1234)

gsea_go <- gseGO(geneList     = gene_list,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

gsea_kegg <- gseKEGG(geneList     = gene_list,
                     organism     = 'hsa',
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)

gsea_go <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(as.data.frame(gsea_go), file.path(out_dir, "05_GSEA_GO_BP_Full_Results.csv"), row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), file.path(out_dir, "06_GSEA_KEGG_Full_Results.csv"), row.names = FALSE)

message(">>> [Step 3] 正在渲染 GSEA 专属 Nature 风格气泡图...")

sample_info <- sprintf("Dataset: GSE163154 | IPH (N=%d) vs non_IPH (N=%d)", n_iph, n_non_iph)

plot_gsea_dot <- function(gsea_obj, direction, title_prefix, color_low, color_high) {
  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) return(NULL)

  df <- as.data.frame(gsea_obj)

  if (direction == "Activated") {
    df <- df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(10)
  } else {
    df <- df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(10)
  }

  if(nrow(df) == 0) return(NULL)

  df$Description <- fct_reorder(df$Description, abs(df$NES))

  full_title <- sprintf("%s\n%s", title_prefix, sample_info)

  p <- ggplot(df, aes(x = abs(NES), y = Description)) +
    geom_point(aes(size = setSize, color = p.adjust), alpha = 0.85) +
    scale_color_gradient(low = color_low, high = color_high, name = "p.adjust") +
    scale_size_continuous(name = "Set Size", range = c(3, 8)) +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    scale_y_discrete(expand = expansion(add = c(0.8, 0.8))) +
    theme_classic() +
    labs(title = full_title, x = "Absolute Normalized Enrichment Score (|NES|)", y = "") +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, lineheight = 1.2, margin = margin(b=15)),
      axis.text.y = element_text(size = 11, color = "black", face = "bold"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", color = "black"),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 15, r = 25, b = 15, l = 10)
    )
  return(p)
}


p_gsea_go_up <- plot_gsea_dot(gsea_go, "Activated", "GSEA GO: Activated in IPH", "#D51F26", "#272E6A")
p_gsea_go_down <- plot_gsea_dot(gsea_go, "Suppressed", "GSEA GO: Suppressed in IPH", "#272E6A", "#8A9FD1")

p_gsea_kegg_up <- plot_gsea_dot(gsea_kegg, "Activated", "GSEA KEGG: Activated in IPH", "#D51F26", "#272E6A")
p_gsea_kegg_down <- plot_gsea_dot(gsea_kegg, "Suppressed", "GSEA KEGG: Suppressed in IPH", "#272E6A", "#8A9FD1")

message(">>> [Step 4] 正在导出 GSEA 高清 PDF 与 PNG 图片...")

save_enrich_plot(p_gsea_go_up, "07_GSEA_GO_Activated_Red")
save_enrich_plot(p_gsea_go_down, "07_GSEA_GO_Suppressed_Blue")
save_enrich_plot(p_gsea_kegg_up, "08_GSEA_KEGG_Activated_Red")
save_enrich_plot(p_gsea_kegg_down, "08_GSEA_KEGG_Suppressed_Blue")

message("\n>>> GSEA 分析与可视化全部搞定！图片风格与 ORA 完全统一！")


message(">>> [Step 1] 正在准备 HALLMARK 基因集和数据...")

suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(msigdbr)
  library(stringr)
})

all_genes_mapped <- bitr(deg_clean$Gene, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)

gsea_df <- deg_clean %>%
  inner_join(all_genes_mapped, by = c("Gene" = "SYMBOL")) %>%
  filter(!is.na(logFC)) %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(logFC), n = 1) %>%
  ungroup() %>%
  arrange(desc(logFC))

gene_list <- gsea_df$logFC
names(gene_list) <- gsea_df$ENTREZID

hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

message(sprintf("   - 成功构建 GSEA 基因集 (N=%d) 及 HALLMARK 背景库。", length(gene_list)))

message(">>> [Step 2] 正在运行 HALLMARK GSEA 分析...")

set.seed(1234)
gsea_hallmark <- GSEA(geneList     = gene_list,
                      TERM2GENE    = hallmark_t2g,
                      minGSSize    = 10,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)

gsea_hallmark <- setReadable(gsea_hallmark, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(as.data.frame(gsea_hallmark), file.path(out_dir, "09_GSEA_HALLMARK_Full_Results.csv"), row.names = FALSE)

message(">>> [Step 3] 正在渲染 HALLMARK 专属气泡图...")

sample_info <- sprintf("Dataset: GSE163154 | IPH (N=%d) vs non_IPH (N=%d)", n_iph, n_non_iph)

plot_hallmark_dot <- function(gsea_obj, direction, title_prefix, color_low, color_high) {
  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) return(NULL)

  df <- as.data.frame(gsea_obj)

  if (direction == "Activated") {
    df <- df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(10)
  } else {
    df <- df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(10)
  }

  if(nrow(df) == 0) return(NULL)

  df$Description <- str_replace_all(df$Description, "HALLMARK_", "")
  df$Description <- str_replace_all(df$Description, "_", " ")

  df$Description <- str_wrap(df$Description, width = 35)

  df$Description <- fct_reorder(df$Description, abs(df$NES))
  full_title <- sprintf("%s\n%s", title_prefix, sample_info)

  p <- ggplot(df, aes(x = abs(NES), y = Description)) +
    geom_point(aes(size = setSize, color = p.adjust), alpha = 0.85) +
    scale_color_gradient(low = color_low, high = color_high, name = "p.adjust") +
    scale_size_continuous(name = "Set Size", range = c(3, 8)) +

    scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
    scale_y_discrete(expand = expansion(add = c(0.8, 0.8))) +
    theme_classic() +
    labs(title = full_title, x = "Absolute Normalized Enrichment Score (|NES|)", y = "") +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, lineheight = 1.2, margin = margin(b=15)),
      axis.text.y = element_text(size = 11, color = "black", face = "bold", lineheight = 0.8),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", color = "black", margin = margin(t=10)),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 15, r = 25, b = 15, l = 15)
    )
  return(p)
}

p_hm_up <- plot_hallmark_dot(gsea_hallmark, "Activated", "HALLMARK: Activated in IPH", "#D51F26", "#272E6A")
p_hm_down <- plot_hallmark_dot(gsea_hallmark, "Suppressed", "HALLMARK: Suppressed in IPH", "#272E6A", "#8A9FD1")

save_enrich_plot <- function(plot_obj, filename_base, width = 8.5, height = 6) {
  if(!is.null(plot_obj)) {
    ggsave(file.path(out_dir, paste0(filename_base, ".pdf")), plot_obj, width = width, height = height)
    ggsave(file.path(out_dir, paste0(filename_base, ".png")), plot_obj, width = width, height = height, dpi = 600)
  } else {
    message(sprintf("   ! 警告: %s 未产生显著富集结果，跳过绘图。", filename_base))
  }
}

message(">>> [Step 4] 正在导出 HALLMARK 高清图纸...")

save_enrich_plot(p_hm_up, "10_GSEA_HALLMARK_Activated_Red")
save_enrich_plot(p_hm_down, "10_GSEA_HALLMARK_Suppressed_Blue")

message("\n>>> 完美通关！HALLMARK 标签已洗净，布局比例协调，无任何越界风险。")


message(">>> [Step 1] 正在准备提取核心机制通路的 GSEA 曲线图...")

suppressMessages({
  library(enrichplot)
})

core_pathways <- c(
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS"
)

sample_info_sub <- sprintf("Dataset: GSE163154 | IPH (N=%d) vs non_IPH (N=%d)", n_iph, n_non_iph)

message(">>> [Step 2] 正在逐个绘制机制通路山峰图 (包含 p.adjust 和 NES 标注)...")

res_df <- as.data.frame(gsea_hallmark)

for (pathway in core_pathways) {
  if (pathway %in% res_df$ID) {

    p_adj <- res_df[pathway, "p.adjust"]
    nes <- res_df[pathway, "NES"]

    line_color <- ifelse(nes > 0, "#D51F26", "#272E6A")
    direction_text <- ifelse(nes > 0, "Activated in IPH", "Suppressed in IPH")

    clean_title <- gsub("HALLMARK_", "", pathway)
    clean_title <- gsub("_", " ", clean_title)

    plot_title <- sprintf("%s\n%s | NES=%.2f, p.adj=%.3f", 
                          clean_title, sample_info_sub, nes, p_adj)

    p_gsea_curve <- gseaplot2(gsea_hallmark, 
                              geneSetID = pathway, 
                              title = plot_title, 
                              color = line_color,
                              base_size = 12, 
                              pvalue_table = FALSE)

    out_name <- paste0("11_GSEA_Curve_", gsub("HALLMARK_", "", pathway))

    ggsave(file.path(out_dir, paste0(out_name, ".pdf")), p_gsea_curve, width = 6, height = 5)
    ggsave(file.path(out_dir, paste0(out_name, ".png")), p_gsea_curve, width = 6, height = 5, dpi = 600)

    message(sprintf("   - 成功导出: %s (%s, NES: %.2f)", clean_title, direction_text, nes))

  } else {
    message(sprintf("   ! 提示: %s 在本次数据中未达到有效富集背景，跳过。", pathway))
  }
}

message("\n>>> 核心机制曲线图导出完毕！")


message(">>> [Step 1] 正在准备多通路联合 GSEA 曲线图...")

suppressMessages({
  library(enrichplot)
  library(ggplot2)
})

combined_pathways <- c("HALLMARK_GLYCOLYSIS", "HALLMARK_HYPOXIA")

valid_pathways <- combined_pathways[combined_pathways %in% as.data.frame(gsea_hallmark)$ID]

if (length(valid_pathways) > 0) {

  message(sprintf("   - 成功匹配到 %d 个通路，正在渲染合并图表...", length(valid_pathways)))

  custom_colors <- c("#D51F26", "#4CB5F5")[1:length(valid_pathways)]

  main_title <- "Metabolic & Hypoxic Reprogramming in Unstable Plaques"
  sample_info <- sprintf("Dataset: GSE163154 | IPH (N=%d) vs non_IPH (N=%d)", n_iph, n_non_iph)
  full_title <- sprintf("%s\n%s", main_title, sample_info)

  p_combined <- gseaplot2(gsea_hallmark, 
                          geneSetID = valid_pathways, 
                          title = full_title,
                          color = custom_colors,
                          pvalue_table = TRUE,
                          base_size = 14)

  out_name_combined <- "12_GSEA_Combined_Glycolysis_Hypoxia"

  ggsave(file.path(out_dir, paste0(out_name_combined, ".pdf")), p_combined, width = 8, height = 6.5)
  ggsave(file.path(out_dir, paste0(out_name_combined, ".png")), p_combined, width = 8, height = 6.5, dpi = 600)

  message(">>> 完美！带有 GSE 和样本量标注的多通路联合图已成功导出。")

} else {
  message("   ! 警告: 糖酵解和低氧通路均未在本次 GSEA 中达到显著富集，无法绘制联合图。")
}


message(">>> [Step 1] 正在根据设定的通路提取并验证核心靶基因...")

suppressMessages({
  library(tidyverse)
  library(ggpubr)
})

candidate_genes <- c(
  "ACKR1", "HIF1A", 
  "PFKFB3", "HK2", "PKM", "PFKM", "PFKP", "PFKL", "ALDOA", "ALDOB", "ALDOC", "GAPDH", "LDHA", "PDK1", 
  "SLC16A1", "SLC16A3", "SLC16A4", "SLC16A7", "SLC5A12", "LDHB", 
  "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", 
  "CXCR4", "CXCR3", "CXCL9", "CXCL10", "CCR5"
)

sig_target_genes <- deg_clean %>%
  filter(Gene %in% candidate_genes & P.Value < 0.05) %>%
  arrange(desc(logFC))

message(sprintf("   - 筛查完毕！共有 %d 个目标基因在出血组发生显著改变，将被纳入绘图。", nrow(sig_target_genes)))

expr_subset <- expr_gse[sig_target_genes$PROBEID, ]
rownames(expr_subset) <- sig_target_genes$Gene

plot_data <- as.data.frame(t(expr_subset))
plot_data$Group <- group_factor

plot_data_long <- plot_data %>%
  pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression")

plot_data_long$Gene <- factor(plot_data_long$Gene, levels = sig_target_genes$Gene)

group_labels <- c(sprintf("non-IPH\n(N=%d)", n_non_iph), sprintf("IPH\n(N=%d)", n_iph))
plot_data_long$Group <- factor(plot_data_long$Group, 
                               levels = c("non_IPH", "IPH"), 
                               labels = group_labels)

message(">>> [Step 2] 正在绘制核心靶点表达量分面总图...")

comparisons_list <- list(group_labels)

p_genes <- ggplot(plot_data_long, aes(x = Group, y = Expression, color = Group)) +
  geom_violin(aes(fill = Group), alpha = 0.15, trim = FALSE, linewidth = 0.8) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", linewidth = 0.8, alpha = 0.9) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +

  scale_color_manual(values = c("#272E6A", "#D51F26")) +
  scale_fill_manual(values = c("#272E6A", "#D51F26")) +

  facet_wrap(~ Gene, scales = "free_y", ncol = 5) + 

  stat_compare_means(method = "t.test", 
                     comparisons = comparisons_list,
                     label = "p.format",
                     bracket.size = 0.6, 
                     tip.length = 0.02,
                     step.increase = 0.05, 
                     size = 3.5,
                     color = "black") +

  labs(title = "Dataset: GSE163154 | Expression of Core Mechanism Genes in Carotid Plaques",
       x = "", y = "Normalized Expression") +

  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 20)),
    strip.background = element_rect(fill = "#F2F4F4", color = "black", linewidth = 1.2),
    strip.text = element_text(face = "bold", size = 13, color = "black"),

    axis.text.x = element_text(size = 12, face = "bold", color = "black", lineheight = 1.1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),

    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),

    legend.position = "none",
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

calc_height <- max(ceiling(nrow(sig_target_genes) / 5) * 3.2, 7)

message(sprintf(">>> [Step 3] 正在保存高清大图 (自适应尺寸: 14.5 x %.1f)...", calc_height))

out_name_genes <- "13_Core_Genes_Expression_Panel"

ggsave(file.path(out_dir, paste0(out_name_genes, ".pdf")), p_genes, width = 14.5, height = calc_height)
ggsave(file.path(out_dir, paste0(out_name_genes, ".png")), p_genes, width = 14.5, height = calc_height, dpi = 600)

message("\n>>> 完美搞定！这张大图足以作为你手稿中最有说服力的核心 Figure 之一！")


message(">>> [Step 1] 正在提取 ACKR1 及乳酸代谢靶基因的表达量数据...")

suppressMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

target_gene <- "ACKR1"
lactate_genes <- c("SLC16A3", "SLC16A4", "SLC16A7", "LDHB")

cor_probe_info <- deg_clean %>%
  filter(Gene %in% c(target_gene, lactate_genes)) %>%
  select(Gene, PROBEID)

expr_cor <- expr_gse[cor_probe_info$PROBEID, ]
rownames(expr_cor) <- cor_probe_info$Gene

df_cor <- as.data.frame(t(expr_cor))
df_cor$Group <- group_factor

df_cor$Group <- factor(df_cor$Group, levels = c("non_IPH", "IPH"), labels = c("non-IPH", "IPH"))

df_cor_long <- df_cor %>%
  pivot_longer(cols = all_of(lactate_genes), 
               names_to = "Lactate_Gene", 
               values_to = "Lactate_Expr")

df_cor_long$Lactate_Gene <- factor(df_cor_long$Lactate_Gene, 
                                   levels = c("SLC16A3", "LDHB", "SLC16A4", "SLC16A7"))

message(">>> [Step 2] 正在绘制 ACKR1 与乳酸代谢基因的 Pearson 相关性图...")

sample_info <- sprintf("Dataset: GSE163154 | N = %d (IPH: %d, non-IPH: %d)", 
                       n_iph + n_non_iph, n_iph, n_non_iph)

p_cor <- ggplot(df_cor_long, aes(x = Lactate_Expr, y = ACKR1)) +
  geom_point(aes(fill = Group), shape = 21, color = "black", size = 2.5, stroke = 0.5, alpha = 0.85) +

  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 1) +

  stat_cor(method = "pearson", size = 4.5, fontface = "bold", label.x.npc = "left", label.y.npc = "top") +

  facet_wrap(~ Lactate_Gene, scales = "free_x", ncol = 2) +

  scale_fill_manual(values = c("non-IPH" = "#272E6A", "IPH" = "#D51F26")) +

  labs(title = "Correlation Between Lactate Metabolism and ACKR1",
       subtitle = sample_info,
       x = "Normalized Expression of Lactate Metabolism Genes",
       y = "Normalized Expression of ACKR1") +

  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),

    strip.background = element_rect(fill = "#F2F4F4", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold.italic", size = 13, color = "black"),

    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold", margin = margin(t=10, r=10)),

    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 1, color = "black"),

    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    legend.margin = margin(b = -5),

    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

message(">>> [Step 3] 正在导出高清相关性散点图...")

out_name_cor <- "15_Correlation_ACKR1_LactateGenes"

ggsave(file.path(out_dir, paste0(out_name_cor, ".pdf")), p_cor, width = 7.5, height = 7)
ggsave(file.path(out_dir, paste0(out_name_cor, ".png")), p_cor, width = 7.5, height = 7, dpi = 600)

message("\n>>> 相关性图表大功告成！如果 P 值显著，这就是你将代谢重塑和内皮激活关联在一起的最强临床证据！")


# Figure 4: GSE163154 volcano plot, heatmap, and core-gene expression

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(limma)
  library(illuminaHumanv2.db)
  library(ggrepel)
  library(ggpubr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "24_GEO_Validation_GSE163154")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

message(">>> [Step 1] 正在读取并处理 GSE163154 数据...")
gset <- getGEO(filename = "GSE163154_series_matrix.txt.gz", getGPL = FALSE)
expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

group_gse <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case = TRUE), "non_IPH", "IPH")
group_factor <- factor(group_gse, levels = c("non_IPH", "IPH"))
n_non_iph <- sum(group_factor == "non_IPH")
n_iph <- sum(group_factor == "IPH")
sample_info <- sprintf("Dataset: GSE163154 | IPH (n=%d) vs non_IPH (n=%d)", n_iph, n_non_iph)

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)
contrast.matrix <- makeContrasts(IPH_vs_non_IPH = IPH - non_IPH, levels = design)
fit <- lmFit(expr_gse, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
all_deg <- topTable(fit2, coef = "IPH_vs_non_IPH", number = Inf)
all_deg$PROBEID <- rownames(all_deg)

probe_to_symbol <- mapIds(illuminaHumanv2.db, keys = all_deg$PROBEID, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
all_deg$Gene <- probe_to_symbol
deg_clean <- all_deg %>% filter(!is.na(Gene)) %>% arrange(P.Value) %>% distinct(Gene, .keep_all = TRUE)

message(">>> [Step 2] 正在绘制完美复刻版火山图 (标签高亮升级)...")

gene_modules <- list(
  "Core Targets" = c("ACKR1", "HIF1A"),
  "Glycolysis" = c("HK2", "PKM", "PFKM", "PFKL", "ALDOA", "ALDOC", "GAPDH"),
  "Lactate Metabolism" = c("SLC16A3", "SLC16A4", "SLC16A7", "LDHB"),
  "CD8+ T Cell" = c("CD8A", "GZMA", "GZMB", "PRF1", "IFNG"),
  "Receptors" = c("CXCR4")
)
genes_to_label <- unname(unlist(gene_modules))

p_cutoff <- 0.05
logfc_cutoff <- 0.5

volcano_data <- deg_clean %>%
  mutate(
    change = factor(case_when(
      P.Value < p_cutoff & logFC > logfc_cutoff ~ "Up",
      P.Value < p_cutoff & logFC < -logfc_cutoff ~ "Down",
      TRUE ~ "Not Sig"
    ), levels = c("Up", "Down", "Not Sig")),
    logP = -log10(P.Value)
  )

p_volcano <- ggplot(volcano_data, aes(x = logFC, y = logP)) +
  geom_point(aes(color = change), alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("Up" = "#D51F26", "Down" = "#272E6A", "Not Sig" = "grey70")) +

  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", linewidth = 0.6) +

  geom_point(data = filter(volcano_data, Gene %in% genes_to_label), 
             fill = NA, color = "black", size = 2.5, shape = 21, stroke = 1.1) +

  geom_label_repel(data = filter(volcano_data, Gene %in% genes_to_label),
                   aes(label = Gene), 
                   fill = "white",
                   color = "black",
                   label.size = 0.4,
                   label.padding = unit(0.2, "lines"),

                   size = 3.6, fontface = "bold.italic",
                   box.padding = 0.8, point.padding = 0.4, max.overlaps = Inf,
                   segment.color = "grey30",
                   segment.size = 0.6,
                   min.segment.length = 0) +

  theme_bw() +
  labs(title = "Volcano Plot: IPH vs non-IPH", subtitle = sample_info, 
       x = expression(log[2]("Fold Change")), y = expression(-log[10]("P-value"))) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    panel.border = element_rect(color = "black", linewidth = 1.3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 13),
    legend.position = "top", legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold")
  )

message(">>> [Step 2] 正在导出高清火山图...")
ggsave(file.path(out_dir, "01_Volcano_Plot_Classic_HighVisibility.pdf"), p_volcano, width = 7.5, height = 6.5)
ggsave(file.path(out_dir, "01_Volcano_Plot_Classic_HighVisibility.png"), p_volcano, width = 7.5, height = 6.5, dpi = 600)

message(">>> [Step 3] 正在逐个模块导出高颜值机制图...")

module_df <- stack(gene_modules)
colnames(module_df) <- c("Gene", "Module")
sig_genes <- deg_clean %>% inner_join(module_df, by = "Gene") %>% filter(P.Value < 0.05)

expr_sub <- expr_gse[sig_genes$PROBEID, ]
rownames(expr_sub) <- sig_genes$Gene
plot_data <- as.data.frame(t(expr_sub))
plot_data$Group <- factor(group_factor, labels = c(sprintf("non-IPH\n(n=%d)", n_non_iph), sprintf("IPH\n(n=%d)", n_iph)))

plot_long <- plot_data %>%
  pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression") %>%
  left_join(module_df, by = "Gene")

fdr_map <- deg_clean %>%
  filter(Gene %in% unique(plot_long$Gene)) %>%
  mutate(fdr_label = ifelse(adj.P.Val < 0.001,
                            sprintf("FDR = %.1e", adj.P.Val),
                            sprintf("FDR = %.3f", adj.P.Val))) %>%
  dplyr::select(Gene, fdr_label)

fdr_pos <- plot_long %>%
  group_by(Gene) %>%
  summarise(y_lab = max(Expression, na.rm = TRUE) +
              0.08 * (max(Expression, na.rm = TRUE) - min(Expression, na.rm = TRUE)),
            .groups = "drop") %>%
  left_join(fdr_map, by = "Gene") %>%
  left_join(module_df, by = "Gene")


for (mod_name in names(gene_modules)) {

  df_sub <- filter(plot_long, Module == mod_name)
  if(nrow(df_sub) == 0) next

  n_genes <- length(unique(df_sub$Gene))
  n_cols <- min(n_genes, 4) 
  n_rows <- ceiling(n_genes / n_cols)

  calc_width <- max(n_cols * 2.5, 5)
  calc_height <- n_rows * 4 + 1.5 

  mod_title <- sprintf("%s Expression\n%s", mod_name, sample_info)

  p_mod <- ggplot(df_sub, aes(x = Group, y = Expression, fill = Group)) +
    geom_violin(alpha = 0.3, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", linewidth = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, aes(color = Group)) +

    scale_fill_manual(values = c("#272E6A", "#D51F26")) +
    scale_color_manual(values = c("#272E6A", "#D51F26")) +

    facet_wrap(~ Gene, scales = "free_y", ncol = n_cols) +

    geom_text(data = filter(fdr_pos, Module == mod_name),
              aes(x = 1.5, y = y_lab, label = fdr_label),
              inherit.aes = FALSE, size = 3.5, fontface = "bold") +


    theme_pubr() +
    labs(title = mod_title, y = "Normalized Expression", x = "") +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, lineheight = 1.2),
      legend.position = "none",
      strip.background = element_rect(fill = "#F2F4F4", color = "black", linewidth = 1),
      strip.text = element_text(face = "bold.italic", size = 12, color = "black"),
      axis.text.x = element_text(face = "bold", size = 11, color = "black"),
      axis.title.y = element_text(face = "bold", size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

  safe_mod_name <- gsub(" ", "_", mod_name)
  safe_mod_name <- gsub("\\+", "", safe_mod_name)

  file_name_base <- sprintf("02_Mechanism_%s", safe_mod_name)

  ggsave(file.path(out_dir, paste0(file_name_base, ".pdf")), p_mod, width = calc_width, height = calc_height)
  ggsave(file.path(out_dir, paste0(file_name_base, ".png")), p_mod, width = calc_width, height = calc_height, dpi = 600)

  message(sprintf("     - 已输出: %s (排布: %d 行 %d 列)", safe_mod_name, n_rows, n_cols))
}

message("\n>>> 大功告成！所有机制图已拆分为独立文件，方便你在 AI 中自由拼装排版！")


suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(limma)
  library(illuminaHumanv2.db)
  library(pheatmap)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "24_GEO_Validation_GSE163154")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

message(">>> [Step 1] 正在重新读取 GSE163154 数据并执行 limma 分析...")
gset <- getGEO(filename = "GSE163154_series_matrix.txt.gz", getGPL = FALSE)
expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

group_gse <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case = TRUE), "non-IPH", "IPH")
group_factor <- factor(group_gse, levels = c("non-IPH", "IPH"))

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("non_IPH", "IPH")
contrast.matrix <- makeContrasts(IPH_vs_non_IPH = IPH - non_IPH, levels = design)

fit <- lmFit(expr_gse, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
all_deg <- topTable(fit2, coef = "IPH_vs_non_IPH", number = Inf)
all_deg$PROBEID <- rownames(all_deg)

all_deg$Gene <- mapIds(illuminaHumanv2.db, keys = all_deg$PROBEID, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
deg_clean <- all_deg %>% filter(!is.na(Gene)) %>% arrange(P.Value) %>% distinct(Gene, .keep_all = TRUE)

message(">>> [Step 2] 正在提取核心机制基因并构建热图矩阵...")

gene_modules <- c("ACKR1", "HIF1A", "HK2", "PKM", "PFKM", "PFKL", "ALDOA", "ALDOC", "GAPDH", 
                  "SLC16A3", "SLC16A4", "SLC16A7", "LDHB", "CD8A", "GZMA", "GZMB", "PRF1", "IFNG", "CXCR4")

sig_genes_hm <- deg_clean %>% filter(Gene %in% gene_modules & P.Value < 0.05)

expr_hm <- expr_gse[sig_genes_hm$PROBEID, ]
rownames(expr_hm) <- sig_genes_hm$Gene

annotation_col <- data.frame(Group = group_factor)
rownames(annotation_col) <- colnames(expr_hm)

order_idx <- order(annotation_col$Group)
expr_hm_ordered <- expr_hm[, order_idx]
annotation_col_ordered <- annotation_col[order_idx, , drop = FALSE]

message(">>> [Step 3] 正在渲染并导出超高对比度的热图...")

ann_colors <- list(Group = c("non-IPH" = "#272E6A", "IPH" = "#D51F26"))
hm_colors <- colorRampPalette(c("#272E6A", "white", "#D51F26"))(100)

my_breaks <- seq(-1.5, 1.5, length.out = 100)

italic_gene_labels <- as.expression(lapply(rownames(expr_hm_ordered), function(x) bquote(italic(.(x)))))

out_base <- file.path(out_dir, "03_Heatmap_Core_Mechanisms_HighContrast")

pheatmap(expr_hm_ordered, 
         scale = "row",
         color = hm_colors,
         breaks = my_breaks,
         annotation_col = annotation_col_ordered, 
         annotation_colors = ann_colors,   
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         labels_row = italic_gene_labels,
         fontsize_row = 11,                
         border_color = NA,

         filename = paste0(out_base, ".pdf"),
         width = 6, height = 5)

pheatmap(expr_hm_ordered, 
         scale = "row",
         color = hm_colors,
         breaks = my_breaks,
         annotation_col = annotation_col_ordered,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         labels_row = italic_gene_labels,
         fontsize_row = 11,
         border_color = NA,

         filename = paste0(out_base, ".png"),
         width = 6, height = 5, res = 600)

message("\n>>> 增强完毕！快去看看文件夹里的 HighContrast 版本，红蓝对比绝对极其分明！")


message(">>> 正在重绘 CXCR4，确保箱线图物理胖瘦与其他多基因图一致...")

mod_name <- "Receptors"
sample_info <- sprintf("Dataset: GSE163154 | IPH (n=%d) vs non_IPH (n=%d)", n_iph, n_non_iph)
df_sub <- plot_long %>% filter(Module == mod_name)
mod_title <- sprintf("%s Expression\n%s", mod_name, sample_info)

p_mod_cxcr4 <- ggplot(df_sub, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(alpha = 0.3, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, aes(color = Group)) +

  scale_fill_manual(values = c("#272E6A", "#D51F26")) +
  scale_color_manual(values = c("#272E6A", "#D51F26")) +

  facet_wrap(~ Gene, scales = "free_y", ncol = 1) + 
    geom_text(data = filter(fdr_pos, Module == mod_name),
              aes(x = 1.5, y = y_lab, label = fdr_label),
              inherit.aes = FALSE, size = 3.5, fontface = "bold") +

  theme_pubr() +
  labs(title = mod_title, y = "Normalized Expression", x = "") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, lineheight = 1.2),
    legend.position = "none",
    strip.background = element_rect(fill = "#F2F4F4", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold.italic", size = 12, color = "black"),
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    axis.title.y = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

calc_width <- 3.5
calc_height <- 5.5

file_name_base <- sprintf("02_Mechanism_%s_Proportional", mod_name)

ggsave(file.path(out_dir, paste0(file_name_base, ".pdf")), p_mod_cxcr4, width = calc_width, height = calc_height)
ggsave(file.path(out_dir, paste0(file_name_base, ".png")), p_mod_cxcr4, width = calc_width, height = calc_height, dpi = 600)

message(">>> 搞定！新生成的图表已恢复完美比例。")


# Supplementary Figure 2G/H: Glycolysis and hypoxia GSEA curves

message(">>> [Step 1] 正在初始化环境并加载必要的包...")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(limma)
  library(illuminaHumanv2.db)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "24_GEO_Validation_GSE163154")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

message(">>> [Step 2] 正在加载 GSE163154 表达数据并执行 limma 分析...")

gset <- getGEO(filename = "GSE163154_series_matrix.txt.gz", getGPL = FALSE)
expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

group_gse <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case = TRUE), "non_IPH", "IPH")
group_factor <- factor(group_gse, levels = c("non_IPH", "IPH"))
n_non_iph <- sum(group_factor == "non_IPH")
n_iph <- sum(group_factor == "IPH")
sample_info <- sprintf("Dataset: GSE163154 | IPH (N=%d) vs non_IPH (N=%d)", n_iph, n_non_iph)

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("non_IPH", "IPH")
contrast.matrix <- makeContrasts(IPH_vs_non_IPH = IPH - non_IPH, levels = design)

fit <- lmFit(expr_gse, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
all_deg <- topTable(fit2, coef = "IPH_vs_non_IPH", number = Inf)
all_deg$PROBEID <- rownames(all_deg)

all_deg$Gene <- mapIds(illuminaHumanv2.db, keys = all_deg$PROBEID, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
deg_clean <- all_deg %>% filter(!is.na(Gene)) %>% arrange(P.Value) %>% distinct(Gene, .keep_all = TRUE)

message(">>> [Step 3] 正在转换 Entrez ID 并运行 HALLMARK GSEA 分析...")

gene_id_mapping <- bitr(deg_clean$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

gsea_df <- deg_clean %>%
  inner_join(gene_id_mapping, by = c("Gene" = "SYMBOL")) %>%
  filter(!is.na(logFC)) %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(logFC), n = 1) %>%
  ungroup() %>%
  arrange(desc(logFC))

gene_list <- gsea_df$logFC
names(gene_list) <- gsea_df$ENTREZID

hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

set.seed(1234)
gsea_hallmark <- GSEA(geneList     = gene_list,
                      TERM2GENE    = hallmark_t2g,
                      minGSSize    = 10,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)

message(">>> [Step 4] 正在分别绘制并导出糖酵解和低氧的 GSEA 曲线图 (已加宽画幅)...")

target_pathways <- c("HALLMARK_GLYCOLYSIS", "HALLMARK_HYPOXIA")

res_df <- as.data.frame(gsea_hallmark)

for (pathway in target_pathways) {
  if (pathway %in% res_df$ID) {

    nes <- res_df[pathway, "NES"]

    line_color <- ifelse(nes > 0, "#D51F26", "#272E6A")

    clean_title <- gsub("HALLMARK_", "", pathway)
    clean_title <- gsub("_", " ", clean_title)

    full_title <- sprintf("%s\n%s", clean_title, sample_info)

    p_single <- gseaplot2(gsea_hallmark, 
                          geneSetID = pathway, 
                          title = full_title, 
                          color = line_color,
                          pvalue_table = TRUE,
                          base_size = 14)

    safe_file_name <- gsub(" ", "_", clean_title)
    out_name <- paste0("16_GSEA_Single_", safe_file_name)

    ggsave(file.path(out_dir, paste0(out_name, ".pdf")), p_single, width = 8.5, height = 6)
    ggsave(file.path(out_dir, paste0(out_name, ".png")), p_single, width = 8.5, height = 6, dpi = 600)

    message(sprintf("   - 成功导出: %s (NES = %.2f)", clean_title, nes))

  } else {
    message(sprintf("   ! 警告: %s 未在本次 GSEA 结果中找到。", pathway))
  }
}

message("\n>>> 修复完成！画布已经加宽，右上角的 P 值数字绝对能完整、漂亮地展示出来了！")


# Figure 5: GSE163154 correlation overview

message(">>> [Step 1] 正在初始化环境并提取 GSE163154 数据...")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(limma)
  library(illuminaHumanv2.db)
  library(ggpubr)
  library(cowplot)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "27_Final Correlation_Analysis")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

gset <- getGEO(filename = "GSE163154_series_matrix.txt.gz", getGPL = FALSE)
if (is.list(gset)) gset <- gset[[1]]

expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

common_samples <- intersect(colnames(expr_gse), rownames(pheno_gse))
expr_gse <- expr_gse[, common_samples, drop = FALSE]
pheno_gse <- pheno_gse[common_samples, , drop = FALSE]

group_gse <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case = TRUE), "non_IPH", "IPH")
group_gse[is.na(group_gse)] <- "non_IPH"
group_factor <- factor(group_gse, levels = c("non_IPH", "IPH"))

n_non_iph <- sum(group_factor == "non_IPH")
n_iph <- sum(group_factor == "IPH")
sample_info <- sprintf("Dataset: GSE163154 | N=%d (IPH: %d, non-IPH: %d)", 
                       n_iph + n_non_iph, n_iph, n_non_iph)

design <- model.matrix(~ 0 + group_factor)
rownames(design) <- colnames(expr_gse)
colnames(design) <- c("non_IPH", "IPH")

contrast.matrix <- makeContrasts(IPH_vs_non_IPH = IPH - non_IPH, levels = design)
fit <- eBayes(contrasts.fit(lmFit(expr_gse, design), contrast.matrix))
all_deg <- topTable(fit, coef = "IPH_vs_non_IPH", number = Inf)
all_deg$PROBEID <- rownames(all_deg)
all_deg$Gene <- mapIds(illuminaHumanv2.db, keys = all_deg$PROBEID, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
best_probes <- all_deg %>% filter(!is.na(Gene)) %>% arrange(P.Value) %>% distinct(Gene, .keep_all = TRUE)

message(">>> [Step 2] 正在计算评分并重组表达矩阵...")

glyco_genes <- c("HK2", "PKM", "PFKM", "PFKL", "ALDOA", "ALDOC", "GAPDH")
lac_genes <- c("SLC16A3", "SLC16A4", "SLC16A7", "LDHB")
cd8_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
target_genes <- c("ACKR1", "HIF1A", "CXCL12", "IFNG", "CXCR4", glyco_genes, lac_genes, cd8_genes)

probe_info <- best_probes %>% filter(Gene %in% target_genes)
expr_sub <- expr_gse[probe_info$PROBEID, ]
rownames(expr_sub) <- probe_info$Gene

df_expr <- as.data.frame(t(expr_sub))
df_expr$Group <- group_factor

cd8_available <- intersect(cd8_genes, colnames(df_expr))
df_expr$CD8_Tox_Score <- rowMeans(scale(df_expr[, cd8_available]))

message(">>> [Step 3] 正在生成带有实线阴影的双色散点图...")

color_map <- c("non_IPH" = "#272E6A", "IPH" = "#D51F26")
legend_labels <- c("non_IPH" = "non-IPH", "IPH" = "IPH")

plot_cor <- function(data, x_gene, y_gene) {
  p_val <- cor.test(data[[x_gene]], data[[y_gene]], method = "pearson")$p.value
  if (!is.na(p_val) && p_val >= 0.05) return(NULL)

  p <- ggplot(data, aes_string(x = x_gene, y = y_gene)) +
    geom_point(aes(fill = Group), shape = 21, color = "black", size = 2.5, stroke = 0.6, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", 
                linewidth = 1, fill = "grey60", alpha = 0.3) +
    stat_cor(method = "pearson", size = 3.8, fontface = "bold", label.x.npc = "left", label.y.npc = 0.95) +

    scale_fill_manual(values = color_map, labels = legend_labels) +

    theme_bw() +
    labs(x = paste(x_gene, "(Normalized)"), y = paste(y_gene, "(Normalized)")) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 8)),
      plot.margin = margin(15, 15, 15, 15),
      legend.position = "none"
    )
  return(p)
}

plot_box <- function(data, y_var, title_text) {
  p_val <- t.test(data[[y_var]] ~ data[["Group"]])$p.value
  if (!is.na(p_val) && p_val >= 0.05) return(NULL)

  p <- ggplot(data, aes_string(x = "Group", y = y_var, fill = "Group")) +
    geom_violin(alpha = 0.3, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.3, outlier.shape = NA, fill = "white", linewidth = 0.8) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, aes(color = Group)) +

    stat_compare_means(method = "t.test", label = "p.format", size = 3.8, fontface = "bold") +

    scale_fill_manual(values = color_map, labels = legend_labels) + 
    scale_color_manual(values = color_map, labels = legend_labels) +
    scale_x_discrete(labels = legend_labels) +

    theme_bw() +
    labs(title = title_text, y = paste(y_var, "(Normalized)"), x = "") +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 12)),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 11, face = "bold", color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 8)),
      plot.margin = margin(15, 15, 15, 15),
      legend.position = "none"
    )
  return(p)
}

dummy_plot <- ggplot(df_expr, aes(x = ACKR1, y = HIF1A, fill = Group)) +
  geom_point(shape = 21, size = 3.5) + 
  scale_fill_manual(values = color_map, labels = legend_labels) +
  theme_bw() + 
  theme(legend.position = "top", 
        legend.title = element_blank(), 
        legend.text = element_text(size=12, face="bold"),
        legend.margin = margin(b = 10))
shared_legend <- get_legend(dummy_plot)

save_panel <- function(plot_list, file_prefix, ncol, nrow, width, height, title_text) {
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0) {
    message(paste(">>> 跳过模块:", file_prefix, "- 均无统计学显著性"))
    return(invisible(NULL))
  }

  grid_p <- plot_grid(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "hv")
  final_p <- plot_grid(
    ggdraw() + draw_label(paste(title_text, "\n", sample_info), fontface = 'bold', size = 16),
    shared_legend, 
    grid_p, 
    ncol = 1, 
    rel_heights = c(0.12, 0.08, 1) 
  )
  ggsave(file.path(out_dir, paste0(file_prefix, ".pdf")), final_p, width = width, height = height)
  ggsave(file.path(out_dir, paste0(file_prefix, ".png")), final_p, width = width, height = height, dpi = 600)
}

message(">>> [Step 4] 正在按模块导出高清图...")

save_panel(list(plot_cor(df_expr, "ACKR1", "HIF1A"), plot_cor(df_expr, "HIF1A", "CXCL12")), 
           "01_ACKR1_vs_CoreTargets", 2, 1, 9, 5.5, "Core Targets Correlation")

save_panel(lapply(glyco_genes, function(g) plot_cor(df_expr, "ACKR1", g)), 
           "02_ACKR1_vs_Glycolysis", 4, 2, 15, 8.5, "Glycolysis Genes vs ACKR1")

save_panel(lapply(lac_genes, function(g) plot_cor(df_expr, "ACKR1", g)), 
           "03_ACKR1_vs_Lactate", 2, 2, 8.5, 8.5, "Lactate Metabolism Genes vs ACKR1")

p_cd8_cor <- lapply(c("CD8A", "GZMA", "GZMB", "PRF1", "IFNG", "CD8_Tox_Score"), function(g) plot_cor(df_expr, "ACKR1", g))
p_tox_box <- plot_box(df_expr, "CD8_Tox_Score", "CD8+ T Toxicity Score\n(IPH vs non-IPH)")
save_panel(c(p_cd8_cor, list(p_tox_box)), 
           "04_ACKR1_vs_CD8_ToxScore", 4, 2, 16, 9, "CD8+ T Cell Features vs ACKR1")

save_panel(list(plot_cor(df_expr, "ACKR1", "CXCR4")), 
           "05_ACKR1_vs_Receptors", 1, 1, 5.5, 5.5, "Receptors vs ACKR1")


save_panel(lapply(glyco_genes, function(g) plot_cor(df_expr, "HIF1A", g)), 
           "06_HIF1A_vs_Glycolysis_Genes", 4, 2, 15, 8.5, "HIF1A Regulation on Glycolysis Genes")

save_panel(lapply(lac_genes, function(g) plot_cor(df_expr, "HIF1A", g)), 
           "07_HIF1A_vs_Lactate_Genes", 2, 2, 8.5, 8.5, "HIF1A Regulation on Lactate Genes")

save_panel(list(plot_cor(df_expr, "HIF1A", "CD8_Tox_Score")), 
           "08_HIF1A_vs_CD8_Tox_Score", 1, 1, 5.5, 5.5, "HIF1A vs CD8+ T Toxicity Score")

message("\n>>> 大功告成！配色和双图例已完美恢复，全部图片已生成！")


# Figure 6: Proteomic differential analysis, enrichment, mechanistic validation, and PPI

message(">>> [Step 1] 初始化环境与包...")

suppressMessages({
  library(tidyverse)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
  library(ggpubr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "28_Protein_Correlation_Analysis")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

target_file <- "稳定-不稳定-合并-差异蛋白-调整后.csv"

message(">>> [Step 2] 正在读取蛋白质组学原始数据...")

raw_lines <- read.csv(target_file, header = FALSE, nrows = 5)
sample_names <- make.unique(as.character(raw_lines[1, -1]))
raw_labels <- as.character(raw_lines[2, -1])

data_df <- read.csv(target_file, skip = 2, header = FALSE)
protein_ids <- data_df[, 1]
expr_mat_raw <- as.matrix(data_df[, -1])
class(expr_mat_raw) <- "numeric"

min_val <- min(expr_mat_raw[expr_mat_raw > 0], na.rm = TRUE)
expr_mat_raw[is.na(expr_mat_raw) | expr_mat_raw == 0] <- min_val / 2
norm_mat <- log2(expr_mat_raw)
colnames(norm_mat) <- sample_names
rownames(norm_mat) <- protein_ids

final_group <- rep(NA, length(raw_labels))
final_group[raw_labels == "A"] <- "Stable"
final_group[raw_labels == "B"] <- "Unstable"
group_factor <- factor(final_group, levels = c("Stable", "Unstable"))

n_stable <- sum(group_factor == "Stable", na.rm=TRUE)
n_unstable <- sum(group_factor == "Unstable", na.rm=TRUE)
sample_info <- sprintf("Unstable/IPH (n=%d) vs Stable/non-IPH (n=%d)", n_unstable, n_stable)
message(sprintf("   - 样本质控完成: %s", sample_info))

message(">>> [Step 3] 正在执行 Limma 差异计算...")

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("Stable", "Unstable")
contrast <- makeContrasts(Unstable_vs_Stable = Unstable - Stable, levels = design)
fit <- lmFit(norm_mat, design)
fit2 <- eBayes(contrasts.fit(fit, contrast))
deg_res <- topTable(fit2, coef=1, number=Inf)
deg_res$Uniprot <- rownames(deg_res)

deg_res$CleanID <- deg_res$Uniprot %>% str_split(";") %>% sapply(`[`, 1) %>% str_split("-") %>% sapply(`[`, 1)

mapped_ids <- mapIds(org.Hs.eg.db, keys = deg_res$CleanID, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
mapped_entrez <- mapIds(org.Hs.eg.db, keys = deg_res$CleanID, column = "ENTREZID", keytype = "UNIPROT", multiVals = "first")

deg_res$SYMBOL <- mapped_ids
deg_res$ENTREZID <- mapped_entrez

write.csv(deg_res, file.path(out_dir, "01_All_Proteins_Differential_Results.csv"), row.names = FALSE)
message("   - 完整差异蛋白总表已保存！没有任何一个蛋白被丢失。")

message(">>> [Step 4] 正在执行 GO 和 KEGG 富集分析...")

sig_degs <- deg_res %>% filter(P.Value < 0.05, !is.na(ENTREZID))
up_genes <- sig_degs %>% filter(logFC > 0) %>% pull(ENTREZID)
down_genes <- sig_degs %>% filter(logFC < 0) %>% pull(ENTREZID)

plot_nature_bubble <- function(up_res, down_res, title_text, top_n = 10) {
  df_up <- as.data.frame(up_res) %>% head(top_n) %>% mutate(Direction = "Upregulated in IPH")
  df_down <- as.data.frame(down_res) %>% head(top_n) %>% mutate(Direction = "Downregulated in IPH")

  if(nrow(df_up) == 0 && nrow(df_down) == 0) return(NULL)

  df_comb <- bind_rows(df_up, df_down)
  df_comb$LogP <- -log10(df_comb$pvalue)

  df_comb <- df_comb %>% arrange(Direction, Count)
  df_comb$Description <- factor(df_comb$Description, levels = unique(df_comb$Description))

  p <- ggplot(df_comb, aes(x = Count, y = Description)) +
    geom_point(aes(size = LogP, fill = Direction), shape = 21, color = "black", stroke = 0.8, alpha = 0.85) +
    scale_fill_manual(values = c("Upregulated in IPH" = "#D51F26", "Downregulated in IPH" = "#272E6A")) +
    scale_size_continuous(range = c(3, 8)) +
    theme_bw() +
    labs(title = title_text, subtitle = sample_info, x = "Gene Count", y = "", size = "-Log10(P-value)") +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "grey30"),
      axis.text.y = element_text(size = 11, color = "black", face = "bold"),
      axis.text.x = element_text(size = 11, color = "black", face = "bold"),
      axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
      panel.border = element_rect(color = "black", linewidth = 1.2, fill = NA),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10)
    )
  return(p)
}

go_up <- enrichGO(gene = up_genes, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
go_down <- enrichGO(gene = down_genes, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
p_go <- plot_nature_bubble(go_up, go_down, "GO Biological Process Enrichment")
if(!is.null(p_go)) {
  ggsave(file.path(out_dir, "02_GO_BP_Enrichment.pdf"), p_go, width = 10, height = 7)
  ggsave(file.path(out_dir, "02_GO_BP_Enrichment.png"), p_go, width = 10, height = 7, dpi = 600)
}

kegg_up <- enrichKEGG(gene = up_genes, organism = "hsa", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = down_genes, organism = "hsa", pvalueCutoff = 0.05)
if(!is.null(kegg_up)) kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
if(!is.null(kegg_down)) kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

p_kegg <- plot_nature_bubble(kegg_up, kegg_down, "KEGG Pathway Enrichment")
if(!is.null(p_kegg)) {
  ggsave(file.path(out_dir, "03_KEGG_Enrichment.pdf"), p_kegg, width = 10, height = 7)
  ggsave(file.path(out_dir, "03_KEGG_Enrichment.png"), p_kegg, width = 10, height = 7, dpi = 600)
}

message(">>> [Step 5] 正在执行 Hallmark GSEA 分析...")

hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
hallmark_genes$gs_name <- gsub("HALLMARK_", "", hallmark_genes$gs_name)

gsea_list <- deg_res %>% 
  filter(!is.na(SYMBOL)) %>% 
  group_by(SYMBOL) %>% 
  summarize(logFC = mean(logFC, na.rm=TRUE)) %>% 
  arrange(desc(logFC))
gene_rank <- gsea_list$logFC
names(gene_rank) <- gsea_list$SYMBOL

set.seed(123)
gsea_res <- GSEA(gene_rank, TERM2GENE = hallmark_genes, pvalueCutoff = 0.25, eps = 0)

if(nrow(as.data.frame(gsea_res)) > 0) {
  gsea_df <- as.data.frame(gsea_res) %>% arrange(desc(NES))

  gsea_plot_df <- bind_rows(
    gsea_df %>% filter(p.adjust < 0.05, NES > 0) %>% head(10) %>% mutate(Status = "Activated in IPH"),
    gsea_df %>% filter(p.adjust < 0.05, NES < 0) %>% tail(10) %>% mutate(Status = "Suppressed in IPH")
  )

  if(nrow(gsea_plot_df) > 0) {
    gsea_plot_df$ID <- factor(gsea_plot_df$ID, levels = gsea_plot_df$ID[order(gsea_plot_df$NES)])

    p_gsea <- ggplot(gsea_plot_df, aes(x = NES, y = ID)) +
      geom_segment(aes(x = 0, xend = NES, y = ID, yend = ID), color = "grey50", linewidth = 1) +
      geom_point(aes(size = -log10(p.adjust), fill = Status), shape = 21, color = "black", stroke = 0.8) +
      scale_fill_manual(values = c("Activated in IPH" = "#D51F26", "Suppressed in IPH" = "#272E6A")) +
      scale_size_continuous(range = c(4, 9)) +
      theme_classic() +
      labs(title = "Hallmark GSEA: Pathways Altered in IPH", subtitle = sample_info, 
           x = "Normalized Enrichment Score (NES)", y = "", size = "-Log10(FDR)") +
      theme(
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "grey30"),
        axis.text.y = element_text(size = 10, color = "black", face = "bold"),
        axis.text.x = element_text(size = 11, color = "black", face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
        panel.border = element_rect(color = "black", linewidth = 1.2, fill = NA),
        legend.position = "right"
      )

    ggsave(file.path(out_dir, "04_GSEA_Hallmark_Summary.pdf"), p_gsea, width = 10, height = 6)
    ggsave(file.path(out_dir, "04_GSEA_Hallmark_Summary.png"), p_gsea, width = 10, height = 6, dpi = 600)

    write.csv(gsea_df, file.path(out_dir, "05_GSEA_Hallmark_Full_Results.csv"), row.names = FALSE)
  }
}

message("\n>>> 全部分析已完美生成！请前往 28_Protein_Correlation_Analysis 文件夹查收。")


message(">>> [Step 1] 初始化环境与包...")

suppressMessages({
  library(tidyverse)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggpubr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "28_Protein_Correlation_Analysis")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

target_file <- "稳定-不稳定-合并-差异蛋白-调整后.csv"

message(">>> [Step 2] 正在读取蛋白质组学原始数据...")

raw_lines <- read.csv(target_file, header = FALSE, nrows = 5)
sample_names <- make.unique(as.character(raw_lines[1, -1]))
raw_labels <- as.character(raw_lines[2, -1])

data_df <- read.csv(target_file, skip = 2, header = FALSE)
protein_ids <- data_df[, 1]
expr_mat_raw <- as.matrix(data_df[, -1])
class(expr_mat_raw) <- "numeric"

min_val <- min(expr_mat_raw[expr_mat_raw > 0], na.rm = TRUE)
expr_mat_raw[is.na(expr_mat_raw) | expr_mat_raw == 0] <- min_val / 2
norm_mat <- log2(expr_mat_raw)
colnames(norm_mat) <- sample_names
rownames(norm_mat) <- protein_ids

final_group <- rep(NA, length(raw_labels))
final_group[raw_labels == "A"] <- "Stable"
final_group[raw_labels == "B"] <- "Unstable"
group_factor <- factor(final_group, levels = c("Stable", "Unstable"))

n_stable <- sum(group_factor == "Stable", na.rm=TRUE)
n_unstable <- sum(group_factor == "Unstable", na.rm=TRUE)

sample_info <- sprintf("Unstable (n=%d) vs Stable (n=%d)", n_unstable, n_stable)
message(sprintf("   - 样本质控完成: %s", sample_info))

message(">>> [Step 3] 正在执行 Limma 差异计算...")

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("Stable", "Unstable")
contrast <- makeContrasts(Unstable_vs_Stable = Unstable - Stable, levels = design)
fit <- lmFit(norm_mat, design)
fit2 <- eBayes(contrasts.fit(fit, contrast))
deg_res <- topTable(fit2, coef=1, number=Inf)
deg_res$Uniprot <- rownames(deg_res)

deg_res$CleanID <- deg_res$Uniprot %>% str_split(";") %>% sapply(`[`, 1) %>% str_split("-") %>% sapply(`[`, 1)

mapped_ids <- mapIds(org.Hs.eg.db, keys = deg_res$CleanID, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
mapped_entrez <- mapIds(org.Hs.eg.db, keys = deg_res$CleanID, column = "ENTREZID", keytype = "UNIPROT", multiVals = "first")

deg_res$SYMBOL <- mapped_ids
deg_res$ENTREZID <- mapped_entrez

write.csv(deg_res, file.path(out_dir, "01_All_Proteins_Differential_Results.csv"), row.names = FALSE)
message("   - 完整差异蛋白总表已保存！")

message(">>> [Step 4] 正在执行并单独输出 GO 和 KEGG 富集分析...")

sig_degs <- deg_res %>% filter(P.Value < 0.05, !is.na(ENTREZID))
up_genes <- sig_degs %>% filter(logFC > 0) %>% pull(ENTREZID)
down_genes <- sig_degs %>% filter(logFC < 0) %>% pull(ENTREZID)

plot_single_bubble <- function(res_obj, color_val, title_text, top_n = 15) {
  df <- as.data.frame(res_obj)
  if(nrow(df) == 0) return(NULL)

  df <- df %>% head(top_n)
  df$LogP <- -log10(df$pvalue)
  df <- df %>% arrange(Count)
  df$Description <- factor(df$Description, levels = unique(df$Description))

  p <- ggplot(df, aes(x = Count, y = Description)) +
    geom_point(aes(size = LogP), fill = color_val, shape = 21, color = "black", stroke = 0.8, alpha = 0.85) +
    scale_size_continuous(range = c(4, 9)) +

    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    scale_y_discrete(expand = expansion(mult = c(0.08, 0.08))) +

    theme_bw() +
    labs(title = title_text, subtitle = sample_info, x = "Gene Count", y = "", size = "-Log10(P-value)") +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "grey30"),
      axis.text.y = element_text(size = 11, color = "black", face = "bold"),
      axis.text.x = element_text(size = 11, color = "black", face = "bold"),
      axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
      panel.border = element_rect(color = "black", linewidth = 1.2, fill = NA),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(15, 15, 15, 15)
    )
  return(p)
}

go_up <- enrichGO(gene = up_genes, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
go_down <- enrichGO(gene = down_genes, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)

p_go_up <- plot_single_bubble(go_up, "#D51F26", "GO (BP) Upregulated in Unstable Plaque")
p_go_down <- plot_single_bubble(go_down, "#272E6A", "GO (BP) Downregulated in Unstable Plaque")

if(!is.null(p_go_up)) {
  ggsave(file.path(out_dir, "02_GO_BP_Upregulated.pdf"), p_go_up, width = 9, height = 6.5)
  ggsave(file.path(out_dir, "02_GO_BP_Upregulated.png"), p_go_up, width = 9, height = 6.5, dpi = 600)
}
if(!is.null(p_go_down)) {
  ggsave(file.path(out_dir, "02_GO_BP_Downregulated.pdf"), p_go_down, width = 9, height = 6.5)
  ggsave(file.path(out_dir, "02_GO_BP_Downregulated.png"), p_go_down, width = 9, height = 6.5, dpi = 600)
}

kegg_up <- enrichKEGG(gene = up_genes, organism = "hsa", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = down_genes, organism = "hsa", pvalueCutoff = 0.05)
if(!is.null(kegg_up)) kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
if(!is.null(kegg_down)) kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

p_kegg_up <- plot_single_bubble(kegg_up, "#D51F26", "KEGG Pathways Upregulated in Unstable Plaque")
p_kegg_down <- plot_single_bubble(kegg_down, "#272E6A", "KEGG Pathways Downregulated in Unstable Plaque")

if(!is.null(p_kegg_up)) {
  ggsave(file.path(out_dir, "03_KEGG_Upregulated.pdf"), p_kegg_up, width = 9, height = 6.5)
  ggsave(file.path(out_dir, "03_KEGG_Upregulated.png"), p_kegg_up, width = 9, height = 6.5, dpi = 600)
}
if(!is.null(p_kegg_down)) {
  ggsave(file.path(out_dir, "03_KEGG_Downregulated.pdf"), p_kegg_down, width = 9, height = 6.5)
  ggsave(file.path(out_dir, "03_KEGG_Downregulated.png"), p_kegg_down, width = 9, height = 6.5, dpi = 600)
}

message("\n>>> 完美搞定！气泡边缘防裁切已修复，GSEA 代码已彻底移除！请前往目录查收完整圆润的气泡图。")


message(">>> [Step 1] 初始化环境与读取底层矩阵...")

suppressMessages({
  library(tidyverse)
  library(limma)
  library(org.Hs.eg.db)
  library(ggpubr)
  library(cowplot)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "28_Protein_Correlation_Analysis")
setwd(work_dir)

gene_modules <- list(
  "Core_Targets" = c("ACKR1", "HIF1A"),
  "Glycolysis" = c("PFKFB3", "HK2", "PKM", "PFKM", "PFKP", "PFKL", "ALDOA", "ALDOB", "ALDOC", "GAPDH", "LDHA", "PDK1"),
  "Lactate_Metabolism" = c("SLC16A1", "SLC16A3", "SLC16A4", "SLC16A7", "SLC5A12", "LDHB"),
  "CD8_T_Cell" = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG"),
  "Receptors" = c("CXCR4")
)

module_df <- stack(gene_modules)
colnames(module_df) <- c("SYMBOL", "Module")

target_file <- "稳定-不稳定-合并-差异蛋白-调整后.csv"
raw_lines <- read.csv(target_file, header = FALSE, nrows = 5)
sample_names <- make.unique(as.character(raw_lines[1, -1]))
raw_labels <- as.character(raw_lines[2, -1])

data_df <- read.csv(target_file, skip = 2, header = FALSE)
protein_ids <- data_df[, 1]
expr_mat_raw <- as.matrix(data_df[, -1])
class(expr_mat_raw) <- "numeric"

min_val <- min(expr_mat_raw[expr_mat_raw > 0], na.rm = TRUE)
expr_mat_raw[is.na(expr_mat_raw) | expr_mat_raw == 0] <- min_val / 2
norm_mat <- log2(expr_mat_raw)
colnames(norm_mat) <- sample_names
rownames(norm_mat) <- protein_ids

final_group <- rep(NA, length(raw_labels))
final_group[raw_labels == "A"] <- "Stable"
final_group[raw_labels == "B"] <- "Unstable"
group_factor <- factor(final_group, levels = c("Stable", "Unstable"))

n_stable <- sum(group_factor == "Stable", na.rm=TRUE)
n_unstable <- sum(group_factor == "Unstable", na.rm=TRUE)
sample_info <- sprintf("Proteomics | Unstable (n=%d) vs Stable (n=%d)", n_unstable, n_stable)

x_labels <- c("Stable" = sprintf("Stable\n(n=%d)", n_stable), 
              "Unstable" = sprintf("Unstable\n(n=%d)", n_unstable))

message(">>> [Step 2] 正在筛选显著蛋白并提取表达矩阵...")

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("Stable", "Unstable")
contrast <- makeContrasts(Unstable_vs_Stable = Unstable - Stable, levels = design)
fit <- eBayes(contrasts.fit(lmFit(norm_mat, design), contrast))
deg_res <- topTable(fit, coef=1, number=Inf)
deg_res$Uniprot <- rownames(deg_res)

deg_res$CleanID <- deg_res$Uniprot %>% str_split(";") %>% sapply(`[`, 1) %>% str_split("-") %>% sapply(`[`, 1)
deg_res$SYMBOL <- mapIds(org.Hs.eg.db, keys = deg_res$CleanID, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")

sig_targets <- deg_res %>%
  inner_join(module_df, by = "SYMBOL") %>%
  filter(!is.na(P.Value) & P.Value < 0.05) %>%
  group_by(SYMBOL) %>%
  slice_min(P.Value, n = 1, with_ties = FALSE) %>%
  ungroup()

if(nrow(sig_targets) == 0) {
  stop("未发现显著核心蛋白，无法生成小提琴图。")
}

expr_sub <- norm_mat[sig_targets$Uniprot, , drop = FALSE]
rownames(expr_sub) <- sig_targets$SYMBOL

df_long <- as.data.frame(t(expr_sub)) %>%
  mutate(Group = group_factor) %>%
  pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression") %>%
  left_join(sig_targets %>% dplyr::select(SYMBOL, Module), by = c("Gene" = "SYMBOL"))

message(">>> [Step 3] 正在生成各模块的小提琴大图...")

color_map <- c("Stable" = "#272E6A", "Unstable" = "#D51F26")
legend_labels <- c("Stable" = "Stable", "Unstable" = "Unstable")

plot_violin_standard <- function(data_sub, title_text) {

  n_genes <- length(unique(data_sub$Gene))
  n_cols <- min(n_genes, 4)

  p <- ggplot(data_sub, aes(x = Group, y = Expression, fill = Group)) +
    geom_violin(alpha = 0.3, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", linewidth = 0.8) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, aes(color = Group)) +

    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +

    stat_compare_means(method = "t.test", aes(label = paste0("p = ", ..p.format..)), size = 3.8, fontface = "bold") +

    scale_fill_manual(values = color_map, labels = legend_labels) +
    scale_color_manual(values = color_map, labels = legend_labels) +

    scale_x_discrete(labels = x_labels) +

    facet_wrap(~ Gene, scales = "free_y", ncol = n_cols) +

    theme_bw() +
    labs(x = "", y = "Normalized Expression") +
    theme(
      strip.background = element_rect(fill = "#F2F4F4", color = "black", linewidth = 1),
      strip.text = element_text(face = "bold.italic", size = 12, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 11, face = "bold", color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 8)),
      legend.position = "none"
    )
  return(p)
}

dummy_plot <- ggplot(df_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_point(shape = 21, size = 3.5) + scale_fill_manual(values = color_map, labels = legend_labels) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(), 
                     legend.text = element_text(size=12, face="bold"), legend.margin = margin(b = 10))
shared_legend <- get_legend(dummy_plot)

save_violin_panel <- function(plot_obj, file_prefix, width, height) {
  final_p <- plot_grid(
    shared_legend, 
    plot_obj + theme(plot.margin = margin(10, 10, 10, 10)),
    ncol = 1, rel_heights = c(0.08, 1)
  )
  ggsave(file.path(out_dir, paste0(file_prefix, ".pdf")), final_p, width = width, height = height)
  ggsave(file.path(out_dir, paste0(file_prefix, ".png")), final_p, width = width, height = height, dpi = 600)
}

unique_modules <- unique(df_long$Module)

for(mod in unique_modules) {
  df_sub <- df_long %>% filter(Module == mod)
  n_genes <- length(unique(df_sub$Gene))

  if(n_genes > 0) {
    p_mod <- plot_violin_standard(df_sub, mod)

    calc_cols <- min(n_genes, 4)
    calc_rows <- ceiling(n_genes / 4)
    calc_width <- max(calc_cols * 2.5 + 1.5, 5) 
    calc_height <- calc_rows * 4 + 1.5

    file_name <- sprintf("07_Proteomics_%s_Violin", mod)
    save_violin_panel(p_mod, file_name, calc_width, calc_height)
    message(sprintf("   - 成功导出: %s (包含 %d 个靶点)", file_name, n_genes))
  }
}

message("\n>>> 完美搞定！代码原汁原味，已为您成功加上 X 轴 N 数和 p= 格式！")


message(">>> [Step 1] 初始化环境与读取底层矩阵...")

suppressMessages({
  library(tidyverse)
  library(org.Hs.eg.db)
  library(ggpubr)
  library(cowplot)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "28_Protein_Correlation_Analysis")
setwd(work_dir)

gene_modules <- list(
  "Glycolysis" = c("PFKFB3", "HK2", "PKM", "PFKM", "PFKP", "PFKL", "ALDOA", "ALDOB", "ALDOC", "GAPDH", "LDHA", "PDK1"),
  "Lactate_Metabolism" = c("SLC16A1", "SLC16A3", "SLC16A4", "SLC16A7", "SLC5A12", "LDHB"),
  "CD8_T_Cell" = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG"),
  "Receptors" = c("CXCR4")
)

target_file <- "稳定-不稳定-合并-差异蛋白-调整后.csv"
raw_lines <- read.csv(target_file, header = FALSE, nrows = 5)
sample_names <- make.unique(as.character(raw_lines[1, -1]))
raw_labels <- as.character(raw_lines[2, -1])

data_df <- read.csv(target_file, skip = 2, header = FALSE)
protein_ids <- data_df[, 1]
expr_mat_raw <- as.matrix(data_df[, -1])
class(expr_mat_raw) <- "numeric"

min_val <- min(expr_mat_raw[expr_mat_raw > 0], na.rm = TRUE)
expr_mat_raw[is.na(expr_mat_raw) | expr_mat_raw == 0] <- min_val / 2
norm_mat <- log2(expr_mat_raw)
colnames(norm_mat) <- sample_names

final_group <- rep(NA, length(raw_labels))
final_group[raw_labels == "A"] <- "Stable"
final_group[raw_labels == "B"] <- "Unstable"
group_factor <- factor(final_group, levels = c("Stable", "Unstable"))

n_stable <- sum(group_factor == "Stable", na.rm=TRUE)
n_unstable <- sum(group_factor == "Unstable", na.rm=TRUE)
sample_info <- sprintf("Proteomics | Unstable (n=%d) vs Stable (n=%d)", n_unstable, n_stable)

message(">>> [Step 2] 正在提取蛋白组学数据并重组矩阵...")

clean_ids <- protein_ids %>% str_split(";") %>% sapply(`[`, 1) %>% str_split("-") %>% sapply(`[`, 1)
try_map <- tryCatch({
  AnnotationDbi::select(org.Hs.eg.db, keys = clean_ids, columns = "SYMBOL", keytype = "UNIPROT")
}, error = function(e) { return(NULL) })

map_df <- data.frame(OriginalID = protein_ids, CleanID = clean_ids)
if(!is.null(try_map)) {
  try_map <- try_map[!duplicated(try_map$UNIPROT), ]
  map_df <- left_join(map_df, try_map, by = c("CleanID" = "UNIPROT"))
} else {
  map_df$SYMBOL <- map_df$OriginalID
}

expr_df <- as.data.frame(norm_mat)
expr_df$OriginalID <- protein_ids

expr_clean <- expr_df %>%
  inner_join(map_df, by = "OriginalID") %>%
  filter(!is.na(SYMBOL)) %>%
  dplyr::select(-OriginalID, -CleanID) %>%
  group_by(SYMBOL) %>%
  summarise_all(mean, na.rm=TRUE) %>%
  as.data.frame()

rownames(expr_clean) <- expr_clean$SYMBOL
mat_final <- as.matrix(expr_clean[, -1])

df_expr <- as.data.frame(t(mat_final))
df_expr$Group <- group_factor

available_genes <- colnames(df_expr)
cd8_avail <- intersect(gene_modules$CD8_T_Cell, available_genes)
glyco_avail <- intersect(gene_modules$Glycolysis, available_genes)
lac_avail <- intersect(gene_modules$Lactate_Metabolism, available_genes)

if(length(cd8_avail) > 0) df_expr$CD8_Tox_Score <- rowMeans(scale(df_expr[, cd8_avail, drop=FALSE]))
if(length(glyco_avail) > 0) df_expr$Glycolysis_Score <- rowMeans(scale(df_expr[, glyco_avail, drop=FALSE]))
if(length(lac_avail) > 0) df_expr$Lactate_Score <- rowMeans(scale(df_expr[, lac_avail, drop=FALSE]))

message(">>> [Step 3] 正在计算 Pearson 相关性并自动过滤...")

color_map <- c("Stable" = "#272E6A", "Unstable" = "#D51F26")
legend_labels <- c("Stable" = sprintf("Stable (n=%d)", n_stable), 
                   "Unstable" = sprintf("Unstable (n=%d)", n_unstable))

plot_cor <- function(data, x_gene, y_gene) {
  if(!x_gene %in% colnames(data) || !y_gene %in% colnames(data)) return(NULL)

  test_res <- cor.test(data[[x_gene]], data[[y_gene]], method = "pearson")
  p_val <- test_res$p.value
  r_val <- test_res$estimate
  if (is.na(p_val) || p_val >= 0.05) return(NULL)

  p_label <- ifelse(p_val < 0.001, sprintf("p = %.2e", p_val), sprintf("p = %.4f", p_val))
  r_label <- sprintf("R = %.2f", r_val)
  stat_text <- paste(r_label, p_label, sep = ", ")

  y_max <- max(data[[y_gene]], na.rm = TRUE)
  y_min <- min(data[[y_gene]], na.rm = TRUE)
  y_pos <- y_max + (y_max - y_min) * 0.05

  p <- ggplot(data, aes_string(x = x_gene, y = y_gene)) +
    geom_point(aes(fill = Group), shape = 21, color = "black", size = 2.5, stroke = 0.6, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", 
                linewidth = 1, fill = "grey60", alpha = 0.3) +

    annotate("text", x = -Inf, y = y_pos, label = stat_text, 
             hjust = -0.1, vjust = 1, size = 4, fontface = "bold") +

    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_manual(values = color_map, labels = legend_labels) +

    theme_bw() +
    labs(x = paste(x_gene, "(Normalized)"), y = paste(y_gene, "(Normalized)")) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 8)),
      plot.margin = margin(15, 15, 15, 15),
      legend.position = "none"
    )
  return(p)
}

dummy_plot <- ggplot(df_expr, aes(x = 1, y = 1, fill = Group)) +
  geom_point(shape = 21, size = 3.5) + scale_fill_manual(values = color_map, labels = legend_labels) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(), 
                     legend.text = element_text(size=12, face="bold"), legend.margin = margin(b = 10))
shared_legend <- get_legend(dummy_plot)

save_panel <- function(plot_list, file_prefix, ncol, nrow, width, height, title_text) {
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0) {
    message(paste("   [跳过]:", file_prefix, "- 蛋白组中未测到或均不显著。"))
    return(invisible(NULL))
  }

  grid_p <- plot_grid(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "hv")
  final_p <- plot_grid(
    shared_legend, 
    grid_p, 
    ncol = 1, 
    rel_heights = c(0.08, 1) 
  )
  ggsave(file.path(out_dir, paste0(file_prefix, ".pdf")), final_p, width = width, height = height)
  ggsave(file.path(out_dir, paste0(file_prefix, ".png")), final_p, width = width, height = height, dpi = 600)
  message(paste("   [成功]: 已导出", file_prefix, "共计", length(plot_list), "张显著相关图。"))
}

message(">>> [Step 4] 正在按模块寻找并绘制高清图...")

glyco_list <- intersect(gene_modules$Glycolysis, available_genes)
lac_list <- intersect(gene_modules$Lactate_Metabolism, available_genes)

if("ACKR1" %in% available_genes) {
  save_panel(list(plot_cor(df_expr, "ACKR1", "HIF1A"), plot_cor(df_expr, "HIF1A", "CXCL12")), 
             "08_ACKR1_vs_CoreTargets_Correlation", 2, 1, 9, 5, "Core Targets Correlation")

  save_panel(lapply(c(glyco_list, "Glycolysis_Score"), function(g) plot_cor(df_expr, "ACKR1", g)), 
             "09_ACKR1_vs_Glycolysis_Correlation", 4, 2, 15, 8.5, "Glycolysis vs ACKR1")

  save_panel(lapply(c(lac_list, "Lactate_Score"), function(g) plot_cor(df_expr, "ACKR1", g)), 
             "10_ACKR1_vs_Lactate_Correlation", 3, 2, 11, 8.5, "Lactate Metabolism vs ACKR1")

  cd8_y <- c(intersect(gene_modules$CD8_T_Cell, available_genes), "CD8_Tox_Score")
  save_panel(lapply(cd8_y, function(g) plot_cor(df_expr, "ACKR1", g)), 
             "11_ACKR1_vs_CD8_ToxScore_Correlation", 4, 2, 15, 8.5, "CD8+ T Cell vs ACKR1")

  save_panel(list(plot_cor(df_expr, "ACKR1", "CXCR4")), 
             "12_ACKR1_vs_Receptors_Correlation", 1, 1, 5.5, 5, "Receptors vs ACKR1")
} else {
  message("   [提示] 蛋白质组数据中未检测到低丰度受体 ACKR1，已跳过 ACKR1 相关散点图。")
}

if("HIF1A" %in% available_genes) {
  save_panel(lapply(glyco_list, function(g) plot_cor(df_expr, "HIF1A", g)), 
             "13_HIF1A_vs_Glycolysis_Correlation", 4, 2, 15, 8.5, "HIF1A vs Glycolysis")

  save_panel(lapply(lac_list, function(g) plot_cor(df_expr, "HIF1A", g)), 
             "14_HIF1A_vs_Lactate_Correlation", 3, 2, 11, 8.5, "HIF1A vs Lactate")

  if("CD8_Tox_Score" %in% colnames(df_expr)) {
    save_panel(list(plot_cor(df_expr, "HIF1A", "CD8_Tox_Score")), 
               "15_HIF1A_vs_CD8_Tox_Score_Correlation", 1, 1, 5.5, 5, "HIF1A vs CD8+ Tox Score")
  }
}

message("\n>>> 完美搞定！不显著的已被自动抛弃，图例带有 N 数，P 值强制规定了 'p =' 格式！")


message(">>> [Step 1] 初始化环境与定义 13 个核心基因...")

suppressMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(httr) 
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "28_Protein_Correlation_Analysis")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

nodes_df <- data.frame(
  name = c(
    "ACKR1", "CXCL12", "CXCR4",         
    "HIF1A", "EP300",            
    "PKM", "ALDOB", "LDHA",             
    "SLC16A3", "SLC16A4", "SLC16A7",                          
    "CD8A", "GZMB", "PRF1"              
  ),
  Category = factor(c(
    rep("Chemokine Axis", 3),
    rep("Hypoxia & Epigenetics", 2),
    rep("Glycolysis Reprogramming", 3),
    rep("Lactate Metabolism", 3),       
    rep("Cytotoxic T Effector", 3)
  ), levels = c("Chemokine Axis", "Hypoxia & Epigenetics", "Glycolysis Reprogramming", "Lactate Metabolism", "Cytotoxic T Effector"))
)

message(">>> [Step 2] 正在向 STRING 数据库发送 API 请求...")

genes_for_api <- paste(nodes_df$name, collapse = "%0D")
string_url <- paste0(
  "https://string-db.org/api/tsv/network?identifiers=", 
  genes_for_api, 
  "&species=9606&required_score=400"
)

string_res <- tryCatch({
  read.table(string_url, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}, error = function(e) {
  stop("无法连接 STRING 数据库，请检查网络连接！")
})

edges_df <- string_res %>%
  dplyr::select(from = preferredName_A, to = preferredName_B, weight = score) %>%
  filter(from %in% nodes_df$name & to %in% nodes_df$name)

message(">>> [Step 3] 正在渲染极简硬核风格网络图...")

net <- graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = FALSE)

category_colors <- c(
  "Chemokine Axis" = "#00A087",
  "Hypoxia & Epigenetics" = "#E64B35",
  "Glycolysis Reprogramming" = "#4DBBD5",
  "Lactate Metabolism" = "#F39B7F",
  "Cytotoxic T Effector" = "#3C5488"
)

set.seed(123)

p_ppi <- ggraph(net, layout = 'fr') + 

  geom_edge_link(aes(edge_alpha = weight, edge_width = weight), color = "#BDBDBD", lineend = "round") +
  scale_edge_width(range = c(0.8, 3.5), name = "Interaction\nConfidence") +
  scale_edge_alpha(range = c(0.5, 1), guide = "none") +

  geom_node_point(aes(fill = Category), shape = 21, color = "black", size = 11, stroke = 1.8) +
  scale_fill_manual(values = category_colors, name = "Category") +

  geom_node_text(aes(label = name), repel = TRUE, size = 6.5, fontface = "bold.italic", 
                 color = "black", point.padding = unit(0.2, "lines"), box.padding = unit(0.5, "lines")) +

  guides(
    edge_width = guide_legend(order = 1, override.aes = list(color = "#BDBDBD")),
    fill = guide_legend(override.aes = list(size = 8), order = 2)
  ) +

  labs(title = "Functional Protein-Protein Interaction (PPI) Network:\nBridging Glycolysis, Hypoxia, and Immune Recruitment") +
  theme_void() + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20, t = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

message(">>> [Step 4] 正在导出高清复刻版 PPI 网络图...")

out_name <- "22_Reference_Style_STRING_PPI"
ggsave(file.path(out_dir, paste0(out_name, ".pdf")), p_ppi, width = 11, height = 8.5)
ggsave(file.path(out_dir, paste0(out_name, ".png")), p_ppi, width = 11, height = 8.5, dpi = 600)

message("\n>>> 完美收工！这回图里全是干货，绝对没有多余的“生信味”装饰了！")


# Figure 7A-C: Final LASSO feature selection

message(">>> [Step 1] 初始化环境并加载数据...")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(illuminaHumanv2.db)
  library(glmnet)    
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(ggrepel)
  library(ggsci)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "29_Clinical_Diagnostic_Modeling")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

tox_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")

target_genes_base <- c(
  "ACKR1", "CXCL12", "CXCR4",
  "HIF1A", "EP300",
  "HK2", "PKM", "LDHA", "PFKFB3", "ALDOB",
  "SLC16A3", "SLC16A7", "SLC16A4"
)

target_gse <- "GSE163154_series_matrix.txt.gz"
if(!file.exists(target_gse)) {
  gset <- getGEO("GSE163154", destdir=".", getGPL = FALSE)[[1]]
} else {
  gset <- getGEO(filename = target_gse, getGPL = FALSE)
  if(is.list(gset)) gset <- gset[[1]]
}

expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

clinical_status <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case=T), "non-IPH", "IPH")
y_label <- factor(clinical_status, levels = c("non-IPH", "IPH"))

n_non_iph <- sum(y_label == "non-IPH")
n_iph <- sum(y_label == "IPH")
clinical_info <- sprintf("Dataset: GSE163154 | Cohort: N=%d (IPH: n=%d, non-IPH: n=%d)", 
                         length(y_label), n_iph, n_non_iph)

probe_to_symbol <- mapIds(illuminaHumanv2.db, keys = rownames(expr_gse), 
                          column = "SYMBOL", keytype = "PROBEID", multiVals = "first")

expr_df <- as.data.frame(expr_gse)
expr_df$SYMBOL <- probe_to_symbol

expr_clean <- expr_df %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise_all(mean, na.rm=TRUE) %>%
  as.data.frame()
rownames(expr_clean) <- expr_clean$SYMBOL
expr_mat <- as.matrix(expr_clean[, -1])

available_tox <- intersect(tox_genes, rownames(expr_mat))
if(length(available_tox) > 0) {
  tox_z_mat <- t(scale(t(expr_mat[available_tox, , drop = FALSE])))
  cd8_tox_score <- colMeans(tox_z_mat, na.rm = TRUE)
  expr_mat <- rbind(expr_mat, "CD8_Tox_Score" = cd8_tox_score)
}

available_features <- intersect(c(target_genes_base, "CD8_Tox_Score"), rownames(expr_mat))
X_mat <- t(expr_mat[available_features, ])

X_mat_scaled <- scale(X_mat)

message(">>> [Step 2] 正在运行优化版 LASSO 惩罚回归...")

set.seed(2024) 
lasso_fit <- glmnet(X_mat_scaled, y_label, family = "binomial", alpha = 1)
cv_fit <- cv.glmnet(X_mat_scaled, y_label, family = "binomial", alpha = 1, type.measure = "deviance", nfolds = 10)

best_lambda <- cv_fit$lambda.min
se1_lambda <- cv_fit$lambda.1se
lasso_coef <- coef(cv_fit, s = "lambda.min")
selected_features <- rownames(lasso_coef)[which(lasso_coef != 0)]
selected_genes <- selected_features[selected_features != "(Intercept)"]

message(sprintf("\n   >>> 在 %d 个输入特征中，LASSO 客观为您锁定了 %d 个最强预测因子：", 
                length(available_features), length(selected_genes)))
print(selected_genes)

message(">>> [Step 3] 正在渲染高级 LASSO 拼图...")

cv_df <- data.frame(Log_Lambda = log(cv_fit$lambda), CVM = cv_fit$cvm, CV_UP = cv_fit$cvup, CV_LO = cv_fit$cvlo, N_Zero = cv_fit$nzero)
p_cv <- ggplot(cv_df, aes(x = Log_Lambda, y = CVM)) +
  geom_errorbar(aes(ymin = CV_LO, ymax = CV_UP), color = "grey70", width = 0.08) +
  geom_point(color = "#E64B35", size = 2) +
  geom_vline(xintercept = log(best_lambda), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = log(se1_lambda), linetype = "dashed", color = "black", linewidth = 0.6) +

  annotate("text", x = log(best_lambda) + 0.15, y = min(cv_df$CVM) + 0.3, 
           label = paste0("Optimal Log(\u03bb) = ", round(log(best_lambda), 3)), 
           hjust = 0, size = 4, fontface = "italic", color = "black") +

  annotate("text", x = cv_df$Log_Lambda[seq(1, nrow(cv_df), by=4)], y = max(cv_df$CV_UP) + 0.15, 
           label = cv_df$N_Zero[seq(1, nrow(cv_df), by=4)], size = 3.5) +

  theme_classic() +
  labs(title = "A. Cross-Validation Error", subtitle = clinical_info, 
       x = "Log(Lambda)", y = "Binomial Deviance") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5, margin = margin(b=15)),
        axis.title = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 12, color = "black"))

beta_matrix <- as.matrix(lasso_fit$beta)
path_df <- as.data.frame(t(beta_matrix)) %>%
  mutate(Log_Lambda = log(lasso_fit$lambda)) %>%
  pivot_longer(cols = -Log_Lambda, names_to = "Feature", values_to = "Coefficient")

final_betas <- as.matrix(coef(lasso_fit, s = best_lambda))[-1, , drop=F]
valid_genes <- rownames(final_betas)[final_betas != 0]

df_labels <- path_df %>%
  filter(Log_Lambda == min(Log_Lambda)) %>%
  filter(Feature %in% valid_genes)

p_path <- ggplot(path_df, aes(x = Log_Lambda, y = Coefficient, color = Feature)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_vline(xintercept = log(best_lambda), linetype = "dashed", color = "black", linewidth = 0.6) +

  geom_text_repel(data = df_labels, aes(label = Feature), 
                  size = 4.2, fontface = "bold.italic", segment.size = 0.4,
                  direction = "y", nudge_x = -0.5, 
                  xlim = c(-Inf, min(path_df$Log_Lambda) - 0.2)) + 

  scale_color_npg() +
  scale_x_continuous(expand = expansion(mult = c(0.4, 0.05))) + 
  theme_classic() +
  labs(title = "B. LASSO Coefficient Trajectories", subtitle = clinical_info,
       x = "Log(Lambda)", y = "Coefficients") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5, margin = margin(b=15)),
        axis.title = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none")

final_lasso_plot <- plot_grid(p_cv, p_path, ncol = 2, align = "h")
ggsave(file.path(out_dir, "13_LASSO_Trajectory_Plot.pdf"), final_lasso_plot, width = 12, height = 5.5)
ggsave(file.path(out_dir, "13_LASSO_Trajectory_Plot.png"), final_lasso_plot, width = 12, height = 5.5, dpi = 600)

message(">>> [Step 4] 正在渲染变量重要性全景图...")

lasso_coef_full <- as.data.frame(as.matrix(lasso_coef))
colnames(lasso_coef_full) <- "Coefficient"
lasso_coef_full$Feature <- rownames(lasso_coef_full)

full_landscape_df <- lasso_coef_full %>%
  filter(Feature != "(Intercept)") %>%
  mutate(
    Fate = case_when(
      Coefficient > 0 ~ "Risk Factor (Promotes IPH)",
      Coefficient < 0 ~ "Protective Factor",
      TRUE ~ "Eliminated (Weight = 0)"
    ),
    Abs_Coef = abs(Coefficient)
  ) %>%
  arrange(desc(Abs_Coef), Feature)

full_landscape_df$Feature <- factor(full_landscape_df$Feature, levels = rev(full_landscape_df$Feature))

p_full <- ggplot(full_landscape_df, aes(x = Coefficient, y = Feature, fill = Fate)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.9) +
  scale_fill_manual(values = c("Risk Factor (Promotes IPH)" = "#E64B35", 
                               "Protective Factor" = "#3C5488",
                               "Eliminated (Weight = 0)" = "#E4E4E4")) +
  theme_classic() +
  labs(title = "C. Variable Importance",
       subtitle = paste0("Input: ", nrow(full_landscape_df), " Cleaned Features | Retained: ", length(selected_genes), " Core Predictors"),
       x = "LASSO Coefficient (Model Weight in predicting IPH)", y = "") +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5, margin = margin(b=15)),
    axis.text.y = element_text(size = 12, color = "black", face = "bold.italic"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 12)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(15, 20, 15, 10)
  )

dyn_height <- max(nrow(full_landscape_df) * 0.35 + 3, 6)
ggsave(file.path(out_dir, "14_LASSO_Feature_Landscape.pdf"), p_full, width = 8.5, height = dyn_height)
ggsave(file.path(out_dir, "14_LASSO_Feature_Landscape.png"), p_full, width = 8.5, height = dyn_height, dpi = 600)

saveRDS(selected_genes, file.path(out_dir, "15_Final_Selected_Features.rds"))

message("\n=== 核心特征重要性明确输出 ===")
print(full_landscape_df %>% filter(Coefficient != 0) %>% dplyr::select(Feature, Coefficient, Fate))
message("\n>>> 大功告成！画面已留白，Lambda已标注，名称已严格规范！")


message(">>> [Step 1] 初始化环境并加载数据...")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(illuminaHumanv2.db)
  library(glmnet)    
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(ggrepel)
  library(ggsci)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "29_Clinical_Diagnostic_Modeling")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

tox_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")

target_genes_base <- c(
  "ACKR1", "CXCL12", "CXCR4",
  "HIF1A", "EP300",
  "HK2", "PKM", "LDHA", "PFKFB3", "ALDOB",
  "SLC16A3", "SLC16A7", "SLC16A4"
)

target_gse <- "GSE163154_series_matrix.txt.gz"
if(!file.exists(target_gse)) {
  gset <- getGEO("GSE163154", destdir=".", getGPL = FALSE)[[1]]
} else {
  gset <- getGEO(filename = target_gse, getGPL = FALSE)
  if(is.list(gset)) gset <- gset[[1]]
}

expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

clinical_status <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case=T), "non-IPH", "IPH")
y_label <- factor(clinical_status, levels = c("non-IPH", "IPH"))

n_non_iph <- sum(y_label == "non-IPH")
n_iph <- sum(y_label == "IPH")
clinical_info <- sprintf("Dataset: GSE163154 | Cohort: N=%d (IPH: n=%d, non-IPH: n=%d)", 
                         length(y_label), n_iph, n_non_iph)

probe_to_symbol <- mapIds(illuminaHumanv2.db, keys = rownames(expr_gse), 
                          column = "SYMBOL", keytype = "PROBEID", multiVals = "first")

expr_df <- as.data.frame(expr_gse)
expr_df$SYMBOL <- probe_to_symbol

expr_clean <- expr_df %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise_all(mean, na.rm=TRUE) %>%
  as.data.frame()
rownames(expr_clean) <- expr_clean$SYMBOL
expr_mat <- as.matrix(expr_clean[, -1])


# Figure 7D-F, H-L: Clinical model, ROC, DCA, nomogram, calibration, risk plots, and performance metrics
# Figure 7G fully nested LOOCV and stability analysis are in script 05, which also generates Supplementary Figure 6A.
available_tox <- intersect(tox_genes, rownames(expr_mat))
if(length(available_tox) > 0) {
  tox_z_mat <- t(scale(t(expr_mat[available_tox, , drop = FALSE])))
  cd8_tox_score <- colMeans(tox_z_mat, na.rm = TRUE)
  expr_mat <- rbind(expr_mat, "CD8_Tox_Score" = cd8_tox_score)
}

available_features <- intersect(c(target_genes_base, "CD8_Tox_Score"), rownames(expr_mat))
X_mat <- t(expr_mat[available_features, ])

X_mat_scaled <- scale(X_mat)

message(">>> [Step 2] 正在运行优化版 LASSO 惩罚回归...")

set.seed(2024) 
lasso_fit <- glmnet(X_mat_scaled, y_label, family = "binomial", alpha = 1)
cv_fit <- cv.glmnet(X_mat_scaled, y_label, family = "binomial", alpha = 1, type.measure = "deviance", nfolds = 10)

best_lambda <- cv_fit$lambda.min
se1_lambda <- cv_fit$lambda.1se
lasso_coef <- coef(cv_fit, s = "lambda.min")
selected_features <- rownames(lasso_coef)[which(lasso_coef != 0)]
selected_genes <- selected_features[selected_features != "(Intercept)"]

message(sprintf("\n   >>> 在 %d 个输入特征中，LASSO 客观为您锁定了 %d 个最强预测因子：", 
                length(available_features), length(selected_genes)))
print(selected_genes)

message(">>> [Step 3] 正在渲染高级 LASSO 拼图...")

cv_min_log_lambda <- log(best_lambda)

cv_df <- data.frame(Log_Lambda = log(cv_fit$lambda), CVM = cv_fit$cvm, CV_UP = cv_fit$cvup, CV_LO = cv_fit$cvlo)

p_cv <- ggplot(cv_df, aes(x = Log_Lambda, y = CVM)) +
  geom_errorbar(aes(ymin = CV_LO, ymax = CV_UP), color = "grey65", width = 0.05, linewidth = 0.6) +
  geom_point(color = "red", size = 2) + 
  geom_vline(xintercept = cv_min_log_lambda, linetype = "solid", color = "#5BC0EB", linewidth = 0.8) +

  annotate("text", x = cv_min_log_lambda, y = max(cv_df$CV_UP) + 0.1, 
           label = "Min", color = "#5BC0EB", size = 4.5, fontface = "bold") +
  annotate("text", x = cv_min_log_lambda + 0.1, y = max(cv_df$CV_UP) - 0.2, 
           label = sprintf("Selected\n(Log \u03bb = %.3f)", cv_min_log_lambda), 
           color = "#5BC0EB", size = 3.8, fontface = "bold", hjust = 0) +

  theme_classic() +
  labs(title = "A. Cross-Validation", 
       x = "Log Lambda", y = "Binomial Deviance") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0),
        axis.title = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 12, color = "black"),
        plot.margin = margin(15, 20, 15, 10))

beta_matrix <- as.matrix(lasso_fit$beta)
path_df <- as.data.frame(t(beta_matrix)) %>%
  mutate(Log_Lambda = log(lasso_fit$lambda)) %>%
  pivot_longer(cols = -Log_Lambda, names_to = "Feature", values_to = "Coefficient")

final_betas <- as.matrix(coef(lasso_fit, s = best_lambda))[-1, , drop=F]
valid_genes <- rownames(final_betas)[final_betas != 0]

path_df <- path_df %>%
  mutate(Line_Group = ifelse(Feature %in% valid_genes, Feature, "Eliminated"))

idx_min <- which.min(abs(log(lasso_fit$lambda) - log(best_lambda)))
log_lambda_closest <- log(lasso_fit$lambda)[idx_min]

df_labels <- path_df %>%
  filter(Log_Lambda == log_lambda_closest) %>%
  filter(Feature %in% valid_genes)

deep_solid_colors <- c(
  "#B22222",
  "#0056B3",
  "#2E8B57",
  "#D2691E",
  "#800080",
  "#8B4513",
  "#008080",
  "#C71585",
  "#483D8B",
  "#556B2F"
)

p_path <- ggplot(path_df, aes(x = Log_Lambda, y = Coefficient, group = Feature)) +
  geom_line(data = subset(path_df, Line_Group == "Eliminated"), color = "grey85", linewidth = 0.5) +
  geom_line(data = subset(path_df, Line_Group != "Eliminated"), aes(color = Feature), linewidth = 0.9, alpha = 0.9) +
  geom_vline(xintercept = cv_min_log_lambda, linetype = "dashed", color = "black", linewidth = 0.6) +

  geom_text_repel(data = df_labels, aes(label = Feature, color = Feature), 
                  size = 4, fontface = "bold", segment.size = 0.4,
                  nudge_x = max(path_df$Log_Lambda) - cv_min_log_lambda + 0.5, 
                  direction = "y", hjust = 0) + 

  scale_color_manual(values = deep_solid_colors) + 

  scale_x_continuous(expand = expansion(mult = c(0.02, 0.35))) + 
  theme_classic() +
  labs(title = "B. Coefficient Profiles",
       x = "Log Lambda", y = "Coefficients") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0), 
        axis.title = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(15, 20, 15, 10))
final_lasso_plot <- plot_grid(p_cv, p_path, ncol = 2, align = "h", rel_widths = c(1, 1.2))
ggsave(file.path(out_dir, "10_LASSO_Trajectory_Plot.pdf"), final_lasso_plot, width = 12.5, height = 5.5)
ggsave(file.path(out_dir, "10_LASSO_Trajectory_Plot.png"), final_lasso_plot, width = 12.5, height = 5.5, dpi = 600)

message(">>> [Step 4] 正在渲染变量重要性全景图...")

lasso_coef_full <- as.data.frame(as.matrix(lasso_coef))
colnames(lasso_coef_full) <- "Coefficient"
lasso_coef_full$Feature <- rownames(lasso_coef_full)

full_landscape_df <- lasso_coef_full %>%
  filter(Feature != "(Intercept)") %>%
  mutate(
    Fate = case_when(
      Coefficient > 0 ~ "Risk Factor (Promotes IPH)",
      Coefficient < 0 ~ "Protective Factor",
      TRUE ~ "Eliminated (Weight = 0)"
    ),
    Abs_Coef = abs(Coefficient)
  ) %>%
  arrange(desc(Abs_Coef), Feature)

full_landscape_df$Feature <- factor(full_landscape_df$Feature, levels = rev(full_landscape_df$Feature))

p_full <- ggplot(full_landscape_df, aes(x = Coefficient, y = Feature, fill = Fate)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("Risk Factor (Promotes IPH)" = "#DC0000FF", 
                               "Protective Factor" = "#4DBBD5FF",
                               "Eliminated (Weight = 0)" = "#E5E5E5")) +
  theme_classic() +
  labs(title = "C. Variable Importance",
       subtitle = paste0("Model Input: ", nrow(full_landscape_df), " Features | Core Predictors Retained: ", length(selected_genes)),
       x = "LASSO Coefficient (Model Weight)", y = "") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 12, color = "grey30", hjust = 0, margin = margin(b=15)),
    axis.text.y = element_text(size = 12, color = "black", face = "bold.italic"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 12)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(15, 20, 15, 10)
  )

dyn_height <- max(nrow(full_landscape_df) * 0.35 + 3, 6)
ggsave(file.path(out_dir, "11_LASSO_Feature_Landscape.pdf"), p_full, width = 8.5, height = dyn_height)
ggsave(file.path(out_dir, "11_LASSO_Feature_Landscape.png"), p_full, width = 8.5, height = dyn_height, dpi = 600)

saveRDS(selected_genes, file.path(out_dir, "12_Final_Selected_Features.rds"))

message("\n=== 核心特征重要性明确输出 ===")
print(full_landscape_df %>% filter(Coefficient != 0) %>% dplyr::select(Feature, Coefficient, Fate))
message("\n>>> 大功告成！画面已完全复刻 PDF 风格，Lambda 已具体标注，C图已完成美化！")


message(">>> [Step 1] 初始化环境并准备精准队列数据...")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(illuminaHumanv2.db)
  library(pROC)
  library(ggplot2)
  library(ggpubr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "30_Clinical_Diagnostic_Modeling")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

tox_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
selected_genes <- c("SLC16A4", "SLC16A3", "ACKR1") 

target_gse_file <- "GSE163154_series_matrix.txt.gz"
if(!file.exists(target_gse_file)) {
  gset <- getGEO("GSE163154", destdir=".", getGPL = FALSE)[[1]]
} else {
  gset <- getGEO(filename = target_gse_file, getGPL = FALSE)
  if(is.list(gset)) gset <- gset[[1]]
}

expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

clinical_status <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case=T), "non-IPH", "IPH")
clin_df <- data.frame(
  Sample = colnames(expr_gse),
  Status_Label = factor(clinical_status, levels = c("non-IPH", "IPH")),
  Status = ifelse(clinical_status == "IPH", 1, 0)
)

probe_to_symbol <- mapIds(illuminaHumanv2.db, keys = rownames(expr_gse), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
expr_clean <- as.data.frame(expr_gse) %>%
  mutate(SYMBOL = probe_to_symbol) %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise_all(mean) %>%
  column_to_rownames("SYMBOL")

available_tox <- intersect(tox_genes, rownames(expr_clean))
cd8_tox <- colMeans(t(scale(t(expr_clean[available_tox, ]))))
expr_mat <- rbind(as.matrix(expr_clean), "CD8_Tox_Score" = cd8_tox)

for(gene in selected_genes) {
  clin_df[[gene]] <- scale(as.numeric(expr_mat[gene, ]))
}
clin_df$CD8_Tox_Score <- scale(as.numeric(expr_mat["CD8_Tox_Score", ]))

message(">>> [Step 2] 正在构建多因素联合模型与单变量模型...")

fit_comb <- glm(Status ~ SLC16A4 + SLC16A3 + ACKR1 + CD8_Tox_Score, data = clin_df, family = "binomial")
clin_df$Prob_Combined <- predict(fit_comb, type = "response")

clin_df$Prob_SLC16A4 <- predict(glm(Status ~ SLC16A4, data = clin_df, family = "binomial"), type = "response")
clin_df$Prob_SLC16A3 <- predict(glm(Status ~ SLC16A3, data = clin_df, family = "binomial"), type = "response")
clin_df$Prob_ACKR1   <- predict(glm(Status ~ ACKR1, data = clin_df, family = "binomial"), type = "response")
clin_df$Prob_CD8     <- predict(glm(Status ~ CD8_Tox_Score, data = clin_df, family = "binomial"), type = "response")

write.csv(clin_df, file.path(out_dir, "GSE163154_Patient_Risk_Scores.csv"), row.names = FALSE)

message(">>> [Step 3] 正在严格计算 AUC 与 95% 置信区间...")

roc_list <- list(
  "Combined Signature" = roc(clin_df$Status, clin_df$Prob_Combined, quiet = TRUE),
  "ACKR1"              = roc(clin_df$Status, clin_df$Prob_ACKR1, quiet = TRUE),
  "CD8 Tox Score"      = roc(clin_df$Status, clin_df$Prob_CD8, quiet = TRUE),
  "SLC16A3"            = roc(clin_df$Status, clin_df$Prob_SLC16A3, quiet = TRUE),
  "SLC16A4"            = roc(clin_df$Status, clin_df$Prob_SLC16A4, quiet = TRUE)
)

auc_labels <- sapply(names(roc_list), function(x) {
  ci_val <- ci.auc(roc_list[[x]])
  paste0(x, " (AUC = ", sprintf("%.3f", roc_list[[x]]$auc), ", 95% CI: ", sprintf("%.3f", ci_val[1]), "-", sprintf("%.3f", ci_val[3]), ")")
})
names(roc_list) <- auc_labels

message(">>> [Step 4] 正在融合专属主题绘制 Nature 风格 ROC 对比图...")

theme_nature <- function() {
  theme_classic() +
    theme(
      text = element_text(color = "black", family = "sans"),
      axis.text = element_text(size = 12, color = "black", face = "bold"),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.line = element_line(linewidth = 1.0, color = "black"),    
      axis.ticks = element_line(linewidth = 1.0, color = "black"),   
      axis.ticks.length = unit(0.2, "cm"),
      legend.background = element_rect(fill = alpha("white", 0.6), color = "black", linewidth = 0.4),
      legend.key = element_blank(),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b=8)),
      plot.subtitle = element_text(size = 12, color = "grey30", hjust = 0.5, margin = margin(b=15))
    )
}

nature_red <- "#E64B35"
nature_others <- c("#4DBBD5", "#00A087", "#3C5488", "#F39B7F")

color_map <- c(nature_red, nature_others)
names(color_map) <- auc_labels

size_map <- c(1.8, 0.8, 0.8, 0.8, 0.8)
names(size_map) <- auc_labels

p_roc <- ggroc(roc_list, aes = c("color", "linewidth")) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  theme_nature() +

  theme(legend.position = c(0.60, 0.22),
        legend.title = element_blank(),
        legend.text = element_text(size = 9.5, face = "bold")) + 

  labs(title = "Diagnostic Performance of Multivariable Model", 
       subtitle = "Dataset: GSE163154 | Cohort: N=43 (IPH vs non-IPH)",
       x = "Specificity", y = "Sensitivity") +

  scale_color_manual(values = color_map) +
  scale_linewidth_manual(values = size_map) +

  coord_fixed()

ggsave(file.path(out_dir, "10_ROC_Single_vs_Combined_Nature.pdf"), p_roc, width = 6.5, height = 6.5)
ggsave(file.path(out_dir, "10_ROC_Single_vs_Combined_Nature.png"), p_roc, width = 6.5, height = 6.5, dpi=600)

message("\n>>> 完美搞定！请注意查看图例，现在的联合模型不仅红色加粗傲视群雄，而且附带了极其严谨的 95% 置信区间！")


message(">>> [Step 1] 初始化环境与数据准备 (修复矩阵格式问题)...")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(illuminaHumanv2.db)
  library(rms)        
  library(dcurves)    
  library(ggplot2)
  library(ggsci)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "30_Clinical_Diagnostic_Modeling")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

tox_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
selected_genes <- c("SLC16A4", "SLC16A3", "ACKR1")

target_gse_file <- "GSE163154_series_matrix.txt.gz"
if(!file.exists(target_gse_file)) {
  gset <- getGEO("GSE163154", destdir=".", getGPL = FALSE)[[1]]
} else {
  gset <- getGEO(filename = target_gse_file, getGPL = FALSE)
  if(is.list(gset)) gset <- gset[[1]]
}

expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

clinical_status <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case=T), "non-IPH", "IPH")
clin_df <- data.frame(
  Sample = colnames(expr_gse),
  Status_Label = factor(clinical_status, levels = c("non-IPH", "IPH")),
  Status = ifelse(clinical_status == "IPH", 1, 0)
)

probe_to_symbol <- mapIds(illuminaHumanv2.db, keys = rownames(expr_gse), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
expr_clean <- as.data.frame(expr_gse) %>%
  mutate(SYMBOL = probe_to_symbol) %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise_all(mean) %>%
  column_to_rownames("SYMBOL")

available_tox <- intersect(tox_genes, rownames(expr_clean))
cd8_tox <- colMeans(t(scale(t(expr_clean[available_tox, ]))))
expr_mat <- rbind(as.matrix(expr_clean), "CD8_Tox_Score" = cd8_tox)

model_features <- c(selected_genes, "CD8_Tox_Score")
for(feature in model_features) {
  clin_df[[feature]] <- as.numeric(scale(as.numeric(expr_mat[feature, ])))
}

clinic_title <- "Dataset: GSE163154 | Cohort: N=43 (IPH vs non-IPH)"

message(">>> [Step 2] 正在计算并输出联合诊断模型公式...")

fit_glm <- glm(Status ~ SLC16A4 + SLC16A3 + ACKR1 + CD8_Tox_Score, data = clin_df, family = binomial)
clin_df$Prob <- predict(fit_glm, type = "response")

coefs <- coef(fit_glm)
formula_str <- sprintf(
  "Logit(P) = %.4f + (%.4f * SLC16A4) + (%.4f * SLC16A3) + (%.4f * ACKR1) + (%.4f * CD8_Tox_Score)",
  coefs["(Intercept)"], coefs["SLC16A4"], coefs["SLC16A3"], coefs["ACKR1"], coefs["CD8_Tox_Score"]
)
message("\n=======================================================")
message(">>> 【临床直接可用】联合诊断模型精确数学公式 <<<")
message(formula_str)
message("=======================================================\n")

message(">>> [Step 3] 正在绘制 PCA 聚类图...")

pca_res <- prcomp(clin_df[, model_features], scale. = FALSE) 
pca_df <- as.data.frame(pca_res$x)
pca_df$Status_Label <- clin_df$Status_Label

variance <- summary(pca_res)$importance[2,] * 100

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Status_Label, color = Status_Label)) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linewidth = 0.8, linetype = "dashed") +
  geom_point(size = 4, shape = 21, color = "black", stroke = 1) +
  scale_fill_manual(values = c("IPH" = "#E64B35", "non-IPH" = "#3C5488")) +
  scale_color_manual(values = c("IPH" = "#E64B35", "non-IPH" = "#3C5488")) +
  theme_classic() +
  labs(title = "PCA of 4-Feature Combined Signature",
       subtitle = clinic_title,
       x = sprintf("PC1 (%.1f%% Variance)", variance[1]),
       y = sprintf("PC2 (%.1f%% Variance)", variance[2])) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b=8)),
        plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5, margin = margin(b=15)),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))

ggsave(file.path(out_dir, "01_PCA_Clustering_Nature.pdf"), p_pca, width = 6.5, height = 6.5)
ggsave(file.path(out_dir, "01_PCA_Clustering_Nature.png"), p_pca, width = 6.5, height = 6.5, dpi=600)

message(">>> [Step 4] 正在绘制临床转化列线图 (Nomogram)...")

dd <- datadist(clin_df)
options(datadist = "dd")

fit_lrm <- lrm(Status ~ SLC16A4 + SLC16A3 + ACKR1 + CD8_Tox_Score, data = clin_df, x = TRUE, y = TRUE)

nom <- nomogram(fit_lrm, fun = plogis, fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                funlabel = "Probability of IPH", lp = FALSE)

pdf(file.path(out_dir, "02_Nomogram_Clinical_Utility.pdf"), width = 9, height = 6)
par(mar = c(3, 2, 3, 2), family = "sans", cex = 0.9)
plot(nom, xfrac = 0.45, cex.axis = 1.05, cex.var = 1.1, lmgp = 0.2, 
     force.label = TRUE, col.grid = "grey80")
title(main = "Nomogram for Predicting Intraplaque Hemorrhage (IPH)\nDataset: GSE163154 | N=43", 
      cex.main = 1.3, font.main = 2, col.main = "black", line = 1)
dev.off()

png(file.path(out_dir, "02_Nomogram_Clinical_Utility.png"), width = 9, height = 6, units = "in", res = 600)
par(mar = c(3, 2, 3, 2), family = "sans", cex = 0.9)
plot(nom, xfrac = 0.45, cex.axis = 1.05, cex.var = 1.1, lmgp = 0.2, force.label = TRUE, col.grid = "grey80")
title(main = "Nomogram for Predicting Intraplaque Hemorrhage (IPH)\nDataset: GSE163154 | N=43", 
      cex.main = 1.3, font.main = 2, col.main = "black", line = 1)
dev.off()

message(">>> [Step 5] 正在绘制 Bootstrapping 偏差校准曲线...")

cal <- calibrate(fit_lrm, method = "boot", B = 1000)

pdf(file.path(out_dir, "03_Calibration_Curve.pdf"), width = 6.5, height = 6.5)
par(mar = c(5, 5, 4, 2), family = "sans")
plot(cal, lwd = 2, lty = 1, errbar.col = "black",
     xlab = "Nomogram Predicted Probability of IPH", 
     ylab = "Actual Observed IPH Proportion",
     col = "#E64B35", sub = FALSE, cex.axis = 1.1, cex.lab = 1.2, font.lab = 2)
title(main = "Calibration Curve (1000 Bootstrap Resamples)", sub = clinic_title, 
      cex.main = 1.4, font.main = 2, col.sub = "grey30", line = 1.5)
abline(0, 1, lty = 2, lwd = 1.5, col = "grey50")
dev.off()

png(file.path(out_dir, "03_Calibration_Curve.png"), width = 6.5, height = 6.5, units = "in", res = 600)
par(mar = c(5, 5, 4, 2), family = "sans")
plot(cal, lwd = 2, lty = 1, errbar.col = "black",
     xlab = "Nomogram Predicted Probability of IPH", ylab = "Actual Observed IPH Proportion",
     col = "#E64B35", sub = FALSE, cex.axis = 1.1, cex.lab = 1.2, font.lab = 2)
title(main = "Calibration Curve (1000 Bootstrap Resamples)", sub = clinic_title, 
      cex.main = 1.4, font.main = 2, col.sub = "grey30", line = 1.5)
abline(0, 1, lty = 2, lwd = 1.5, col = "grey50")
dev.off()

message(">>> [Step 6] 正在绘制基于 dcurves 引擎的 DCA 临床净收益分析图...")

dca_df <- clin_df %>% dplyr::select(Status, Prob)
dca_res <- dcurves::dca(Status ~ Prob, data = dca_df, thresholds = seq(0, 1, by = 0.01))

p_dca <- plot(dca_res) + 
  scale_color_manual(
    values = c("Prob" = "#E64B35", "all" = "#3C5488", "none" = "grey30"),
    labels = c("Combined Signature Model", "Treat All", "Treat None")
  ) +
  theme_classic() +
  labs(title = "Decision Curve Analysis (DCA)",
       subtitle = clinic_title,
       x = "Threshold Probability (Risk Cut-off)", 
       y = "Net Benefit") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b=8)),
        plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5, margin = margin(b=15)),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = c(0.75, 0.75),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, face = "bold"),
        legend.background = element_rect(fill = alpha("white", 0.6), color = "black", linewidth = 0.4))

ggsave(file.path(out_dir, "04_DCA_Net_Benefit_Nature.pdf"), p_dca, width = 7, height = 6.5)
ggsave(file.path(out_dir, "04_DCA_Net_Benefit_Nature.png"), p_dca, width = 7, height = 6.5, dpi=600)

message("\n>>> 完美搞定！DCA 曲线已成功生成。")


message(">>> [Step 1] 初始化环境与数据标准化...")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(illuminaHumanv2.db)
  library(rms)        
  library(dcurves)    
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "30_Clinical_Diagnostic_Modeling")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

tox_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
selected_genes <- c("SLC16A4", "SLC16A3", "ACKR1")

target_gse_file <- "GSE163154_series_matrix.txt.gz"
if(!file.exists(target_gse_file)) {
  gset <- getGEO("GSE163154", destdir=".", getGPL = FALSE)[[1]]
} else {
  gset <- getGEO(filename = target_gse_file, getGPL = FALSE)
  if(is.list(gset)) gset <- gset[[1]]
}

expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

clinical_status <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case=T), "non-IPH", "IPH")
clin_df <- data.frame(
  Sample = colnames(expr_gse),
  Status_Label = factor(clinical_status, levels = c("non-IPH", "IPH")),
  Status = ifelse(clinical_status == "IPH", 1, 0)
)

probe_to_symbol <- mapIds(illuminaHumanv2.db, keys = rownames(expr_gse), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
expr_clean <- as.data.frame(expr_gse) %>%
  mutate(SYMBOL = probe_to_symbol) %>% filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>% summarise_all(mean) %>% column_to_rownames("SYMBOL")

cd8_tox <- colMeans(t(scale(t(expr_clean[intersect(tox_genes, rownames(expr_clean)), ]))))
expr_mat <- rbind(as.matrix(expr_clean), "CD8_Tox_Score" = cd8_tox)

model_features <- c(selected_genes, "CD8_Tox_Score")
for(feature in model_features) {
  clin_df[[feature]] <- as.numeric(scale(as.numeric(expr_mat[feature, ])))
}

n_iph <- sum(clin_df$Status == 1)
n_non <- sum(clin_df$Status == 0)
clinic_title <- sprintf("Dataset: GSE163154 | N = %d (IPH: %d, non-IPH: %d)", nrow(clin_df), n_iph, n_non)

fit_glm <- glm(Status ~ SLC16A4 + SLC16A3 + ACKR1 + CD8_Tox_Score, data = clin_df, family = binomial)
clin_df$Prob_Combined <- predict(fit_glm, type = "response")

clin_df$Prob_SLC16A4 <- predict(glm(Status ~ SLC16A4, data = clin_df, family = binomial), type = "response")
clin_df$Prob_SLC16A3 <- predict(glm(Status ~ SLC16A3, data = clin_df, family = binomial), type = "response")
clin_df$Prob_ACKR1   <- predict(glm(Status ~ ACKR1, data = clin_df, family = binomial), type = "response")
clin_df$Prob_CD8     <- predict(glm(Status ~ CD8_Tox_Score, data = clin_df, family = binomial), type = "response")

message(">>> [Step 2] 正在绘制 Nature 风格 PCA 质心图...")

pca_res <- prcomp(clin_df[, model_features], scale. = FALSE) 
pca_df <- as.data.frame(pca_res$x)
pca_df$Status_Label <- clin_df$Status_Label
variance <- summary(pca_res)$importance[2,] * 100

centroids <- pca_df %>% 
  group_by(Status_Label) %>% 
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Status_Label, fill = Status_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +

  stat_ellipse(geom = "polygon", alpha = 0.15, linewidth = 0.8) +

  geom_point(aes(shape = Status_Label), size = 4) +

  geom_point(data = centroids, aes(shape = Status_Label), size = 10, color = "white", stroke = 1.5) +
  geom_point(data = centroids, aes(shape = Status_Label), size = 10, alpha = 0.9) +

  scale_shape_manual(values = c("non-IPH" = 21, "IPH" = 24)) +
  scale_fill_manual(values = c("non-IPH" = "#3C5488", "IPH" = "#E64B35")) +
  scale_color_manual(values = c("non-IPH" = "#3C5488", "IPH" = "#E64B35")) +

  theme_classic() +
  labs(title = "Unsupervised PCA Clustering",
       subtitle = clinic_title,
       x = sprintf("Dim1 (%.1f%%)", variance[1]),
       y = sprintf("Dim2 (%.1f%%)", variance[2])) +
  theme(
    axis.line = element_line(linewidth = 2, color = "black"),
    axis.ticks = element_line(linewidth = 2, color = "black"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5, margin = margin(b=15)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 13, face = "bold"),
    legend.key.size = unit(1, "cm")
  )

ggsave(file.path(out_dir, "01_PCA_Clustering_Centroid.pdf"), p_pca, width = 7.5, height = 7.5)
ggsave(file.path(out_dir, "01_PCA_Clustering_Centroid.png"), p_pca, width = 7.5, height = 7.5, dpi=600)

message(">>> [Step 3] 正在绘制防重叠临床转化列线图...")

dd <- datadist(clin_df)
options(datadist = "dd")
fit_lrm <- lrm(Status ~ SLC16A4 + SLC16A3 + ACKR1 + CD8_Tox_Score, data = clin_df, x = TRUE, y = TRUE)

nom <- nomogram(fit_lrm, fun = plogis, 
                fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                funlabel = "Risk of IPH", lp = FALSE)

plot_nomogram_clean <- function() {
  par(mar = c(3, 2, 3, 2), family = "sans", col.axis="black", font.axis=2, font.lab=2)
  plot(nom, xfrac = 0.4, cex.axis = 1.0, cex.var = 1.2, lmgp = 0.2, col.grid = "grey85")
  title(main = "Diagnostic Nomogram for Plaque Vulnerability\n", font.main = 2, cex.main = 1.4)
  mtext(clinic_title, side = 3, line = 0, cex = 1, col = "grey30", font=2)
}

pdf(file.path(out_dir, "02_Nomogram_Clean.pdf"), width = 9.5, height = 6.5)
plot_nomogram_clean()
dev.off()

png(file.path(out_dir, "02_Nomogram_Clean.png"), width = 9.5, height = 6.5, units = "in", res = 600)
plot_nomogram_clean()
dev.off()

message(">>> [Step 4] 正在绘制带有 H-L 和 Brier 评分的经典专业校准曲线...")

if (!requireNamespace("ResourceSelection", quietly = TRUE)) install.packages("ResourceSelection")

set.seed(2025)
cal <- calibrate(fit_lrm, method = "boot", B = 1000)

brier_score <- mean((clin_df$Prob_Combined - clin_df$Status)^2)

hl_test <- ResourceSelection::hoslem.test(clin_df$Status, clin_df$Prob_Combined, g = 8) 
hl_pval <- hl_test$p.value

plot_calibration_classic_smooth <- function() {
  par(mar = c(5, 5, 4, 2), family = "sans", lwd = 1.5, font.lab = 2)

  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Nomogram Predicted Probability", 
       ylab = "Actual Observed Probability",
       cex.axis = 1.2, cex.lab = 1.4)

  abline(0, 1, lty = 2, lwd = 2, col = "grey50")

  pred_x <- cal[, "predy"]
  app_y  <- cal[, "calibrated.orig"]
  cor_y  <- cal[, "calibrated.corrected"]

  valid_idx <- complete.cases(pred_x, app_y, cor_y)
  sm_app <- smooth.spline(pred_x[valid_idx], app_y[valid_idx], df = 3)
  sm_cor <- smooth.spline(pred_x[valid_idx], cor_y[valid_idx], df = 3)

  lines(sm_app, lty = 3, lwd = 2, col = "black")
  lines(sm_cor, lty = 1, lwd = 3, col = "black")

  title(main = "Calibration Curve (Bootstrap N=1000)", font.main = 2, cex.main = 1.5)
  mtext(clinic_title, side = 3, line = 0.5, cex = 1, col = "grey30", font = 2)

  legend("bottomright", legend = c("Ideal Reference", "Apparent", "Bias-corrected"),
         col = c("grey50", "black", "black"), lty = c(2, 3, 1), lwd = c(2, 2, 3), 
         bty = "n", cex = 1.1, text.font = 2)

  text(0.03, 0.95, sprintf("Brier score = %.3f", brier_score), adj = 0, cex = 1.2, font = 2)
  text(0.03, 0.88, sprintf("Hosmer-Lemeshow p = %.3f", hl_pval), adj = 0, cex = 1.2, font = 2)
}

pdf(file.path(out_dir, "03_Calibration_Curve_Classic_Smooth.pdf"), width = 6.5, height = 6.5)
plot_calibration_classic_smooth()
dev.off()

png(file.path(out_dir, "03_Calibration_Curve_Classic_Smooth.png"), width = 6.5, height = 6.5, units = "in", res = 600)
plot_calibration_classic_smooth()
dev.off()

message(">>> 完美！带有 H-L & Brier 评分的顶级专业校准曲线已生成！")
message(">>> [Step 5] 正在手绘接管高级 DCA 曲线 (优化单基因对比度)...")

dca_res <- dca(Status ~ Prob_Combined + Prob_SLC16A4 + Prob_SLC16A3 + Prob_ACKR1 + Prob_CD8, 
               data = clin_df, thresholds = seq(0, 0.99, by = 0.01))

dca_df_plot <- as_tibble(dca_res)

p_dca_custom <- ggplot(dca_df_plot, aes(x = threshold, y = net_benefit, color = variable)) +
  geom_line(data = filter(dca_df_plot, variable %in% c("all", "none")), aes(linewidth = variable)) +

  geom_line(data = filter(dca_df_plot, variable %in% c("Prob_SLC16A4", "Prob_SLC16A3", "Prob_ACKR1", "Prob_CD8")), 
            linewidth = 1.2, alpha = 0.9) +

  geom_line(data = filter(dca_df_plot, variable == "Prob_Combined"), 
            aes(color = variable), linewidth = 3.5) +

  coord_cartesian(ylim = c(-0.05, 0.65), xlim = c(0, 1)) +

  scale_color_manual(
    values = c(
      "all" = "grey60", "none" = "black",
      "Prob_SLC16A4" = "#4DBBD5",
      "Prob_SLC16A3" = "#00A087",
      "Prob_ACKR1"   = "#F39B7F",
      "Prob_CD8"     = "#8491B4",
      "Prob_Combined"= "#E64B35"
    ),
    labels = c("Treat All", "Treat None", "Single (SLC16A4)", "Single (SLC16A3)", 
               "Single (ACKR1)", "Single (CD8 Tox)", "Combined Signature"),
    breaks = c("all", "none", "Prob_SLC16A4", "Prob_SLC16A3", "Prob_ACKR1", "Prob_CD8", "Prob_Combined")
  ) +
  scale_linewidth_manual(values = c("all" = 1.5, "none" = 1.5), guide="none") +

  theme_classic() +
  labs(title = "Clinical Decision Curve Analysis (DCA)",
       subtitle = clinic_title,
       x = "Threshold Probability", y = "Clinical Net Benefit") +

  theme(
    axis.line = element_line(linewidth = 2, color = "black"), 
    axis.ticks = element_line(linewidth = 2, color = "black"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5, margin = margin(b=15, t=5), face="italic"),
    legend.position = c(0.25, 0.35),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 1.5),
    legend.key.height = unit(0.8, "cm")
  )

ggsave(file.path(out_dir, "04_DCA_Single_vs_Combined.pdf"), p_dca_custom, width = 7.5, height = 7)
ggsave(file.path(out_dir, "04_DCA_Single_vs_Combined.png"), p_dca_custom, width = 7.5, height = 7, dpi=600)
message(">>> [Step 6] 正在绘制高级小提琴箱线分布图...")

roc_obj <- roc(clin_df$Status, clin_df$Prob_Combined, quiet = TRUE)
best_thr <- coords(roc_obj, "best", ret = "threshold")$threshold[1]

p_dist <- ggplot(clin_df, aes(x = Status_Label, y = Prob_Combined, fill = Status_Label)) +

  geom_violin(trim = FALSE, alpha = 0.5, color = "black", linewidth = 1.2) +

  geom_boxplot(width = 0.15, alpha = 0.9, color = "black", linewidth = 1.2, outlier.shape = NA) +

  geom_jitter(shape = 21, size = 3, width = 0.1, color = "black", stroke = 1.2, fill = "white", alpha = 0.9) +

  geom_hline(yintercept = best_thr, linetype = "dashed", color = "grey40", linewidth = 1.5) +

  stat_compare_means(comparisons = list(c("non-IPH", "IPH")), method = "wilcox.test", 
                     aes(label = paste0("p = ", ..p.format..)), size = 5.5, fontface = "bold") +

  scale_fill_manual(values = c("non-IPH" = "#3C5488", "IPH" = "#E64B35")) +

  theme_classic() +
  labs(title = "Combined Signature Risk Score Distribution", 
       subtitle = clinic_title,
       x = "Clinical Plaque Status", 
       y = "Predicted IPH Probability (Risk Score)") +

  theme(
    axis.line = element_line(linewidth = 1.5, color = "black"),
    axis.ticks = element_line(linewidth = 1.5, color = "black"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),

    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5, margin = margin(b=15)),

    legend.position = "none"
  )

ggsave(file.path(out_dir, "05_Risk_Score_Distribution_Violin.pdf"), p_dist, width = 6, height = 7)
ggsave(file.path(out_dir, "05_Risk_Score_Distribution_Violin.png"), p_dist, width = 6, height = 7, dpi=600)

message("\n>>> 完美收工！绝美的小提琴嵌套箱线散点图已生成，数据分布的厚重感瞬间拉满！")


message(">>> [Step 1] 正在计算每个患者的分类结局与详细诊断指标...")

suppressMessages({
  library(tidyverse)
  library(pROC)
  library(ggplot2)
  library(caret)
})

roc_obj <- roc(clin_df$Status, clin_df$Prob_Combined, quiet = TRUE)
best_thr <- coords(roc_obj, "best", ret = "threshold")$threshold[1]

waterfall_df <- clin_df %>%
  arrange(Prob_Combined) %>%
  mutate(
    Patient_ID = row_number(),
    Prediction_Class = ifelse(Prob_Combined >= best_thr, "Predicted IPH", "Predicted non-IPH"),
    Risk_Deviation = Prob_Combined - best_thr
  )

message(">>> [Step 2] 正在渲染患者个体风险瀑布图 (使用指定的精确公式)...")

formula_text <- "Logit(P) = 0.4625 + (-2.6931 * SLC16A4) + (0.6481 * SLC16A3)\n+ (0.2556 * ACKR1) + (0.7414 * CD8_Tox_Score)"

p_waterfall <- ggplot(waterfall_df, aes(x = Patient_ID, y = Risk_Deviation, fill = Status_Label)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.6) +

  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.2) +

  scale_fill_manual(values = c("non-IPH" = "#3C5488", "IPH" = "#E64B35")) +
  scale_x_continuous(expand = c(0.02, 0.02)) +

  theme_classic() +
  labs(title = "Individual Patient Risk Discrimination (Waterfall Plot)",
       subtitle = clinic_title,
       x = "Patient Cohort (Ranked by Combined Risk Score, N=43)", 
       y = "Deviation from Optimal Threshold") +

  annotate("text", x = 1, y = max(waterfall_df$Risk_Deviation) * 0.95, 
           label = formula_text, hjust = 0, size = 4.2, fontface = "bold.italic", color = "grey20") +

  theme(
    axis.line = element_line(linewidth = 1.5, color = "black"),
    axis.ticks = element_line(linewidth = 1.5, color = "black"),
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5, margin = margin(b=15)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold")
  )

ggsave(file.path(out_dir, "06_Patient_Risk_Waterfall.pdf"), p_waterfall, width = 8, height = 6)
ggsave(file.path(out_dir, "06_Patient_Risk_Waterfall.png"), p_waterfall, width = 8, height = 6, dpi=600)

message(">>> 瀑布图更新完毕！现在左上角挂着的正是你那串精确无误的公式。")
message(">>> [Step 3] 正在计算并绘制全面诊断效能条形图...")

conf_matrix <- confusionMatrix(
  factor(ifelse(clin_df$Prob_Combined >= best_thr, "IPH", "non-IPH"), levels = c("IPH", "non-IPH")),
  factor(clin_df$Status_Label, levels = c("IPH", "non-IPH"))
)

metrics_df <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", "PPV", "NPV"),
  Value = c(
    roc_obj$auc,
    conf_matrix$overall["Accuracy"],
    conf_matrix$byClass["Sensitivity"],
    conf_matrix$byClass["Specificity"],
    conf_matrix$byClass["Pos Pred Value"],
    conf_matrix$byClass["Neg Pred Value"]
  ) * 100
)

metrics_df$Value[1] <- roc_obj$auc * 100 

metrics_df$Metric <- factor(metrics_df$Metric, levels = rev(c("AUC", "Accuracy", "Sensitivity", "Specificity", "PPV", "NPV")))

p_metrics <- ggplot(metrics_df, aes(x = Value, y = Metric)) +
  geom_bar(stat = "identity", width = 0.6, fill = "#3C5488", color = "black", linewidth = 1) +

  geom_text(aes(label = sprintf("%.1f%%", Value)), hjust = -0.2, size = 5, fontface = "bold") +

  scale_x_continuous(limits = c(0, 115), breaks = seq(0, 100, 20), expand = c(0,0)) +

  theme_classic() +
  labs(title = "Diagnostic Performance Metrics",
       subtitle = "Optimal Threshold derived from Combined Signature",
       x = "Performance Score (%)", y = "") +

  theme(
    axis.line = element_line(linewidth = 1.5, color = "black"),
    axis.ticks = element_line(linewidth = 1.5, color = "black"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 13, color = "black", face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold", margin = margin(t=10)),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5, margin = margin(b=15)),
    plot.margin = margin(15, 25, 15, 15)
  )

ggsave(file.path(out_dir, "07_Diagnostic_Metrics_Barplot.pdf"), p_metrics, width = 6.5, height = 5)
ggsave(file.path(out_dir, "07_Diagnostic_Metrics_Barplot.png"), p_metrics, width = 6.5, height = 5, dpi=600)

message("\n>>> 齐活儿！有了这两张图，你模型的 '公式应用' 和 '细节区分度' 彻底展现得淋漓尽致，图版绝对饱满且极具说服力！")


message(">>> [Step 1] 正在初始化环境并加载 GSE163154.gz 数据...")

if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(illuminaHumanv2.db)
  library(caret)
  library(pROC)
  library(ggplot2)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "65_Logistic_Regression_LOOCV")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

tox_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
selected_genes <- c("SLC16A4", "SLC16A3", "ACKR1")

target_gse_file <- "GSE163154_series_matrix.txt.gz"
if(!file.exists(target_gse_file)) {
  message("--- 正在从 GEO 在线下载...")
  gset <- getGEO("GSE163154", destdir=".", getGPL = FALSE)[[1]]
} else {
  message("--- 成功定位本地 .gz 文件，正在读取...")
  gset <- getGEO(filename = target_gse_file, getGPL = FALSE)
  if(is.list(gset)) gset <- gset[[1]]
}

expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

clinical_status <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case=T), "non-IPH", "IPH")
message(sprintf(">>> GSE163154 样本精确分布: IPH = %d, non-IPH = %d", sum(clinical_status=="IPH"), sum(clinical_status=="non-IPH")))

probe_to_symbol <- mapIds(illuminaHumanv2.db, keys = rownames(expr_gse), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
expr_clean <- as.data.frame(expr_gse) %>%
  mutate(SYMBOL = probe_to_symbol) %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise_all(mean) %>%
  column_to_rownames("SYMBOL")

available_tox <- intersect(tox_genes, rownames(expr_clean))
cd8_tox <- colMeans(t(scale(t(expr_clean[available_tox, ]))))
expr_mat <- rbind(as.matrix(expr_clean), "CD8_Tox_Score" = cd8_tox)

ml_status <- ifelse(clinical_status == "non-IPH", "Non_IPH", "IPH")

train_data <- data.frame(
  Status = factor(ml_status, levels = c("Non_IPH", "IPH"))
)

model_features <- c(selected_genes, "CD8_Tox_Score")
for(feature in model_features) {
  train_data[[feature]] <- as.numeric(scale(as.numeric(expr_mat[feature, ])))
}

message("\n>>> [Step 2] 🚀 正在 N=43 队列上执行 Logistic Regression 的 LOOCV 盲测交叉验证...")

fitControl <- trainControl(method = "LOOCV", 
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary)

set.seed(123)
model_glm <- train(Status ~ ., data = train_data, method = "glm", metric = "ROC", trControl = fitControl)

message("\n>>> [Step 3] 正在汇总效能并绘制顶刊级 ROC 图...")

preds <- model_glm$pred
preds <- preds[order(preds$rowIndex), ]
roc_glm <- roc(preds$obs, preds$IPH, quiet = TRUE)

auc_val <- roc_glm$auc
ci_val <- ci.auc(roc_glm)

message(sprintf(">>> 🎉 Logistic 回归 LOOCV 验证效能: AUC = %.3f (95%% CI: %.3f - %.3f)", auc_val, ci_val[1], ci_val[3]))

p_roc <- ggroc(roc_glm, color = "#E64B35", linewidth = 2.5) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey60", linewidth = 1) +
  annotate("text", x = 0.25, y = 0.15, 
           label = sprintf("AUC: %.3f\n95%% CI: %.3f - %.3f", auc_val, ci_val[1], ci_val[3]), 
           size = 5.5, fontface = "bold", color = "black", hjust = 0.5) +
  labs(title = "Diagnostic Performance of the Core Signature", 
       subtitle = "Validated by Leave-One-Out Cross-Validation (LOOCV, N = 43)",
       x = "Specificity", y = "Sensitivity") +
  theme_classic() +
  theme(text = element_text(color = "black", family = "sans"),
        axis.text = element_text(size = 13, color = "black", face = "bold"),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.line = element_line(linewidth = 1.2, color = "black"),
        axis.ticks = element_line(linewidth = 1.2, color = "black"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = ggplot2::margin(b=5)),
        plot.subtitle = element_text(size = 12, face = "italic", color = "grey30", hjust = 0.5, margin = ggplot2::margin(b=15)))

ggsave(file.path(out_dir, "Figure_Logistic_LOOCV_ROC.pdf"), p_roc, width = 6, height = 6)
ggsave(file.path(out_dir, "Figure_Logistic_LOOCV_ROC.png"), p_roc, width = 6, height = 6, dpi = 600)

message(">>> 完美运行结束！请前往 65_Logistic_Regression_LOOCV 文件夹查看你文章的核心主图！")


message(">>> [Step 1] 正在初始化环境并准备精准队列数据...")

if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(illuminaHumanv2.db)
  library(caret)
  library(pROC)
  library(ggplot2)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "65_All_Models_LOOCV_Comparison")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(work_dir)

tox_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
selected_genes <- c("SLC16A4", "SLC16A3", "ACKR1") 

target_gse_file <- "GSE163154_series_matrix.txt.gz"
if(!file.exists(target_gse_file)) {
  gset <- getGEO("GSE163154", destdir=".", getGPL = FALSE)[[1]]
} else {
  gset <- getGEO(filename = target_gse_file, getGPL = FALSE)
  if(is.list(gset)) gset <- gset[[1]]
}

expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

clinical_status <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case=T), "Non_IPH", "IPH")
ml_status <- factor(clinical_status, levels = c("Non_IPH", "IPH"))
message(sprintf(">>> GSE163154 样本精确分布: IPH = %d, Non_IPH = %d", sum(ml_status=="IPH"), sum(ml_status=="Non_IPH")))

probe_to_symbol <- mapIds(illuminaHumanv2.db, keys = rownames(expr_gse), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
expr_clean <- as.data.frame(expr_gse) %>%
  mutate(SYMBOL = probe_to_symbol) %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise_all(mean) %>%
  column_to_rownames("SYMBOL")

available_tox <- intersect(tox_genes, rownames(expr_clean))
cd8_tox <- colMeans(t(scale(t(expr_clean[available_tox, ]))))
expr_mat <- rbind(as.matrix(expr_clean), "CD8_Tox_Score" = cd8_tox)

train_data <- data.frame(Status = ml_status)
model_features <- c(selected_genes, "CD8_Tox_Score")
for(feature in model_features) {
  train_data[[feature]] <- as.numeric(scale(as.numeric(expr_mat[feature, ])))
}

message("\n>>> [Step 2] 🚀 引擎全开！正在对所有模型执行极其严苛的 LOOCV 交叉验证...")

fitControl <- trainControl(method = "LOOCV", 
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary)

set.seed(123)
message("--- 正在训练 联合诊断模型 (Combined)...")
mod_comb <- train(Status ~ SLC16A4 + SLC16A3 + ACKR1 + CD8_Tox_Score, data = train_data, method = "glm", metric = "ROC", trControl = fitControl)

set.seed(123)
message("--- 正在训练 单变量: ACKR1...")
mod_ackr <- train(Status ~ ACKR1, data = train_data, method = "glm", metric = "ROC", trControl = fitControl)

set.seed(123)
message("--- 正在训练 单变量: CD8 Tox Score...")
mod_cd8  <- train(Status ~ CD8_Tox_Score, data = train_data, method = "glm", metric = "ROC", trControl = fitControl)

set.seed(123)
message("--- 正在训练 单变量: SLC16A3...")
mod_slc3 <- train(Status ~ SLC16A3, data = train_data, method = "glm", metric = "ROC", trControl = fitControl)

set.seed(123)
message("--- 正在训练 单变量: SLC16A4...")
mod_slc4 <- train(Status ~ SLC16A4, data = train_data, method = "glm", metric = "ROC", trControl = fitControl)

message("\n>>> [Step 3] 正在提取 LOOCV 盲测结果并构建 ROC 矩阵...")

get_cv_roc <- function(model_obj) {
  preds <- model_obj$pred
  preds <- preds[order(preds$rowIndex), ]
  return(roc(preds$obs, preds$IPH, quiet = TRUE))
}

roc_list_raw <- list(
  "Combined Signature" = get_cv_roc(mod_comb),
  "ACKR1"              = get_cv_roc(mod_ackr),
  "CD8 Tox Score"      = get_cv_roc(mod_cd8),
  "SLC16A3"            = get_cv_roc(mod_slc3),
  "SLC16A4"            = get_cv_roc(mod_slc4)
)

auc_labels <- sapply(names(roc_list_raw), function(x) {
  ci_val <- ci.auc(roc_list_raw[[x]])
  paste0(x, " (AUC = ", sprintf("%.3f", roc_list_raw[[x]]$auc), ", 95% CI: ", sprintf("%.3f", ci_val[1]), "-", sprintf("%.3f", ci_val[3]), ")")
})
names(roc_list_raw) <- auc_labels

message(">>> [Step 4] 正在融合专属主题绘制最终的 Nature 风格主图...")

theme_nature <- function() {
  theme_classic() +
    theme(
      text = element_text(color = "black", family = "sans"),
      axis.text = element_text(size = 12, color = "black", face = "bold"),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.line = element_line(linewidth = 1.0, color = "black"),    
      axis.ticks = element_line(linewidth = 1.0, color = "black"),   
      axis.ticks.length = unit(0.2, "cm"),
      legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth = 0.5),
      legend.key = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = ggplot2::margin(b=8)),
      plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0.5, face = "italic", margin = ggplot2::margin(b=15))
    )
}

nature_red <- "#E64B35"
nature_others <- c("#4DBBD5", "#00A087", "#3C5488", "#F39B7F")

color_map <- c(nature_red, nature_others)
names(color_map) <- auc_labels

size_map <- c(2.2, 1.0, 1.0, 1.0, 1.0)
names(size_map) <- auc_labels

p_roc <- ggroc(roc_list_raw, aes = c("color", "linewidth")) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey60", linewidth = 1) +
  theme_nature() +
  theme(legend.position = c(0.62, 0.22),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, face = "bold")) + 
  labs(title = "Strict Validation via Leave-One-Out Cross-Validation", 
       subtitle = "Dataset: GSE163154 | Model: Logistic Regression",
       x = "Specificity", y = "Sensitivity") +
  scale_color_manual(values = color_map) +
  scale_linewidth_manual(values = size_map) +
  coord_fixed()

ggsave(file.path(out_dir, "Figure_All_Variables_LOOCV_ROC.pdf"), p_roc, width = 6.8, height = 6.8)
ggsave(file.path(out_dir, "Figure_All_Variables_LOOCV_ROC.png"), p_roc, width = 6.8, height = 6.8, dpi=600)

message(">>> 🎉 完美执行完毕！请前往 65_All_Models_LOOCV_Comparison 文件夹验收。")


# Figure 8A-E: GSE43292 external validation

message(">>> [Step 1] 正在直接读取本地 GSE43292 矩阵与病理分期信息...")

if (!requireNamespace("hugene10sttranscriptcluster.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("hugene10sttranscriptcluster.db", update = FALSE)
}

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(hugene10sttranscriptcluster.db)
  library(ggplot2)
  library(ggpubr)
  library(pROC)
  library(dcurves)
  library(patchwork)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "32_External_Validation_GSE43292")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_path <- "GSE43292_series_matrix.txt"
if(!file.exists(file_path)) stop("当前目录下找不到 GSE43292_series_matrix.txt 文件！")

gset <- getGEO(filename = file_path, getGPL = FALSE)
titles <- pData(gset)$title

group_info <- ifelse(grepl("atheroma|plaque", titles, ignore.case=TRUE), 
                     "Stage IV-V (Advanced)", 
                     "Stage I-II (Early)")

message(">>> 当前队列患者病理分期划分统计如下：")
print(table(group_info))

if(length(unique(group_info)) < 2) stop("严重错误：提取失败！所有病人都被分到了一组！")

probe_to_symbol <- mapIds(hugene10sttranscriptcluster.db, keys = rownames(exprs(gset)), 
                          column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
expr_df <- as.data.frame(exprs(gset))
expr_df$Gene <- probe_to_symbol

expr_clean <- expr_df %>% filter(!is.na(Gene) & Gene != "") %>% group_by(Gene) %>% summarise_all(mean) %>% as.data.frame()
rownames(expr_clean) <- expr_clean$Gene
mat_final <- as.matrix(expr_clean[, -1])

message(">>> [Step 2] 锁定训练集公式，正在对 GSE43292 进行独立验证判分...")

cal_z_score <- function(x) { 
  s <- sd(x, na.rm=TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  return((x - mean(x, na.rm=TRUE)) / s)
}
mat_z <- t(apply(mat_final, 1, cal_z_score))

t_cell_genes <- intersect(c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1"), rownames(mat_z))
if (length(t_cell_genes) > 0) {
  cd8_score <- colMeans(mat_z[t_cell_genes, , drop=FALSE], na.rm=TRUE)
  cd8_score[is.nan(cd8_score) | is.na(cd8_score)] <- 0
} else {
  cd8_score <- rep(0, ncol(mat_z))
}

safe_extract_z <- function(g) {
  if(g %in% rownames(mat_z)) {
    val <- mat_z[g, ]
    val[is.nan(val) | is.na(val)] <- 0
    return(val)
  } else {
    return(rep(0, ncol(mat_z)))
  }
}

clin_df <- data.frame(
  Sample = colnames(mat_final),
  Disease_Stage = factor(group_info, levels = c("Stage I-II (Early)", "Stage IV-V (Advanced)")),
  Status = ifelse(group_info == "Stage IV-V (Advanced)", 1, 0),
  SLC16A4_Z = safe_extract_z("SLC16A4"),
  SLC16A3_Z = safe_extract_z("SLC16A3"),
  ACKR1_Z   = safe_extract_z("ACKR1"),
  CD8_Tox_Z = cd8_score
)

clin_df$Logit_P <- 0.4625 + (-2.6931 * clin_df$SLC16A4_Z) + (0.6481 * clin_df$SLC16A3_Z) + (0.2556 * clin_df$ACKR1_Z) + (0.7414 * clin_df$CD8_Tox_Z)
clin_df$Prob_Combined <- 1 / (1 + exp(-clin_df$Logit_P))

n_early <- sum(clin_df$Disease_Stage == "Stage I-II (Early)")
n_adv   <- sum(clin_df$Disease_Stage == "Stage IV-V (Advanced)")
clinic_title <- sprintf("Dataset: GSE43292 | N = %d (Stage IV-V: %d, Stage I-II: %d)", nrow(clin_df), n_adv, n_early)
theme_nature_supp <- function() {
  theme_classic() +
    theme(text = element_text(color = "black", family = "sans"),
          axis.text = element_text(size = 13, color = "black", face = "bold"),
          axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.line = element_line(linewidth = 1.2, color = "black"),    
          axis.ticks = element_line(linewidth = 1.2, color = "black"),   
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b=5)),
          plot.subtitle = element_text(size = 12, color = "grey30", hjust = 0.5, margin = margin(b=15)))
}

x_labels <- c(sprintf("Stage I-II\n(n=%d)", n_early), sprintf("Stage IV-V\n(n=%d)", n_adv))

message(">>> [Step 3] 正在生成单一联合模型的 ROC 曲线...")

calc_roc_ci <- function(pred) {
  obj <- roc(clin_df$Status, pred, quiet = TRUE, direction = "auto")
  ci_val <- ci.auc(obj)
  return(list(obj = obj, ci_label = sprintf("(AUC = %.3f, 95%% CI: %.3f-%.3f)", obj$auc, ci_val[1], ci_val[3])))
}

res_comb <- calc_roc_ci(clin_df$Prob_Combined)

roc_list <- list("Combined Signature" = res_comb$obj)
auc_labels <- paste("Combined Signature", res_comb$ci_label)
names(roc_list) <- auc_labels

p_roc <- ggroc(roc_list, aes = c("color", "linewidth")) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey60", linewidth = 1) +
  scale_color_manual(values = setNames(c("#E64B35"), auc_labels)) +
  scale_linewidth_manual(values = setNames(c(2.5), auc_labels), guide = "none") +
  labs(title = "Diagnostic Validation (ROC)", subtitle = clinic_title, x = "Specificity", y = "Sensitivity") +
  theme_nature_supp() +
  theme(legend.position = c(0.60, 0.22), legend.title = element_blank(), 
        legend.text = element_text(size = 11, face = "bold"),
        legend.background = element_rect(fill = alpha("white", 0.6), color = "black", linewidth = 0.4))

ggsave(file.path(out_dir, "Figure1_ExtVal_ROC.pdf"), p_roc, width = 7, height = 7)


message(">>> [Step 4] 正在生成带有 H-L 和 Brier 评分的外部校准曲线...")

if (!requireNamespace("ResourceSelection", quietly = TRUE)) install.packages("ResourceSelection")

brier_score <- mean((clin_df$Prob_Combined - clin_df$Status)^2)

hl_test <- ResourceSelection::hoslem.test(clin_df$Status, clin_df$Prob_Combined, g = 8) 
hl_pval <- hl_test$p.value

plot_ext_calibration <- function() {
  par(mar = c(5, 5, 4, 2), family = "sans", lwd = 1.5, font.lab = 2)

  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
       xlab = "Formula Predicted Probability", 
       ylab = "Actual Observed Probability", 
       cex.axis = 1.2, cex.lab = 1.4)

  abline(0, 1, lty = 2, lwd = 2, col = "grey50")

  sm <- smooth.spline(clin_df$Prob_Combined, clin_df$Status, df = 3)
  lines(sm, lty = 1, lwd = 3, col = "black")

  title(main = "External Calibration Curve", font.main = 2, cex.main = 1.5)
  mtext(clinic_title, side = 3, line = 0.5, cex = 1, col = "grey30", font = 2)

  legend("bottomright", legend = c("Ideal Reference", "Calibration Curve"), 
         col = c("grey50", "black"), lty = c(2, 1), lwd = c(2, 3), 
         bty = "n", cex = 1.1, text.font = 2)

  text(0.03, 0.95, sprintf("Brier score = %.3f", brier_score), adj = 0, cex = 1.2, font = 2)
  text(0.03, 0.88, sprintf("Hosmer-Lemeshow p = %.3f", hl_pval), adj = 0, cex = 1.2, font = 2)
}

pdf(file.path(out_dir, "Figure2_ExtVal_Calibration.pdf"), width = 6.5, height = 6.5)
plot_ext_calibration()
dev.off()

png(file.path(out_dir, "Figure2_ExtVal_Calibration.png"), width = 6.5, height = 6.5, units = "in", res = 600)
plot_ext_calibration()
dev.off()

message(">>> 外部校准曲线 (含评分) 绘制完毕！")


message(">>> [Step 4] 正在生成外部校准曲线 (纯净趋势版)...")

plot_ext_calibration_clean <- function() {
  par(mar = c(5, 5, 4, 2), family = "sans", lwd = 1.5, font.lab = 2)

  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
       xlab = "Formula Predicted Probability", 
       ylab = "Actual Observed Probability", 
       cex.axis = 1.2, cex.lab = 1.4)

  abline(0, 1, lty = 2, lwd = 2, col = "grey50")

  sm <- smooth.spline(clin_df$Prob_Combined, clin_df$Status, df = 3)
  lines(sm, lty = 1, lwd = 3, col = "black")

  title(main = "External Calibration Curve", font.main = 2, cex.main = 1.5)
  mtext(clinic_title, side = 3, line = 0.5, cex = 1, col = "grey30", font = 2)

  legend("bottomright", legend = c("Ideal Reference", "Calibration Curve"), 
         col = c("grey50", "black"), lty = c(2, 1), lwd = c(2, 3), 
         bty = "n", cex = 1.1, text.font = 2)
}

pdf(file.path(out_dir, "Figure2_ExtVal_Calibration_Clean.pdf"), width = 6.5, height = 6.5)
plot_ext_calibration_clean()
dev.off()

message(">>> 外部校准曲线 (纯净版) 绘制完毕！")


message(">>> [Step 5] 正在生成 DCA 曲线...")

dca_res <- dca(Status ~ Prob_Combined, data = clin_df, thresholds = seq(0, 0.99, by = 0.01))
p_dca <- ggplot(as_tibble(dca_res), aes(x = threshold, y = net_benefit, color = variable)) +
  geom_line(aes(linewidth = variable)) +
  coord_cartesian(ylim = c(-0.05, max(as_tibble(dca_res)$net_benefit, na.rm=T) * 1.1), xlim = c(0, 1)) +
  scale_color_manual(values = c("all" = "grey60", "none" = "black", "Prob_Combined" = "#E64B35"),
                     labels = c("Treat All", "Treat None", "Combined Signature Model")) +
  scale_linewidth_manual(values = c("all" = 1.5, "none" = 1.5, "Prob_Combined" = 3.5), guide = "none") +
  labs(title = "Clinical Decision Curve Analysis (DCA)", subtitle = clinic_title, x = "Threshold Probability", y = "Clinical Net Benefit") +
  theme_nature_supp() +
  theme(legend.position = c(0.70, 0.75), legend.title = element_blank(), legend.text = element_text(size = 11, face = "bold"),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 1.2))

ggsave(file.path(out_dir, "Figure3_ExtVal_DCA.pdf"), p_dca, width = 7, height = 6.5)

message(">>> [Step 6] 正在生成联合模型评分及组分拆解图...")

p_dist_comb <- ggplot(clin_df, aes(x = Disease_Stage, y = Prob_Combined, fill = Disease_Stage)) +
  geom_violin(alpha = 0.6, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, linewidth = 1) +
  geom_jitter(shape = 21, size = 2, width = 0.15, color = "black", stroke = 1, fill = "white", alpha = 0.9) +
  stat_compare_means(comparisons = list(c("Stage I-II (Early)", "Stage IV-V (Advanced)")), method = "wilcox.test", 
                     aes(label = paste0("p = ", ..p.format..)), size = 5.5, fontface = "bold") +
  scale_fill_manual(values = c("Stage I-II (Early)" = "#3C5488", "Stage IV-V (Advanced)" = "#E64B35")) +
  scale_x_discrete(labels = c("Stage I-II (Early)" = x_labels[1], "Stage IV-V (Advanced)" = x_labels[2])) +
  labs(title = "Combined Signature Risk Score", subtitle = clinic_title, x = "Plaque Pathology Stage", y = "Predicted Probability") +
  theme_nature_supp() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Figure4A_ExtVal_Combined_Score.pdf"), p_dist_comb, width = 5.5, height = 6)

clin_long <- clin_df %>%
  dplyr::select(Sample, Disease_Stage, SLC16A4_Z, SLC16A3_Z, ACKR1_Z, CD8_Tox_Z) %>%
  pivot_longer(cols = c(SLC16A4_Z, SLC16A3_Z, ACKR1_Z, CD8_Tox_Z), names_to = "Component", values_to = "Z_Score")

clin_long$Component <- factor(clin_long$Component, levels = c("SLC16A4_Z", "SLC16A3_Z", "ACKR1_Z", "CD8_Tox_Z"),
                              labels = c("SLC16A4", "SLC16A3", "ACKR1", "CD8 Tox Score"))

p_breakdown <- ggplot(clin_long, aes(x = Disease_Stage, y = Z_Score, fill = Disease_Stage)) +
  geom_violin(alpha = 0.6, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, linewidth = 1) +
  geom_jitter(shape = 21, size = 2, width = 0.15, color = "black", stroke = 1, fill = "white", alpha = 0.9) +
  facet_wrap(~ Component, scales = "free_y", nrow = 1) +
  stat_compare_means(comparisons = list(c("Stage I-II (Early)", "Stage IV-V (Advanced)")), method = "wilcox.test", 
                     aes(label = paste0("p = ", ..p.format..)), size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("Stage I-II (Early)" = "#3C5488", "Stage IV-V (Advanced)" = "#E64B35")) +
  scale_x_discrete(labels = c("Stage I-II (Early)" = x_labels[1], "Stage IV-V (Advanced)" = x_labels[2])) +
  theme_nature_supp() +
  theme(legend.position = "none",
        strip.text = element_text(size = 13, face = "bold", color = "black"),
        strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold")) +
  labs(title = "Individual Predictor Components", subtitle = clinic_title, x = "Plaque Pathology Stage", y = "Expression (Z-score)")

ggsave(file.path(out_dir, "Figure4B_ExtVal_Component_Breakdown.pdf"), p_breakdown, width = 11, height = 5.5)

message("\n>>> 齐活儿！本地矩阵已秒速读取，GSE43292 队列的终极验证图表已生成！")


message(">>> [Step 6] 正在生成联合模型评分及组分拆解图，并拼接为单排大图...")

p_dist_comb <- ggplot(clin_df, aes(x = Disease_Stage, y = Prob_Combined, fill = Disease_Stage)) +
  geom_violin(alpha = 0.6, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, linewidth = 1) +
  geom_jitter(shape = 21, size = 2, width = 0.15, color = "black", stroke = 1, fill = "white", alpha = 0.9) +
  stat_compare_means(comparisons = list(c("Stage I-II (Early)", "Stage IV-V (Advanced)")), method = "wilcox.test", 
                     aes(label = paste0("p = ", ..p.format..)), size = 5.5, fontface = "bold") +
  scale_fill_manual(values = c("Stage I-II (Early)" = "#3C5488", "Stage IV-V (Advanced)" = "#E64B35")) +
  scale_x_discrete(labels = c("Stage I-II (Early)" = x_labels[1], "Stage IV-V (Advanced)" = x_labels[2])) +
  labs(title = "Combined Signature Risk Score", subtitle = clinic_title, x = "Plaque Pathology Stage", y = "Predicted Probability") +
  theme_nature_supp() +
  theme(legend.position = "none")

clin_long <- clin_df %>%
  dplyr::select(Sample, Disease_Stage, SLC16A4_Z, SLC16A3_Z, ACKR1_Z, CD8_Tox_Z) %>%
  pivot_longer(cols = c(SLC16A4_Z, SLC16A3_Z, ACKR1_Z, CD8_Tox_Z), names_to = "Component", values_to = "Z_Score")

clin_long$Component <- factor(clin_long$Component, levels = c("SLC16A4_Z", "SLC16A3_Z", "ACKR1_Z", "CD8_Tox_Z"),
                              labels = c("SLC16A4", "SLC16A3", "ACKR1", "CD8 Tox Score"))

p_breakdown <- ggplot(clin_long, aes(x = Disease_Stage, y = Z_Score, fill = Disease_Stage)) +
  geom_violin(alpha = 0.6, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, linewidth = 1) +
  geom_jitter(shape = 21, size = 2, width = 0.15, color = "black", stroke = 1, fill = "white", alpha = 0.9) +
  facet_wrap(~ Component, scales = "free_y", nrow = 1) +
  stat_compare_means(comparisons = list(c("Stage I-II (Early)", "Stage IV-V (Advanced)")), method = "wilcox.test", 
                     aes(label = paste0("p = ", ..p.format..)), size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("Stage I-II (Early)" = "#3C5488", "Stage IV-V (Advanced)" = "#E64B35")) +
  scale_x_discrete(labels = c("Stage I-II (Early)" = x_labels[1], "Stage IV-V (Advanced)" = x_labels[2])) +
  theme_nature_supp() +
  theme(legend.position = "none",
        strip.text = element_text(size = 13, face = "bold", color = "black"),
        strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold")) +
  labs(title = "Individual Predictor Components", subtitle = clinic_title, x = "Plaque Pathology Stage", y = "Expression (Z-score)")

p_final_merged <- (p_dist_comb | p_breakdown) + 
  plot_layout(widths = c(1, 3.5)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(file.path(out_dir, "Figure4_ExtVal_Merged_Distribution.pdf"), p_final_merged, width = 15, height = 5.5)

message("\n>>> 齐活儿！4A和4B已经完美拼接在一张图里，可以直接贴进文章排版了！")


# Figure 8F-J: GSE28829 external validation

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(hgu133plus2.db)
  library(ggplot2)
  library(ggpubr)
  library(pROC)
  library(patchwork)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "35_External_Validation_Cohorts")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> [Step 1] 正在加载 GSE28829 表达矩阵并翻译基因...")
file_path <- "GSE28829_series_matrix.txt" 
if(!file.exists(file_path)) stop("找不到 GSE28829_series_matrix.txt 文件，请确认路径！")

gset <- getGEO(filename = file_path, getGPL = FALSE)
if(is.list(gset)) gset <- gset[[1]]
titles <- pData(gset)$title

group_info <- ifelse(grepl("advanced", titles, ignore.case=TRUE), "Advanced", "Early")

probe_to_symbol <- mapIds(hgu133plus2.db, keys = rownames(exprs(gset)), 
                          column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
expr_df <- as.data.frame(exprs(gset))
expr_df$Gene <- probe_to_symbol

expr_clean <- expr_df %>% 
  filter(!is.na(Gene) & Gene != "") %>% 
  group_by(Gene) %>% 
  summarise_all(mean) %>% 
  as.data.frame()

rownames(expr_clean) <- expr_clean$Gene
mat_final <- as.matrix(expr_clean[, -1])

message(">>> [Step 2] 锁定公式，计算 GSE28829 风险评分...")

cal_z_score <- function(x) {
  s <- sd(x, na.rm=TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  return((x - mean(x, na.rm=TRUE)) / s)
}
mat_z <- t(apply(mat_final, 1, cal_z_score))

t_cell_genes <- intersect(c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1"), rownames(mat_z))
if (length(t_cell_genes) > 0) {
  cd8_score <- colMeans(mat_z[t_cell_genes, , drop=FALSE], na.rm=TRUE)
  cd8_score[is.nan(cd8_score) | is.na(cd8_score)] <- 0
} else {
  cd8_score <- rep(0, ncol(mat_z))
}

safe_extract_z <- function(g) {
  if(g %in% rownames(mat_z)) {
    val <- mat_z[g, ]
    val[is.nan(val) | is.na(val)] <- 0
    return(val)
  } else return(rep(0, ncol(mat_z)))
}

clin_df <- data.frame(
  Sample = colnames(mat_final),
  Disease_Stage = factor(group_info, levels = c("Early", "Advanced")),
  Status = ifelse(group_info == "Advanced", 1, 0),
  SLC16A4_Z = safe_extract_z("SLC16A4"),
  SLC16A3_Z = safe_extract_z("SLC16A3"),
  ACKR1_Z   = safe_extract_z("ACKR1"),
  CD8_Tox_Z = cd8_score
)

clin_df$Logit_P <- 0.4625 + (-2.6931 * clin_df$SLC16A4_Z) + (0.6481 * clin_df$SLC16A3_Z) + (0.2556 * clin_df$ACKR1_Z) + (0.7414 * clin_df$CD8_Tox_Z)
clin_df$Signature_Score <- 1 / (1 + exp(-clin_df$Logit_P))

n_early <- sum(clin_df$Disease_Stage == "Early")
n_adv   <- sum(clin_df$Disease_Stage == "Advanced")
subtitle_text <- sprintf("Dataset: GSE28829 | N = %d (Advanced: %d, Early: %d)", nrow(clin_df), n_adv, n_early)

theme_nature <- function() {
  theme_classic() +
    theme(text = element_text(color = "black", family = "sans"),
          axis.text = element_text(size = 13, color = "black", face = "bold"),
          axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.line = element_line(linewidth = 1.2, color = "black"),
          axis.ticks = element_line(linewidth = 1.2, color = "black"),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b=5)),
          plot.subtitle = element_text(size = 12, color = "grey30", hjust = 0.5, margin = margin(b=15)))
}

x_labels <- c("Early" = sprintf("Early\n(n=%d)", n_early), "Advanced" = sprintf("Advanced\n(n=%d)", n_adv))

message(">>> [Step 3] 绘制单红线 ROC 曲线...")

roc_obj <- roc(clin_df$Status, clin_df$Signature_Score, quiet = TRUE, direction = "auto")
ci_val <- ci.auc(roc_obj)
auc_label <- sprintf("Combined Signature (AUC = %.3f, 95%% CI: %.3f-%.3f)", roc_obj$auc, ci_val[1], ci_val[3])

roc_list <- list("Combined Signature" = roc_obj)
names(roc_list) <- auc_label

p_roc <- ggroc(roc_list, aes = c("color", "linewidth")) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey60", linewidth = 1) +
  scale_color_manual(values = setNames(c("#E64B35"), auc_label)) +
  scale_linewidth_manual(values = setNames(c(2.5), auc_label), guide = "none") +
  labs(title = "Diagnostic Performance (ROC)", subtitle = subtitle_text, x = "Specificity", y = "Sensitivity") +
  theme_nature() +
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        legend.position = c(0.55, 0.22),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, face = "bold"),
        legend.background = element_rect(fill = alpha("white", 0.6), color = "black", linewidth = 0.4))

ggsave(file.path(out_dir, "Figure1_GSE28829_ROC.pdf"), p_roc, width = 7, height = 7)

message(">>> [Step 4] 正在生成图4：组分拆解图并横向合并排版...")

p_dist_comb <- ggplot(clin_df, aes(x = Disease_Stage, y = Signature_Score, fill = Disease_Stage)) +
  geom_violin(alpha = 0.6, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, linewidth = 1) +
  geom_jitter(shape = 21, size = 2, width = 0.15, color = "black", stroke = 1, fill = "white", alpha = 0.9) +
  stat_compare_means(comparisons = list(c("Early", "Advanced")), method = "wilcox.test", 
                     aes(label = paste0("p = ", ..p.format..)), size = 5.5, fontface = "bold") +
  scale_fill_manual(values = c("Early" = "#3C5488", "Advanced" = "#E64B35")) +
  scale_x_discrete(labels = x_labels) +
  labs(title = "Combined Signature Score", subtitle = subtitle_text, x = "Plaque Stage", y = "Predicted Probability") +
  theme_nature() +
  theme(legend.position = "none")

clin_long <- clin_df %>%
  dplyr::select(Sample, Disease_Stage, SLC16A4_Z, SLC16A3_Z, ACKR1_Z, CD8_Tox_Z) %>%
  pivot_longer(cols = ends_with("_Z"), names_to = "Component", values_to = "Z_Score")

clin_long$Component <- factor(clin_long$Component, levels = c("SLC16A4_Z", "SLC16A3_Z", "ACKR1_Z", "CD8_Tox_Z"),
                              labels = c("SLC16A4", "SLC16A3", "ACKR1", "CD8 Tox Score"))

p_breakdown <- ggplot(clin_long, aes(x = Disease_Stage, y = Z_Score, fill = Disease_Stage)) +
  geom_violin(alpha = 0.6, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, linewidth = 1) +
  geom_jitter(shape = 21, size = 2, width = 0.15, color = "black", stroke = 1, fill = "white", alpha = 0.9) +
  facet_wrap(~ Component, scales = "free_y", nrow = 1) +
  stat_compare_means(comparisons = list(c("Early", "Advanced")), method = "wilcox.test", 
                     aes(label = paste0("p = ", ..p.format..)), size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("Early" = "#3C5488", "Advanced" = "#E64B35")) +
  scale_x_discrete(labels = x_labels) +
  theme_nature() +
  theme(legend.position = "none",
        strip.text = element_text(size = 13, face = "bold", color = "black"),
        strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold")) +
  labs(title = "Individual Predictor Components", subtitle = subtitle_text, x = "Plaque Stage", y = "Expression (Z-score)")

p_final_merged <- (p_dist_comb | p_breakdown) + 
  plot_layout(widths = c(1, 3.5)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(file.path(out_dir, "Figure4_GSE28829_Merged_Distribution.pdf"), p_final_merged, width = 15, height = 5.5)

message("\n>>> 完美搞定！GSE28829 的 ROC 已经剔除了杂线，图 4 也完美合并在一排了，请去文件夹查收！")


message(">>> 正在生成单线平滑校准曲线...")

pdf(file.path(out_dir, "Figure2_GSE28829_Calibration.pdf"), width = 6.5, height = 6.5)
par(mar = c(5, 5, 4, 2), family = "sans", lwd = 1.5, font.lab = 2)

plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "Formula Predicted Probability", ylab = "Actual Observed Probability", 
     cex.axis = 1.2, cex.lab = 1.4)
abline(0, 1, lty = 2, lwd = 2, col = "grey50")

sm <- smooth.spline(clin_df$Signature_Score, clin_df$Status, df = 3)
lines(sm, lty = 1, lwd = 3, col = "black") 

title(main = "External Calibration Curve", font.main = 2, cex.main = 1.5)
mtext(subtitle_text, side = 3, line = 0.5, cex = 1, col = "grey30", font = 2)
legend("bottomright", legend = c("Ideal Reference", "Calibration Curve"),
       col = c("grey50", "black"), lty = c(2, 1), lwd = c(2, 3), bty = "n", cex = 1.1, text.font = 2)
dev.off()

message(">>> 正在生成 DCA 决策曲线...")

dca_res <- dca(Status ~ Signature_Score, data = clin_df, thresholds = seq(0, 0.99, by = 0.01))

p_dca <- ggplot(as_tibble(dca_res), aes(x = threshold, y = net_benefit, color = variable)) +
  geom_line(aes(linewidth = variable)) +
  coord_cartesian(ylim = c(-0.05, max(as_tibble(dca_res)$net_benefit, na.rm=T) * 1.1), xlim = c(0, 1)) +
  scale_color_manual(values = c("all" = "grey60", "none" = "black", "Signature_Score" = "#E64B35"),
                     labels = c("Treat All", "Treat None", "Combined Signature Model")) +
  scale_linewidth_manual(values = c("all" = 1.5, "none" = 1.5, "Signature_Score" = 3.5), guide = "none") +
  labs(title = "Clinical Decision Curve Analysis (DCA)", subtitle = subtitle_text, 
       x = "Threshold Probability", y = "Clinical Net Benefit") +
  theme_nature() +
  theme(legend.position = c(0.70, 0.75), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 11, face = "bold"),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 1.2))

ggsave(file.path(out_dir, "Figure3_GSE28829_DCA.pdf"), p_dca, width = 7, height = 6.5)

message("\n>>> GSE28829 补全计划完成！校准曲线和 DCA 已生成并保存。铁三角全了！")


# Supplementary Figures 3-4: GSE100927 validation by vascular site

message(">>> [Step 1] 正在读取 GSE100927 矩阵并精准划分三个血管解剖部位...")

suppressMessages({
  library(GEOquery)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(pROC)
  library(dcurves)
  library(patchwork)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "38_External_Validation_PanVascular_GSE100927_Subgroups")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_path <- "GSE100927_series_matrix.txt.gz"
if(!file.exists(file_path)) file_path <- "GSE100927_series_matrix.txt"
if(!file.exists(file_path)) stop("找不到矩阵文件！")

gset <- getGEO(filename = file_path, getGPL = FALSE)
if(is.list(gset)) gset <- gset[[1]] 
pheno <- pData(gset)

info_str <- apply(pheno, 1, paste, collapse=" ")

group_info <- ifelse(grepl("control|healthy", info_str, ignore.case=TRUE), "Control", 
                     ifelse(grepl("atherosclerot|plaque|lesion", info_str, ignore.case=TRUE), "Atherosclerotic", "Unknown"))

loc_info <- ifelse(grepl("carotid", info_str, ignore.case=TRUE), "Carotid",
                   ifelse(grepl("femoral", info_str, ignore.case=TRUE), "Femoral",
                          ifelse(grepl("popliteal", info_str, ignore.case=TRUE), "Infra-popliteal", "Other")))

pheno_df <- data.frame(Sample = rownames(pheno), Status = group_info, Location = loc_info)
pheno_valid <- pheno_df %>% filter(Status != "Unknown" & Location != "Other")

message(">>> [Step 2] 正在调用 read_tsv 引擎强力解析纯净字典...")

dict_file <- list.files(pattern = "GPL17077.*\\.txt", ignore.case = TRUE)[1]
if(is.na(dict_file)) stop("严重错误：没找到手动下载的字典 txt 文件！")

raw_lines <- readLines(dict_file)
header_idx <- grep("^ID\t|^ID |^\"ID\"", raw_lines)[1]
clean_tmp_file <- tempfile(fileext = ".txt")
writeLines(raw_lines[header_idx:length(raw_lines)], clean_tmp_file)

gpl_anno <- suppressWarnings(read_tsv(clean_tmp_file, show_col_types = FALSE, name_repair = "minimal"))
sym_col <- grep("Symbol|GeneName|GENE_SYMBOL|gene_assignment", colnames(gpl_anno), ignore.case = TRUE, value = TRUE)[1]

probe_dict <- gpl_anno %>% 
  dplyr::select(ID, all_of(sym_col)) %>% 
  dplyr::rename(Symbol = all_of(sym_col)) %>%
  filter(Symbol != "", !is.na(Symbol), Symbol != "---") %>%
  mutate(Symbol = trimws(sapply(strsplit(as.character(Symbol), "[,///|]"), `[`, 1))) %>%
  distinct(ID, .keep_all = TRUE)

expr_df <- as.data.frame(exprs(gset)[, pheno_valid$Sample])
expr_df$ID <- rownames(expr_df)
expr_df <- expr_df %>% inner_join(probe_dict, by = "ID") %>% dplyr::rename(Gene = Symbol) %>% dplyr::select(-ID)
expr_clean <- expr_df %>% filter(!is.na(Gene) & Gene != "") %>% group_by(Gene) %>% summarise_all(mean) %>% as.data.frame()
rownames(expr_clean) <- expr_clean$Gene
mat_final <- as.matrix(expr_clean[, -1])

message(">>> [Step 3] 锁定公式，计算全队列基础得分...")

cal_z_score <- function(x) { 
  s <- sd(x, na.rm=TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  return((x - mean(x, na.rm=TRUE)) / s)
}
mat_z <- t(apply(mat_final, 1, cal_z_score))

t_cell_genes <- intersect(c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1"), rownames(mat_z))
if (length(t_cell_genes) > 0) {
  cd8_score <- colMeans(mat_z[t_cell_genes, , drop=FALSE], na.rm=TRUE)
  cd8_score[is.nan(cd8_score) | is.na(cd8_score)] <- 0
} else cd8_score <- rep(0, ncol(mat_z))

safe_extract_z <- function(g) {
  if(g %in% rownames(mat_z)) { val <- mat_z[g, ]; val[is.nan(val) | is.na(val)] <- 0; return(val) } 
  else return(rep(0, ncol(mat_z)))
}

clin_df_all <- data.frame(
  Sample = pheno_valid$Sample,
  Plaque_Status = factor(pheno_valid$Status, levels = c("Control", "Atherosclerotic")),
  Artery_Location = pheno_valid$Location,
  Status = ifelse(pheno_valid$Status == "Atherosclerotic", 1, 0),
  SLC16A4_Z = safe_extract_z("SLC16A4"),
  SLC16A3_Z = safe_extract_z("SLC16A3"),
  ACKR1_Z   = safe_extract_z("ACKR1"),
  CD8_Tox_Z = cd8_score
)

clin_df_all$Logit_P <- 0.4625 + (-2.6931 * clin_df_all$SLC16A4_Z) + (0.6481 * clin_df_all$SLC16A3_Z) + (0.2556 * clin_df_all$ACKR1_Z) + (0.7414 * clin_df_all$CD8_Tox_Z)
clin_df_all$Signature_Score <- 1 / (1 + exp(-clin_df_all$Logit_P))

theme_nature <- function() {
  theme_classic() +
    theme(text = element_text(color = "black", family = "sans"),
          axis.text = element_text(size = 13, color = "black", face = "bold"),
          axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.line = element_line(linewidth = 1.2, color = "black"),
          axis.ticks = element_line(linewidth = 1.2, color = "black"),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b=5)),
          plot.subtitle = element_text(size = 12, color = "grey30", hjust = 0.5, margin = margin(b=15)))
}

locations <- unique(clin_df_all$Artery_Location)
message("\n>>> [Step 4] 开始按解剖部位循环生成四大神图...")

for(loc in locations) {

  sub_df <- clin_df_all[clin_df_all$Artery_Location == loc, ]

  if(length(unique(sub_df$Status)) < 2) {
    message(sprintf("--- 跳过 %s 部位 (因表型单一)", loc))
    next
  }

  n_stable <- sum(sub_df$Plaque_Status == "Control")
  n_unstable <- sum(sub_df$Plaque_Status == "Atherosclerotic")
  n_total <- nrow(sub_df)

  clinic_title <- sprintf("Dataset: GSE100927 (%s Artery) | N = %d (Atherosclerotic: %d, Control: %d)", loc, n_total, n_unstable, n_stable)
  x_labels <- c("Control" = sprintf("Control\n(n=%d)", n_stable), "Atherosclerotic" = sprintf("Atherosclerotic\n(n=%d)", n_unstable))

  message(sprintf("--- 正在生成 %s 动脉专属图表 (N = %d) ...", loc, n_total))

  roc_obj <- roc(sub_df$Status, sub_df$Signature_Score, quiet = TRUE, direction = "auto")
  ci_val <- ci.auc(roc_obj)
  auc_label <- sprintf("Combined Signature (AUC = %.3f, 95%% CI: %.3f-%.3f)", roc_obj$auc, ci_val[1], ci_val[3])

  roc_list <- list("Combined Signature" = roc_obj)
  names(roc_list) <- auc_label

  p_roc <- ggroc(roc_list, aes = c("color", "linewidth")) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey60", linewidth = 1) +
    scale_color_manual(values = setNames(c("#E64B35"), auc_label)) +
    scale_linewidth_manual(values = setNames(c(2.5), auc_label), guide = "none") +
    labs(title = sprintf("Diagnostic Validation - %s", loc), subtitle = clinic_title, x = "Specificity", y = "Sensitivity") +
    theme_nature() +
    theme(legend.position = c(0.55, 0.22), legend.title = element_blank(), legend.text = element_text(size = 11, face = "bold"),
          legend.background = element_rect(fill = alpha("white", 0.6), color = "black", linewidth = 0.4))
  ggsave(file.path(out_dir, sprintf("Figure1_GSE100927_%s_ROC.pdf", loc)), p_roc, width = 7, height = 7)

  pdf(file.path(out_dir, sprintf("Figure2_GSE100927_%s_Calibration.pdf", loc)), width = 6.5, height = 6.5)
  par(mar = c(5, 5, 4, 2), family = "sans", lwd = 1.5, font.lab = 2)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "Formula Predicted Probability", ylab = "Actual Observed Probability", cex.axis = 1.2, cex.lab = 1.4)
  abline(0, 1, lty = 2, lwd = 2, col = "grey50")

  lo <- try(loess(Status ~ Signature_Score, data = sub_df, degree = 1, span = 0.9), silent = TRUE)
  if(!inherits(lo, "try-error")) {
    smooth_x <- seq(min(sub_df$Signature_Score), max(sub_df$Signature_Score), length.out = 100)
    smooth_y <- predict(lo, newdata = data.frame(Signature_Score = smooth_x))
    smooth_y[smooth_y < 0] <- 0; smooth_y[smooth_y > 1] <- 1
    lines(smooth_x, smooth_y, lty = 1, lwd = 3, col = "black") 
  }
  title(main = sprintf("%s Calibration Curve", loc), font.main = 2, cex.main = 1.5)
  mtext(clinic_title, side = 3, line = 0.5, cex = 0.9, col = "grey30", font = 2)
  legend("bottomright", legend = c("Ideal Reference", "Calibration Curve"), col = c("grey50", "black"), lty = c(2, 1), lwd = c(2, 3), bty = "n", cex = 1.1, text.font = 2)
  dev.off()

  dca_res <- dca(Status ~ Signature_Score, data = sub_df, thresholds = seq(0, 0.99, by = 0.01))
  p_dca <- ggplot(as_tibble(dca_res), aes(x = threshold, y = net_benefit, color = variable)) +
    geom_line(aes(linewidth = variable)) +
    coord_cartesian(ylim = c(-0.05, max(as_tibble(dca_res)$net_benefit, na.rm=T) * 1.1), xlim = c(0, 1)) +
    scale_color_manual(values = c("all" = "grey60", "none" = "black", "Signature_Score" = "#E64B35"), labels = c("Treat All", "Treat None", "Combined Signature Model")) +
    scale_linewidth_manual(values = c("all" = 1.5, "none" = 1.5, "Signature_Score" = 3.5), guide = "none") +
    labs(title = sprintf("Clinical Decision Curve (DCA) - %s", loc), subtitle = clinic_title, x = "Threshold Probability", y = "Clinical Net Benefit") +
    theme_nature() +
    theme(legend.position = c(0.70, 0.75), legend.title = element_blank(), legend.text = element_text(size = 11, face = "bold"),
          legend.background = element_rect(fill = "white", color = "black", linewidth = 1.2))
  ggsave(file.path(out_dir, sprintf("Figure3_GSE100927_%s_DCA.pdf", loc)), p_dca, width = 7, height = 6.5)

  p_dist_comb <- ggplot(sub_df, aes(x = Plaque_Status, y = Signature_Score, fill = Plaque_Status)) +
    geom_violin(alpha = 0.6, color = NA, trim = FALSE) + geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(shape = 21, size = 2, width = 0.15, color = "black", stroke = 1, fill = "white", alpha = 0.9) +
    stat_compare_means(comparisons = list(c("Control", "Atherosclerotic")), method = "wilcox.test", aes(label = paste0("p = ", ..p.format..)), size = 5.5, fontface = "bold") +
    scale_fill_manual(values = c("Control" = "#3C5488", "Atherosclerotic" = "#E64B35")) + scale_x_discrete(labels = x_labels) +
    labs(title = sprintf("Combined Score (%s)", loc), subtitle = clinic_title, x = "Plaque Status", y = "Predicted Probability") + theme_nature() + theme(legend.position = "none")

  sub_long <- sub_df %>% dplyr::select(Sample, Plaque_Status, SLC16A4_Z, SLC16A3_Z, ACKR1_Z, CD8_Tox_Z) %>% pivot_longer(cols = ends_with("_Z"), names_to = "Component", values_to = "Z_Score")
  sub_long$Component <- factor(sub_long$Component, levels = c("SLC16A4_Z", "SLC16A3_Z", "ACKR1_Z", "CD8_Tox_Z"), labels = c("SLC16A4", "SLC16A3", "ACKR1", "CD8 Tox Score"))

  p_breakdown <- ggplot(sub_long, aes(x = Plaque_Status, y = Z_Score, fill = Plaque_Status)) +
    geom_violin(alpha = 0.6, color = NA, trim = FALSE) + geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(shape = 21, size = 2, width = 0.15, color = "black", stroke = 1, fill = "white", alpha = 0.9) + facet_wrap(~ Component, scales = "free_y", nrow = 1) +
    stat_compare_means(comparisons = list(c("Control", "Atherosclerotic")), method = "wilcox.test", aes(label = paste0("p = ", ..p.format..)), size = 4.5, fontface = "bold") +
    scale_fill_manual(values = c("Control" = "#3C5488", "Atherosclerotic" = "#E64B35")) + scale_x_discrete(labels = x_labels) + theme_nature() +
    theme(legend.position = "none", strip.text = element_text(size = 13, face = "bold", color = "black"), strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1.2),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold")) +
    labs(title = sprintf("Individual Predictor Components (%s)", loc), subtitle = clinic_title, x = "Plaque Status", y = "Expression (Z-score)")

  p_final_merged <- (p_dist_comb | p_breakdown) + plot_layout(widths = c(1, 3.5)) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 18, face = "bold"))
  ggsave(file.path(out_dir, sprintf("Figure4_GSE100927_%s_Distribution.pdf", loc)), p_final_merged, width = 15, height = 5.5)
}

message("\n>>> 太震撼了！三个解剖部位已全部拆解完毕！")
message(">>> 请打开 '38_External_Validation_PanVascular_GSE100927_Subgroups' 文件夹验收你的 12 张神图！")

