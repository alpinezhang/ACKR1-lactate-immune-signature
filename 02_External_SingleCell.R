# Supplementary Figure 1: ACKR1 external single-cell analysis (English; current revision order)
# Corresponds to current Supplementary Figure 1 (GSE253903 external single-cell panels A-J).


suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(harmony)    
  library(DoubletFinder) 
})

work_dir <- "PATH_TO_DATA/GSE253903"
out_dir  <- file.path(work_dir, "01_QC_and_Clustering")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 当前工作目录已设置为: ", getwd())
message(">>> 输出图表将保存在: ", out_dir)

group_info <- data.frame(
  GSM = c("GSM8029963", "GSM8029967", "GSM8029968", "GSM8029969", "GSM8029971", "GSM8029972",
          "GSM8029962", "GSM8029964", "GSM8029965", "GSM8029966", "GSM8029970", "GSM8029973"),
  Condition = c(rep("Asymptomatic", 6), rep("Symptomatic", 6)),
  stringsAsFactors = FALSE
)

seurat_list <- list()
message(">>> 正在批量读取 GSE253903 单细胞数据...")

for (i in 1:nrow(group_info)) {
  gsm_id <- group_info$GSM[i]
  folder_path <- file.path(work_dir, gsm_id)

  if (dir.exists(folder_path)) {
    counts <- tryCatch({
      Read10X(data.dir = folder_path)
    }, error = function(e) {
      tryCatch({
        Read10X(data.dir = folder_path, gene.column = 1)
      }, error = function(e2) {
        return(NULL)
      })
    })

    if (is.null(counts) || length(counts) == 0) {
      message("   ! 错误: 无法解析文件夹 ", gsm_id, " 的表达矩阵，已跳过。请检查内部文件命名！")
      next
    }

    if (is.list(counts) && "Gene Expression" %in% names(counts)) {
      counts <- counts$`Gene Expression`
    }

    sobj <- CreateSeuratObject(counts, project = gsm_id)
    sobj$Orig_Sample <- gsm_id
    sobj$Condition   <- group_info$Condition[i]

    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
    sobj[["percent.hb"]] <- PercentageFeatureSet(sobj, pattern = "^HB[^P]")

    seurat_list[[gsm_id]] <- sobj
    message("   - 成功读取样本: ", gsm_id)
  } else {
    message("   ! 警告: 找不到文件夹 ", gsm_id, "，已跳过。")
  }
}

sc_raw <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "GSE253903")

message(">>> 正在生成全局质控图 (Pre-QC)...")

nature_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

meta <- sc_raw@meta.data
features_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb")
features_to_plot <- intersect(features_to_plot, colnames(meta))

plot_data <- data.frame()
for (feat in features_to_plot) {
  temp <- data.frame(
    Sample = meta$Orig_Sample,
    Condition = meta$Condition,
    Value = meta[[feat]],
    Feature = feat
  )
  plot_data <- rbind(plot_data, temp)
}
plot_data$Feature <- factor(plot_data$Feature, levels = features_to_plot)

p_qc <- ggplot(plot_data, aes(x = Sample, y = Value, fill = Condition)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.5, alpha = 0.8) + 
  facet_wrap(~ Feature, scales = "free_y", ncol = length(features_to_plot)) +
  scale_fill_npg() + 
  nature_theme +
  labs(x = "Samples", y = "Value", title = "Pre-QC Metrics across Samples (GSE253903)")

ggsave(file.path(out_dir, "00_Pre_QC_Violin.pdf"), plot = p_qc, width = 14, height = 5)
ggsave(file.path(out_dir, "00_Pre_QC_Violin.png"), plot = p_qc, width = 14, height = 5, dpi = 300)

message(">>> 提示: 过滤前的质控图已保存至 ", file.path(out_dir, "00_Pre_QC_Violin.pdf"))

min_features <- 500
max_features <- 6000
max_mt       <- 15
max_hb       <- 5

message(sprintf(">>> 正在应用过滤条件: nFeature %d-%d, 线粒体 < %d%%, 红细胞 < %d%%", 
                min_features, max_features, max_mt, max_hb))

for(i in seq_along(seurat_list)){
  if("percent.hb" %in% colnames(seurat_list[[i]]@meta.data)) {
    seurat_list[[i]] <- subset(seurat_list[[i]], 
                               subset = nFeature_RNA > min_features & 
                                 nFeature_RNA < max_features & 
                                 percent.mt < max_mt &
                                 percent.hb < max_hb)
  } else {
    seurat_list[[i]] <- subset(seurat_list[[i]], 
                               subset = nFeature_RNA > min_features & 
                                 nFeature_RNA < max_features & 
                                 percent.mt < max_mt)
  }
}

message(">>> 开始执行 DoubletFinder 去双细胞...")

for (i in seq_along(seurat_list)) {
  sample_name <- names(seurat_list)[i]
  sobj <- seurat_list[[i]]

  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(sobj, verbose = FALSE)
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = 1:15, verbose = FALSE)

  doublet_rate <- (ncol(sobj) / 1000) * 0.008 
  nExp_poi <- round(doublet_rate * ncol(sobj))

  sobj <- doubletFinder_v3(sobj, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  df_col <- grep("DF.classifications", colnames(sobj@meta.data), value = TRUE)
  sobj$Doublet_Status <- sobj@meta.data[[df_col]]

  sobj <- subset(sobj, subset = Doublet_Status == "Singlet")
  seurat_list[[i]] <- sobj
  message(sprintf("   - 样本 %s 去双细胞完成，保留 %d 细胞", sample_name, ncol(sobj)))
}

message(">>> 正在合并数据，并运行 Harmony 去除 GSM 之间的批次效应...")

sc_filtered <- merge(seurat_list[[1]], y = seurat_list[-1], 
                     add.cell.ids = names(seurat_list), project = "GSE253903")
rm(seurat_list, sc_raw); gc() 

if (inherits(sc_filtered[["RNA"]], "Assay5")) {
  message("   - 检测到 Seurat V5，正在合并图层 (JoinLayers)...")
  sc_filtered <- JoinLayers(sc_filtered)
}

sc_filtered <- NormalizeData(sc_filtered, verbose = FALSE)
sc_filtered <- FindVariableFeatures(sc_filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sc_filtered <- ScaleData(sc_filtered, verbose = FALSE)
sc_filtered <- RunPCA(sc_filtered, npcs = 30, verbose = FALSE)

sc_filtered <- RunHarmony(sc_filtered, group.by.vars = "Orig_Sample", plot_convergence = FALSE, verbose = FALSE)

message(">>> 正在进行 UMAP 降维与分群...")
sc_filtered <- RunUMAP(sc_filtered, reduction = "harmony", dims = 1:20, verbose = FALSE)
sc_filtered <- FindNeighbors(sc_filtered, reduction = "harmony", dims = 1:20, verbose = FALSE)
sc_filtered <- FindClusters(sc_filtered, resolution = 0.6, verbose = FALSE)

message(">>> 正在生成 UMAP 聚类图与分组对比图...")

num_clusters <- length(levels(sc_filtered$seurat_clusters))
cluster_colors <- colorRampPalette(pal_npg("nrc")(10))(num_clusters)

p_cluster <- DimPlot(sc_filtered, reduction = "umap", group.by = "seurat_clusters", 
                     label = TRUE, label.size = 5, pt.size = 0.1) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("UMAP by Clusters") +
  nature_theme + 
  theme(legend.position = "right")

p_group <- DimPlot(sc_filtered, reduction = "umap", group.by = "Condition", 
                   pt.size = 0.1) +
  scale_color_npg() + 
  ggtitle("UMAP by Clinical Condition") +
  nature_theme + 
  theme(legend.position = "right")

p_combined <- p_cluster + p_group + plot_layout(ncol = 2)

ggsave(file.path(out_dir, "01_UMAP_Combined.pdf"), p_combined, width = 12, height = 5)
ggsave(file.path(out_dir, "01_UMAP_Combined.png"), p_combined, width = 12, height = 5, dpi = 300)

message(">>> 正在保存 Metadata 表格与 RDS...")
write.csv(sc_filtered@meta.data, file.path(out_dir, "02_Filtered_Metadata.csv"), row.names = TRUE)

saveRDS(sc_filtered, file.path(work_dir, "GSE253903_Integrated_Unannotated.rds"))
message(">>> GSE253903 单细胞处理流程全部完成！")


# Supplementary Figure 1A: External-cohort UMAP
message(">>> 正在进行高级 UMAP 图形渲染 (带图例版本)...")

umap_data <- as.data.frame(Embeddings(sc_filtered, "umap"))
umap_data$Cluster <- Idents(sc_filtered)

cluster_centers <- umap_data %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = 'drop')

min_x <- min(umap_data$UMAP_1)
min_y <- min(umap_data$UMAP_2)
range_x <- diff(range(umap_data$UMAP_1))
range_y <- diff(range(umap_data$UMAP_2))

len_x <- range_x * 0.15
len_y <- range_y * 0.15

strong_colors <- c(
  "#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", 
  "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4", 
  "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", 
  "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", "#3D3D3D",
  "#1F77B4"
)
num_clusters <- length(levels(umap_data$Cluster))
custom_colors <- strong_colors[1:num_clusters]

p_umap_nature <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = 0.4, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = custom_colors) +

  geom_text_repel(data = cluster_centers, aes(label = Cluster),
                  size = 5, fontface = "bold", color = "black",
                  bg.color = "white", bg.r = 0.15, 
                  segment.color = "black", segment.size = 0.5, 
                  max.overlaps = Inf, show.legend = FALSE) +

  theme_void() +
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15))) + 

  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +

  annotate("segment", x = min_x, xend = min_x + len_x, y = min_y, yend = min_y,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed"), linewidth = 0.8) +
  annotate("segment", x = min_x, xend = min_x, y = min_y, yend = min_y + len_y,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed"), linewidth = 0.8) +

  annotate("text", x = min_x + len_x/2, y = min_y - range_y*0.04, 
           label = "UMAP_1", fontface = "italic", size = 4) +
  annotate("text", x = min_x - range_x*0.04, y = min_y + len_y/2, 
           label = "UMAP_2", fontface = "italic", size = 4, angle = 90) +

  ggtitle("UMAP of GSE253903 Clusters (Unannotated)")

ggsave(file.path(out_dir, "02_UMAP_Clusters_NatureStyle.pdf"), p_umap_nature, width = 9, height = 8)

# Supplementary Figure 1B: Cell-type annotation and marker dot plot

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

work_dir <- "PATH_TO_DATA/GSE253903"
out_dir  <- file.path(work_dir, "02_Annotation_and_Visualization")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 正在加载单细胞数据...")
sc_filtered <- readRDS("GSE253903_Integrated_Unannotated.rds")

message(">>> 正在为 Cluster 添加细胞类型注释...")
Idents(sc_filtered) <- "seurat_clusters" 

cluster2celltype <- c(
  "0"  = "0 T cells",
  "1"  = "1 Monocytes",
  "2"  = "2 CD8+T cells",
  "3"  = "3 SMCs",
  "4"  = "4 SMCs",
  "5"  = "5 Macrophages", 
  "6"  = "6 Monocytes",
  "7"  = "7 ACKR1+ ECs",
  "8"  = "8 cDC2",
  "9"  = "9 Proliferating",
  "10" = "10 ACKR1- ECs",
  "11" = "11 Plasma cells",
  "12" = "12 Mast cells",
  "13" = "13 SMCs",
  "14" = "14 SMCs",
  "15" = "15 B cells",
  "16" = "16 pDCs"
)

sc_filtered <- RenameIdents(sc_filtered, cluster2celltype)
sc_filtered$cell_type_annotated <- Idents(sc_filtered)

grouped_cluster_order <- c(
  "0 T cells", "2 CD8+T cells",
  "1 Monocytes", "6 Monocytes", 
  "5 Macrophages", 
  "8 cDC2", "16 pDCs",
  "7 ACKR1+ ECs", "10 ACKR1- ECs",
  "3 SMCs", "4 SMCs", "13 SMCs","14 SMCs",
  "15 B cells", "11 Plasma cells",
  "12 Mast cells", 
  "9 Proliferating"
)

Idents(sc_filtered) <- factor(Idents(sc_filtered), levels = rev(grouped_cluster_order))

marker_list <- list(
  T_cells       = c("CD3D", "CD3E", "CD8A", "CD8B"),
  NKT_cells     = c("CD3D", "CD3E", "CD8A", "CD8B", "GNLY", "NKG7", "PRF1"),
  Monocytes     = c("FCGR3A", "FCN1"),
  Macrophages   = c("CD68", "MSR1", "APOE", "C1QA"),
  cDC1          = c("XCR1", "CLEC9A"),
  cDC2          = c("CD1C", "CLEC10A", "FCER1A"),
  pDCs          = c("LILRA4", "SPIB", "IRF7"),
  Endothelial   = c("PECAM1", "CDH5", "VWF"),
  SMCs          = c("COL1A2", "COL3A1", "ACTA2", "MYH11"),
  B_cells       = c("CD79A", "CD79B", "MS4A1"),
  Plasma_cells  = c("JCHAIN", "IGHG1", "IGHA1"),
  Mast_cells    = c("KIT", "CPA3"),
  Proliferative = c("MKI67"),
  Target_Gene   = c("ACKR1") 
)

all_markers <- unlist(marker_list)
valid_markers <- unique(all_markers[all_markers %in% rownames(sc_filtered)])

seen_genes <- c()
category_lengths <- numeric(length(marker_list))

for (i in seq_along(marker_list)) {
  genes_in_cat <- intersect(marker_list[[i]], rownames(sc_filtered))
  new_genes <- setdiff(genes_in_cat, seen_genes) 
  category_lengths[i] <- length(new_genes)
  seen_genes <- c(seen_genes, new_genes)
}

valid_lengths <- category_lengths[category_lengths > 0]
x_lines <- cumsum(valid_lengths)[-length(valid_lengths)] + 0.5

message(">>> 正在绘制分组气泡图...")
p_dot <- DotPlot(sc_filtered, features = valid_markers, cluster.idents = FALSE, dot.scale = 5.5) +
  scale_color_gradient2(low = "#E0F3F8", mid = "#FFFFFF", high = "#B2182B", midpoint = 0) + 
  geom_vline(xintercept = x_lines, color = "black", linetype = "dashed", linewidth = 0.5) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold.italic", size = 11, color="black"),
    axis.text.y = element_text(face = "bold", size = 12, color="black"),
    axis.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  ggtitle("Expression of Cell-Type Specific Markers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)))

ggsave(file.path(out_dir, "03_DotPlot_Annotated_Grouped.pdf"), p_dot, width = 16, height = 7)
ggsave(file.path(out_dir, "03_DotPlot_Annotated_Grouped.png"), p_dot, width = 16, height = 7, dpi = 600)

message(">>> 图表优化渲染完成！同类型细胞已完美组合在一起！")


message("\n>>> 正在提取 Cluster 7 和 10 (Endothelial Cells)...")

sc_Endo <- subset(sc_filtered, subset = seurat_clusters %in% c("7", "10"))

message("\n--- 提取完成！检查细胞数量 ---")
print(table(sc_Endo$seurat_clusters))
message(paste0(">>> 提取出的内皮细胞总数为: ", ncol(sc_Endo), " 个"))

save_path <- file.path(work_dir, "GSE253903_Endothelial_Subpopulation_C7_C10.rds")
message(paste0("\n>>> 正在将内皮亚群数据保存至: ", save_path))

saveRDS(sc_Endo, file = save_path)

# Supplementary Figure 1C-H: Differential analysis and GO/KEGG enrichment

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

work_dir <- "PATH_TO_DATA/GSE253903"
out_dir_de <- file.path(work_dir, "03_Endo_DE_Analysis")
out_dir_en <- file.path(work_dir, "04_Endo_Enrichment")

setwd(work_dir)
if(!dir.exists(out_dir_de)) dir.create(out_dir_de, recursive = TRUE)
if(!dir.exists(out_dir_en)) dir.create(out_dir_en, recursive = TRUE)

message(">>> 正在加载内皮亚群数据...")
sc_Endo <- readRDS("GSE253903_Endothelial_Subpopulation_C7_C10.rds")
Idents(sc_Endo) <- "seurat_clusters"

message(">>> 正在计算包含所有背景基因的 DEGs (C7 vs C10, logfc.threshold = 0)...")
deg_results_full <- FindMarkers(sc_Endo, 
                                ident.1 = "7", 
                                ident.2 = "10", 
                                logfc.threshold = 0.0,
                                min.pct = 0.1)

deg_results_full$gene <- rownames(deg_results_full)

logFC_cutoff <- 0.5   
padj_cutoff  <- 0.05  

deg_results_full <- deg_results_full %>%
  mutate(
    Significance = case_when(
      p_val_adj < padj_cutoff & avg_log2FC > logFC_cutoff  ~ "Up",
      p_val_adj < padj_cutoff & avg_log2FC < -logFC_cutoff ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

de_file <- file.path(out_dir_de, "DEG_C7_vs_C10_Full.csv")
write.csv(deg_results_full, de_file, row.names = FALSE)
message(">>> 全背景差异基因表格已保存至: ", de_file)


message(">>> 正在筛选火山图需要标注的显著基因...")

top_genes_df <- deg_results_full %>%
  filter(Significance %in% c("Up", "Down")) %>%
  group_by(Significance) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 15) %>%
  ungroup()

top_genes <- top_genes_df$gene

core_topic_genes <- c(
  "ACKR1", "CXCL12", "CXCR4",
  "SLC16A3", "SLC16A1", "LDHA", "LDHB", "HK2", "HIF1A",
  "EP300", "KAT2A", "HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3"
)

sig_target_genes <- deg_results_full %>%
  filter(gene %in% core_topic_genes & Significance == "Up") %>%
  pull(gene)

if (length(sig_target_genes) > 0) {
  message(">>> 🎯 命中！成功匹配到显著上调的课题核心基因：", paste(sig_target_genes, collapse = ", "))
} else {
  message(">>> 提示：在本次比较中，预设的课题核心基因未达到显著上调标准。")
}

final_label_genes <- unique(c(top_genes, sig_target_genes))

deg_results_full <- deg_results_full %>%
  mutate(
    Label = ifelse(gene %in% final_label_genes, gene, NA),
    Is_Core_Target = ifelse(gene %in% sig_target_genes, "Target", "Background")
  )


message(">>> 正在渲染 Nature 风格火山图 (高亮课题靶基因)...")

volcano_colors <- c("Up" = "#D51F26", "Not Sig" = "#D3D3D3", "Down" = "#272E6A")

p_volcano <- ggplot(deg_results_full, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(data = filter(deg_results_full, Significance == "Not Sig"), 
             aes(color = Significance), size = 1.5, alpha = 0.5, stroke = 0) +

  geom_point(data = filter(deg_results_full, Significance != "Not Sig" & Is_Core_Target == "Background"), 
             aes(color = Significance), size = 2, alpha = 0.85, stroke = 0) +

  geom_point(data = filter(deg_results_full, Is_Core_Target == "Target"), 
             fill = "#D51F26", color = "black", shape = 21, size = 4, stroke = 1.2) +

  scale_color_manual(values = volcano_colors) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  coord_cartesian(clip = "off") +

  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black", linewidth = 0.5) +

  geom_text_repel(aes(label = Label),
                  color = ifelse(deg_results_full$Is_Core_Target == "Target", "#8B0000", "black"),
                  size = ifelse(deg_results_full$Is_Core_Target == "Target", 5.5, 4.5),
                  fontface = "bold.italic", 
                  bg.color = "white", bg.r = 0.15,
                  segment.color = "black", segment.size = 0.5,
                  min.segment.length = 0,         
                  box.padding = 0.8,              
                  point.padding = 0.4,            
                  force = 5,                      
                  max.iter = 50000,               
                  max.overlaps = Inf,             
                  seed = 42) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    plot.margin = margin(t = 30, r = 30, b = 10, l = 30)
  ) +
  labs(
    title = "Volcano Plot: C7 (ACKR1+ ECs) vs C10 (ACKR1- ECs)",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

ggsave(file.path(out_dir_de, "02_Volcano_C7_vs_C10_NatureStyle.pdf"), p_volcano, width = 9, height = 7)
ggsave(file.path(out_dir_de, "02_Volcano_C7_vs_C10_NatureStyle.png"), p_volcano, width = 9, height = 7, dpi = 600)
message(">>> 绝杀！带核心靶点高亮的火山图已成功绘制！")


message(">>> 正在进行基因 ID 转换 (Symbol to EntrezID)...")

gene_id <- bitr(deg_results_full$gene, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)

deg_merged <- merge(deg_results_full, gene_id, by.x = "gene", by.y = "SYMBOL")

up_genes <- deg_merged %>% filter(Significance == "Up") %>% pull(ENTREZID)
down_genes <- deg_merged %>% filter(Significance == "Down") %>% pull(ENTREZID)

message(sprintf("   - 成功提取 %d 个上调基因 和 %d 个下调基因 用于 GO/KEGG 分析", length(up_genes), length(down_genes)))

nature_theme <- theme_classic() +
  theme(
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black")
  )

message(">>> 正在执行 C7 上调基因的 GO-BP 富集分析...")
if(length(up_genes) > 0) {
  ego_up <- enrichGO(gene = up_genes, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego_up) && nrow(ego_up) > 0) {
    write.csv(as.data.frame(ego_up), file.path(out_dir_en, "01_GO_BP_Enrichment_C7_Up.csv"), row.names = FALSE)
    p_go_up <- dotplot(ego_up, showCategory = 15) + 
      scale_color_gradient(low = "#D51F26", high = "#272E6A") + nature_theme + ggtitle("GO Biological Process: C7 (ACKR1+) Up-regulated")
    ggsave(file.path(out_dir_en, "01_GO_BP_Dotplot_C7_Up.pdf"), p_go_up, width = 8, height = 7)
    ggsave(file.path(out_dir_en, "01_GO_BP_Dotplot_C7_Up.png"), p_go_up, width = 8, height = 7, dpi = 600)
  }

  message(">>> 正在执行 C7 上调基因的 KEGG 富集分析...")
  ekegg_up <- enrichKEGG(gene = up_genes, organism = 'hsa', pvalueCutoff = 0.05)
  if (!is.null(ekegg_up) && nrow(ekegg_up) > 0) {
    ekegg_up <- setReadable(ekegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(ekegg_up), file.path(out_dir_en, "02_KEGG_Enrichment_C7_Up.csv"), row.names = FALSE)
    p_kegg_up <- dotplot(ekegg_up, showCategory = 15) + 
      scale_color_gradient(low = "#D51F26", high = "#272E6A") + nature_theme + ggtitle("KEGG Pathways: C7 (ACKR1+) Up-regulated")
    ggsave(file.path(out_dir_en, "02_KEGG_Dotplot_C7_Up.pdf"), p_kegg_up, width = 8, height = 7)
    ggsave(file.path(out_dir_en, "02_KEGG_Dotplot_C7_Up.png"), p_kegg_up, width = 8, height = 7, dpi = 600)
  }
}

message(">>> 正在执行 C7 下调 (C10 高表达) 基因的 GO-BP 富集分析...")
if(length(down_genes) > 0) {
  ego_down <- enrichGO(gene = down_genes, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego_down) && nrow(ego_down) > 0) {
    write.csv(as.data.frame(ego_down), file.path(out_dir_en, "03_GO_BP_Enrichment_C7_Down.csv"), row.names = FALSE)
    p_go_down <- dotplot(ego_down, showCategory = 15) + 
      scale_color_gradient(low = "#272E6A", high = "#8A9FD1") + nature_theme + ggtitle("GO Biological Process: C7 Down-regulated")
    ggsave(file.path(out_dir_en, "03_GO_BP_Dotplot_C7_Down.pdf"), p_go_down, width = 8, height = 7)
    ggsave(file.path(out_dir_en, "03_GO_BP_Dotplot_C7_Down.png"), p_go_down, width = 8, height = 7, dpi = 600)
  }

  message(">>> 正在执行 C7 下调 (C10 高表达) 基因的 KEGG 富集分析...")
  ekegg_down <- enrichKEGG(gene = down_genes, organism = 'hsa', pvalueCutoff = 0.05)
  if (!is.null(ekegg_down) && nrow(ekegg_down) > 0) {
    ekegg_down <- setReadable(ekegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(ekegg_down), file.path(out_dir_en, "04_KEGG_Enrichment_C7_Down.csv"), row.names = FALSE)
    p_kegg_down <- dotplot(ekegg_down, showCategory = 15) + 
      scale_color_gradient(low = "#272E6A", high = "#8A9FD1") + nature_theme + ggtitle("KEGG Pathways: C7 Down-regulated")
    ggsave(file.path(out_dir_en, "04_KEGG_Dotplot_C7_Down.pdf"), p_kegg_down, width = 8, height = 7)
    ggsave(file.path(out_dir_en, "04_KEGG_Dotplot_C7_Down.png"), p_kegg_down, width = 8, height = 7, dpi = 600)
  }
}

message("\n>>> 完美！所有差异基因分析、火山图绘制及富集分析流程已全部完成！")


# Supplementary Figure 1D: Circulation reference gene-set validation

out_dir_de <- file.path(work_dir, "03_Endo_DE_Analysis") 

message(">>> 1. 正在提取文献靶基因集...")

target_genes <- c(
  "GJA4", "GJA5", "GATA2", "MECOM",  

  "ACKR1", "NR2F2", "PLVAP",         

  "CXCL12"                           
)

valid_genes <- target_genes[target_genes %in% rownames(sc_Endo)]
message("   - 测序矩阵中成功检测到基因: ", paste(valid_genes, collapse = ", "))

Idents(sc_Endo) <- factor(sc_Endo$seurat_clusters, levels = c("7", "10"))

sc_Endo_rename <- RenameIdents(sc_Endo, 
                               "7"  = "C7 (Venular/ACKR1+ ECs)", 
                               "10" = "C10 (Arterial/ACKR1- ECs)")

message(">>> 2. 正在绘制同款 DotPlot...")

p_dot <- DotPlot(sc_Endo_rename, features = valid_genes, dot.scale = 8) +
  scale_color_gradient2(low = "#272E6A", mid = "#FFFFFF", high = "#D51F26", midpoint = 0) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold.italic", size = 13, color = "black"),
    axis.text.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title = element_blank(),

    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),

    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),

    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 15)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  ggtitle("Validation of EC Subtypes (Ref: Circulation 2025)")

pdf_path_dot <- file.path(out_dir_de, "02_Validation_Circulation_Fig1G.pdf")
png_path_dot <- file.path(out_dir_de, "02_Validation_Circulation_Fig1G.png")

ggsave(pdf_path_dot, p_dot, width = 8.5, height = 4.5)
ggsave(png_path_dot, p_dot, width = 8.5, height = 4.5, dpi = 600)

message(">>> 完美！致敬版气泡图已生成！请前往 03_Endo_DE_Analysis 查看结果。")


message(">>> 正在绘制同款 DotPlot (低表达蓝 - 高表达红)...")

p_dot <- DotPlot(sc_Endo_rename, features = valid_genes, dot.scale = 8) +
  scale_color_gradient(low = "#272E6A", high = "#D51F26") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold.italic", size = 13, color = "black"),
    axis.text.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title = element_blank(),

    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),

    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),

    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 15)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  ggtitle("Validation of EC Subtypes (Ref: Circulation 2025)")

pdf_path_dot <- file.path(out_dir_de, "02_Validation_Circulation_Fig1G.pdf")
png_path_dot <- file.path(out_dir_de, "02_Validation_Circulation_Fig1G.png")

ggsave(pdf_path_dot, p_dot, width = 8.5, height = 4.5)
ggsave(png_path_dot, p_dot, width = 8.5, height = 4.5, dpi = 600)

message(">>> 完美！实事求是版气泡图已生成并覆盖原文件！")


# Supplementary Figure 1I: Core endothelial-marker and metabolic-gene expression

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

work_dir <- "PATH_TO_DATA/GSE253903"
out_dir_ec <- file.path(work_dir, "03_Endo_DE_Analysis")
setwd(work_dir)

message(">>> 1. 正在加载内皮亚群 (C7 & C10) 单细胞数据...")
sc_Endo <- readRDS("GSE253903_Endothelial_Subpopulation_C7_C10.rds")

DefaultAssay(sc_Endo) <- "RNA"
if (inherits(sc_Endo[["RNA"]], "Assay5")) {
  sc_Endo <- JoinLayers(sc_Endo)
}

message(">>> 2. 正在匹配基因与别名...")

gene_mapping <- c(
  "ACKR1"   = "ACKR1",
  "HIF1A"   = "HIF1A",
  "CXCL12"  = "CXCL12",
  "EGLN1"   = "EGLN1 (PHD2)",
  "KDR"     = "KDR (VEGFR2)",
  "SLC16A7" = "SLC16A7 (MCT2)",
  "ENO1"    = "ENO1",
  "SLC2A3"  = "SLC2A3 (GLUT3)"
)

query_genes <- names(gene_mapping)

valid_genes <- query_genes[query_genes %in% rownames(sc_Endo)]
message("   - 成功检测到靶基因: ", paste(valid_genes, collapse = ", "))

valid_labels <- gene_mapping[valid_genes]

Idents(sc_Endo) <- factor(sc_Endo$seurat_clusters, levels = c("10", "7"))

sc_Endo_rename <- RenameIdents(sc_Endo, 
                               "10" = "C10 (ACKR1- ECs)", 
                               "7"  = "C7 (ACKR1+ ECs)")

message(">>> 3. 正在绘制对比气泡图 (带别名标签)...")

p_dot_alias <- DotPlot(sc_Endo_rename, features = valid_genes, dot.scale = 10) +
  scale_color_gradient(low = "#272E6A", high = "#D51F26") +

  scale_x_discrete(labels = valid_labels) +

  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold.italic", size = 13, color = "black"),
    axis.text.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title = element_blank(),

    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),

    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),

    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 15)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  ggtitle("Expression of Core Markers & Metabolic Genes")

pdf_path <- file.path(out_dir_ec, "06_C7_vs_C10_TargetGenes_Alias_DotPlot.pdf")
png_path <- file.path(out_dir_ec, "06_C7_vs_C10_TargetGenes_Alias_DotPlot.png")

ggsave(pdf_path, p_dot_alias, width = 9.5, height = 4.5)
ggsave(png_path, p_dot_alias, width = 9.5, height = 4.5, dpi = 600)

message(">>> 完美！带有别名标注的气泡图已生成！请前往 03_Endo_DE_Analysis 文件夹查看结果。")


# Supplementary Figure 1J: Cluster 0 vs Cluster 2 marker comparison

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

work_dir <- "PATH_TO_DATA/GSE253903"
out_dir_vis <- file.path(work_dir, "02_Annotation_and_Visualization")
setwd(work_dir)

message(">>> 正在加载单细胞数据...")
sc_obj <- readRDS("GSE253903_Integrated_Unannotated.rds") 

Idents(sc_obj) <- "seurat_clusters"

message(">>> 正在提取 0 群和 2 群并重命名...")
sc_subset <- subset(sc_obj, idents = c("0", "2"))

sc_subset <- RenameIdents(sc_subset, 
                          "0" = "0 T cells", 
                          "2" = "2 CD8+T cells")

Idents(sc_subset) <- factor(Idents(sc_subset), levels = c("0 T cells", "2 CD8+T cells"))

molecules <- c("CD3D", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7")

receptors <- c("CXCR4", "CXCR3", "CXCR6", "CCR5", "CCR7", "CCR2", "CX3CR1")

all_queries <- c(molecules, receptors)
valid_genes <- all_queries[all_queries %in% rownames(sc_subset)]

missing_genes <- setdiff(all_queries, valid_genes)
if(length(missing_genes) > 0){
  message("! 提示：以下基因在矩阵中未检测到: ", paste(missing_genes, collapse = ", "))
}

message(">>> 正在生成分子与受体表达对比气泡图...")

p_dot <- DotPlot(sc_subset, features = valid_genes, dot.scale = 10) +
  scale_color_gradient(low = "#272E6A", high = "#D51F26") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               face = "bold.italic", size = 13, color = "black"),
    axis.text.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title = element_blank(),

    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 11),

    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),

    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 15)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  geom_vline(xintercept = length(molecules[molecules %in% valid_genes]) + 0.5, 
             linetype = "dashed", color = "grey50", linewidth = 0.8) +
  labs(title = "Expression Comparison: Cluster 0 vs Cluster 2")

ggsave(file.path(out_dir_vis, "04_Cluster0_vs_2_Markers_DotPlot.pdf"), p_dot, width = 10, height = 4.5)
ggsave(file.path(out_dir_vis, "04_Cluster0_vs_2_Markers_DotPlot.png"), p_dot, width = 10, height = 4.5, dpi = 600)

message(">>> 完美！对比气泡图已生成并保存至 02_Annotation_and_Visualization 文件夹。")

