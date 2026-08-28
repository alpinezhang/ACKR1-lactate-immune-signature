# Figures 2-3: ACKR1 internal single-cell analysis (English; current revision order)
# Corresponds to current main Figures 2 (internal single-cell overview) and 3 (mechanistic, immune, and clinical-association analyses).
# Code for Figure 2C and Figures 3G-J follows generation of the endothelial/T-cell objects; panel numbers are retained in the section titles.


suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggsci)      
  library(patchwork)
  library(harmony)    
  library(DoubletFinder) 
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "01_QC_and_Clustering")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 当前工作目录已设置为: ", getwd())
message(">>> 输出图表将保存在: ", out_dir)

sample_info <- list(
  "Sample_XXXX01" = c("Patient_XXXX01", "Plaque", "Stable"),
  "Sample_XXXX02" = c("Patient_XXXX01", "Plaque", "Stable"),
  "Sample_XXXX03" = c("Patient_XXXX02", "Plaque", "Unstable"),
  "Sample_XXXX04" = c("Patient_XXXX02", "Plaque", "Unstable"),
  "Sample_XXXX05" = c("Patient_XXXX03", "Plaque", "Unstable"),
  "Sample_XXXX06" = c("Patient_XXXX03", "Plaque", "Unstable"),
  "Sample_XXXX07" = c("Patient_XXXX04", "Plaque", "Stable"),
  "Sample_XXXX08" = c("Patient_XXXX04", "Plaque", "Stable"),
  "Sample_XXXX09" = c("Patient_XXXX05", "Plaque", "Stable"),
  "Sample_XXXX10" = c("Patient_XXXX05", "Plaque", "Stable"),
  "Sample_XXXX11" = c("Patient_XXXX06", "Plaque", "Unstable"),
  "Sample_XXXX12" = c("Patient_XXXX06", "Plaque", "Unstable"),
  "Sample_XXXX13" = c("Patient_XXXX07", "Plaque", "Stable"),
  "Sample_XXXX14" = c("Patient_XXXX07", "Plaque", "Stable"),
  "Sample_XXXX15" = c("Patient_XXXX08", "Plaque", "Stable"),
  "Sample_XXXX16" = c("Patient_XXXX08", "Plaque", "Stable"),
  "Sample_XXXX17" = c("Patient_XXXX09", "Plaque", "Stable"),
  "Sample_XXXX18" = c("Patient_XXXX09", "Plaque", "Stable"),
  "Sample_XXXX19" = c("Patient_XXXX10", "Plaque", "Stable"),
  "Sample_XXXX20" = c("Patient_XXXX10", "Plaque", "Stable")
)

seurat_list <- list()
folder_names <- names(sample_info)
message(">>> 正在批量读取斑块 (Plaque) 数据...")

for (folder in folder_names) {
  if (dir.exists(folder)) {
    counts <- tryCatch({
      Read10X(data.dir = folder)
    }, error = function(e) {
      tryCatch({
        Read10X(data.dir = folder, gene.column = 1)
      }, error = function(e2) {
        return(NULL)
      })
    })

    if (is.null(counts) || length(counts) == 0) {
      message("   ! 错误: 无法解析文件夹 ", folder, " 的表达矩阵，已跳过。")
      next
    }

    if (is.list(counts) && "Gene Expression" %in% names(counts)) {
      counts <- counts$`Gene Expression`
    }

    sobj <- CreateSeuratObject(counts, project = folder)

    meta <- sample_info[[folder]]
    sobj$PatientID <- meta[1]
    sobj$Tissue    <- meta[2]
    sobj$Stability <- meta[3]
    sobj$Orig_Sample <- folder

    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
    sobj[["percent.hb"]] <- PercentageFeatureSet(sobj, pattern = "^HB[^P]")

    seurat_list[[folder]] <- sobj
    message("   - 成功读取样本: ", folder)
  } else {
    message("   ! 警告: 找不到文件夹 ", folder, "，已跳过。")
  }
}

sc_raw <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "Plaque_Project")

message(">>> 正在提取 Metadata，使用纯 ggplot2 生成质控图...")

meta <- sc_raw@meta.data
features_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb")
features_to_plot <- intersect(features_to_plot, colnames(meta))

plot_data <- data.frame()
for (feat in features_to_plot) {
  temp <- data.frame(
    Sample = meta$Orig_Sample,
    Value = meta[[feat]],
    Feature = feat
  )
  plot_data <- rbind(plot_data, temp)
}

plot_data$Feature <- factor(plot_data$Feature, levels = features_to_plot)

p_qc <- ggplot(plot_data, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.5) + 
  facet_wrap(~ Feature, scales = "free_y", ncol = length(features_to_plot)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  labs(x = "Sample ID", y = "", title = "Pre-QC Metrics across Samples")

ggsave(file.path(out_dir, "00_Pre_QC_Violin.pdf"), plot = p_qc, width = 16, height = 6)
ggsave(file.path(out_dir, "00_Pre_QC_Violin.png"), plot = p_qc, width = 16, height = 6, dpi = 300)

message(">>> 提示: 过滤前的质控图已保存至 ", file.path(out_dir, "00_Pre_QC_Violin.pdf"))
message(">>> 【请在这里暂停】打开图表，评估需要设置的 nFeature_RNA 和 percent.mt 阈值！")

min_features <- 300
max_features <- 6000
max_mt       <- 15
max_hb       <- 5

message(sprintf(">>> 正在应用针对内皮细胞优化的过滤条件: nFeature %d-%d, 线粒体 < %d%%, 红细胞 < %d%%", 
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

message(">>> 过滤完成，准备进入 DoubletFinder 去双细胞环节...")
message(">>> 正在为每个样本单独执行 DoubletFinder 去双细胞流程 (可能需要一些时间)...")

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
  message(sprintf("   - 样本 %s 去双细胞完成，保留了 %d 个单细胞", sample_name, ncol(sobj)))
}

message(">>> 正在合并去双细胞后的数据，并运行 Harmony 批次效应校正...")

sc_filtered <- merge(seurat_list[[1]], y = seurat_list[-1], 
                     add.cell.ids = names(seurat_list), project = "Plaque_Project")
rm(seurat_list, sc_raw) 
gc() 

if (inherits(sc_filtered[["RNA"]], "Assay5")) {
  message(">>> 检测到 Seurat V5，正在合并图层 (JoinLayers)...")
  sc_filtered <- JoinLayers(sc_filtered)
}

sc_filtered <- NormalizeData(sc_filtered, verbose = FALSE)
sc_filtered <- FindVariableFeatures(sc_filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sc_filtered <- ScaleData(sc_filtered, verbose = FALSE)
sc_filtered <- RunPCA(sc_filtered, npcs = 30, verbose = FALSE)

sc_filtered <- RunHarmony(sc_filtered, group.by.vars = "PatientID", plot_convergence = FALSE, verbose = FALSE)

message(">>> 正在进行 UMAP 降维与分群...")
sc_filtered <- RunUMAP(sc_filtered, reduction = "harmony", dims = 1:20, verbose = FALSE)
sc_filtered <- FindNeighbors(sc_filtered, reduction = "harmony", dims = 1:20, verbose = FALSE)
sc_filtered <- FindClusters(sc_filtered, resolution = 0.6, verbose = FALSE)

message(">>> 正在生成 SCI 级别的可视化图表...")

my_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

num_clusters <- length(levels(sc_filtered$seurat_clusters))
large_palette <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
  "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
)
custom_colors <- large_palette[1:num_clusters]

p_cluster <- DimPlot(sc_filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, pt.size = 0.1) +
  scale_color_manual(values = custom_colors) +
  ggtitle("UMAP of Plaque Cells by Cluster") +
  my_theme + theme(legend.position = "right")

ggsave(file.path(out_dir, "01_UMAP_Clusters.pdf"), p_cluster, width = 8, height = 6)
ggsave(file.path(out_dir, "01_UMAP_Clusters.png"), p_cluster, width = 8, height = 6, dpi = 300)

saveRDS(sc_filtered, file.path(work_dir, "Plaque_Integrated_Unannotated.rds"))

# Figure 2A: Nature-style UMAP

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "01_QC_and_Clustering")
setwd(work_dir)

message(">>> 正在加载单细胞数据...")
sc_filtered <- readRDS("Plaque_Integrated_Unannotated.rds")

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

message(">>> 正在进行高级图形渲染 (带图例版本)...")

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

  ggtitle("UMAP of Atherosclerotic Plaque Cells")

ggsave(file.path(out_dir, "01_UMAP_Clusters_NatureStyle_V3_Legend.pdf"), p_umap_nature, width = 9, height = 8)
ggsave(file.path(out_dir, "01_UMAP_Clusters_NatureStyle_V3_Legend.png"), p_umap_nature, width = 9, height = 8, dpi = 600)

message(">>> 带图例的 Nature 级 UMAP 渲染完成！请查看 V3 图片。")


# Figure 2B: Cell-type annotation and marker dot plot

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "01_QC_and_Clustering")
setwd(work_dir)

sc_filtered <- readRDS("Plaque_Integrated_Unannotated.rds")

message(">>> 正在为 Cluster 添加细胞类型注释...")
Idents(sc_filtered) <- "seurat_clusters" 

cluster2celltype <- c(
  "0"  = "0 T cells",
  "1"  = "1 NKT cells",
  "2"  = "2 SMCs",
  "3"  = "3 Macrophages",
  "4"  = "4 T cells",
  "5"  = "5 Macrophages",
  "6"  = "6 Monocytes",
  "7"  = "7 ACKR1- ECs",
  "8"  = "8 B cells",
  "9"  = "9 cDC2",
  "10" = "10 T cells",
  "11" = "11 ACKR1+ ECs",
  "12" = "12 Proliferating",
  "13" = "13 Mast cells",
  "14" = "14 Plasma cells",
  "15" = "15 T cells",
  "16" = "16 SMCs",
  "17" = "17 pDCs"
)

sc_filtered <- RenameIdents(sc_filtered, cluster2celltype)
sc_filtered$cell_type_annotated <- Idents(sc_filtered)

grouped_cluster_order <- c(
  "0 T cells", "4 T cells", "10 T cells", "15 T cells",
  "1 NKT cells",
  "6 Monocytes",
  "3 Macrophages", "5 Macrophages", 
  "9 cDC2",
  "17 pDCs",
  "7 ACKR1- ECs", "11 ACKR1+ ECs",
  "2 SMCs", "16 SMCs",
  "8 B cells",
  "14 Plasma cells",
  "13 Mast cells",
  "12 Proliferating"
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


suppressMessages({
  library(Seurat)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
setwd(work_dir)

sc_filtered <- readRDS("Plaque_Integrated_Unannotated.rds")

message(">>> 正在提取 Cluster 7 和 11 (Endothelial Cells)...")

sc_Endo <- subset(sc_filtered, subset = seurat_clusters %in% c("7", "11"))

message("\n--- 提取完成！检查细胞数量 ---")
print(table(sc_Endo$seurat_clusters))
message(paste0(">>> 提取出的内皮细胞总数为: ", ncol(sc_Endo), " 个"))

save_path <- file.path(work_dir, "Endothelial_Subpopulation_C7_C11.rds")
message(paste0("\n>>> 正在将内皮亚群数据保存至: ", save_path))

saveRDS(sc_Endo, file = save_path)

# Figure 2D: Differential analysis and volcano plot (C11 vs C7)
message(">>> 正在重新计算包含所有背景基因的 DEGs (logfc.threshold = 0)...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")
Idents(sc_Endo) <- "seurat_clusters"

deg_results_full <- FindMarkers(sc_Endo, 
                                ident.1 = "11", 
                                ident.2 = "7", 
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

write.csv(deg_results_full, "02_Endo_DE_Analysis/DEG_C11_vs_C7_Full.csv", row.names = FALSE)
message(">>> 全背景差异基因表格已更新！")


library(dplyr)
library(ggplot2)
library(ggrepel)

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "02_Endo_DE_Analysis")
setwd(work_dir)

logFC_cutoff <- 0.5   
padj_cutoff  <- 0.05  

message(">>> 正在加载全背景差异基因表格...")
deg_results <- read.csv("02_Endo_DE_Analysis/DEG_C11_vs_C7_Full.csv")

message(">>> 正在筛选火山图需要标注的显著基因 (Up Top 15 + Down Top 15 + 核心靶基因)...")

top_genes_df <- deg_results %>%
  filter(Significance %in% c("Up", "Down")) %>%
  group_by(Significance) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 15) %>%
  ungroup()

top_genes <- top_genes_df$gene

target_genes_to_check <- c("HIF1A", "CXCL12")

for (target_gene in target_genes_to_check) {
  if (target_gene %in% deg_results$gene) {
    gene_stats <- deg_results %>% filter(gene == target_gene)

    if (gene_stats$Significance == "Up" && !(target_gene %in% top_genes)) {
      top_genes <- c(top_genes, target_gene)
      message(sprintf("   * 检测到 %s 在 C11 中显著上调，已强制加入火山图标注队列！", target_gene))
    } else if (target_gene %in% top_genes) {
      message(sprintf("   * %s 已包含在上调的 Top 15 显著基因中。", target_gene))
    } else {
      message(sprintf("   * 提示：%s 在本次比较中未达到显著上调的阈值，将不予标注。", target_gene))
    }
  }
}

deg_results <- deg_results %>%
  mutate(Label = ifelse(gene %in% top_genes, gene, NA))

message(sprintf("   - 最终准备在火山图上标注 %d 个基因", length(top_genes)))

message(">>> 正在渲染 Nature 风格火山图...")

volcano_colors <- c("Up" = "#D51F26", "Not Sig" = "#D3D3D3", "Down" = "#272E6A")

p_volcano <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(data = filter(deg_results, Significance == "Not Sig"), 
             aes(color = Significance), size = 1.5, alpha = 0.5, stroke = 0) +
  geom_point(data = filter(deg_results, Significance != "Not Sig"), 
             aes(color = Significance), size = 2, alpha = 0.85, stroke = 0) +

  scale_color_manual(values = volcano_colors) +

  scale_x_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  coord_cartesian(clip = "off") +

  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black", linewidth = 0.5) +

  geom_text_repel(aes(label = Label),
                  size = 4.5, 
                  fontface = "bold.italic", 
                  color = "black",
                  bg.color = "white", bg.r = 0.15,
                  segment.color = "black", segment.size = 0.4,
                  min.segment.length = 0,         
                  box.padding = 0.8,              
                  point.padding = 0.4,            
                  force = 4,                      
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
    title = "Volcano Plot: C11 (ACKR1+ ECs) vs C7 (ACKR1- ECs)",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

pdf_path <- file.path(out_dir, "01_Volcano_C11_vs_C7_NatureStyle.pdf")
png_path <- file.path(out_dir, "01_Volcano_C11_vs_C7_NatureStyle.png")

ggsave(pdf_path, p_volcano, width = 9, height = 7)

# Figure 2E/G: GO and KEGG enrichment of upregulated C11 genes

suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
in_file  <- file.path(work_dir, "02_Endo_DE_Analysis", "DEG_C11_vs_C7_Full.csv")
out_dir  <- file.path(work_dir, "03_Endo_Enrichment")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 当前工作目录: ", getwd())
message(">>> 富集分析结果将保存在: ", out_dir)

message(">>> 正在加载差异基因表格...")
deg_results <- read.csv(in_file)

message(">>> 正在进行基因 ID 转换 (Symbol to EntrezID)...")

gene_id <- bitr(deg_results$gene, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)

deg_merged <- merge(deg_results, gene_id, by.x = "gene", by.y = "SYMBOL")

up_genes <- deg_merged %>% 
  filter(Significance == "Up") %>% 
  pull(ENTREZID)

message(sprintf("   - 成功提取 %d 个 C11 显著上调基因用于 GO/KEGG 分析", length(up_genes)))

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

message(">>> 正在执行 GO-BP 富集分析...")
ego <- enrichGO(gene          = up_genes,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

write.csv(as.data.frame(ego), file.path(out_dir, "01_GO_BP_Enrichment_C11_Up.csv"), row.names = FALSE)

p_go <- dotplot(ego, showCategory = 15) + 
  scale_color_gradient(low = "#D51F26", high = "#272E6A") +
  nature_theme +
  ggtitle("GO Biological Process: C11 Up-regulated")

ggsave(file.path(out_dir, "01_GO_BP_Dotplot_NatureStyle.pdf"), p_go, width = 8, height = 7)
ggsave(file.path(out_dir, "01_GO_BP_Dotplot_NatureStyle.png"), p_go, width = 8, height = 7, dpi = 600)

message(">>> 正在执行 KEGG 通路富集分析...")
ekegg <- enrichKEGG(gene          = up_genes,
                    organism      = 'hsa',
                    pvalueCutoff  = 0.05)

ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(as.data.frame(ekegg), file.path(out_dir, "02_KEGG_Enrichment_C11_Up.csv"), row.names = FALSE)

p_kegg <- dotplot(ekegg, showCategory = 15) + 
  scale_color_gradient(low = "#D51F26", high = "#272E6A") + 
  nature_theme +
  ggtitle("KEGG Pathways: C11 Up-regulated")

ggsave(file.path(out_dir, "02_KEGG_Dotplot_NatureStyle.pdf"), p_kegg, width = 8, height = 7)
ggsave(file.path(out_dir, "02_KEGG_Dotplot_NatureStyle.png"), p_kegg, width = 8, height = 7, dpi = 600)


# Figure 2F/H: GO and KEGG enrichment of downregulated C11 genes

suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
in_file  <- file.path(work_dir, "02_Endo_DE_Analysis", "DEG_C11_vs_C7_Full.csv")
out_dir  <- file.path(work_dir, "03_Endo_Enrichment")
setwd(work_dir)

message(">>> 正在加载全景差异基因表格...")
deg_results <- read.csv(in_file)

message(">>> 正在进行下调基因 ID 转换 (Symbol to EntrezID)...")

gene_id <- bitr(deg_results$gene, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)

deg_merged <- merge(deg_results, gene_id, by.x = "gene", by.y = "SYMBOL")

down_genes <- deg_merged %>% 
  filter(Significance == "Down") %>% 
  pull(ENTREZID)

message(sprintf("   - 成功提取 %d 个 C11 显著下调基因用于 GO/KEGG 分析", length(down_genes)))

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

message(">>> 正在执行下调基因的 GO-BP 富集分析...")
ego_down <- enrichGO(gene          = down_genes,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

if (!is.null(ego_down) && nrow(ego_down) > 0) {
  write.csv(as.data.frame(ego_down), file.path(out_dir, "04_GO_BP_Enrichment_C11_Down.csv"), row.names = FALSE)

  p_go_down <- dotplot(ego_down, showCategory = 15) + 
    scale_color_gradient(low = "#272E6A", high = "#8A9FD1") + 
    nature_theme +
    ggtitle("GO Biological Process: C11 Down-regulated")

  ggsave(file.path(out_dir, "04_GO_BP_Dotplot_C11_Down.pdf"), p_go_down, width = 8, height = 7)
  ggsave(file.path(out_dir, "04_GO_BP_Dotplot_C11_Down.png"), p_go_down, width = 8, height = 7, dpi = 600)
  message("   - 下调基因 GO 分析完成并已出图！(深蓝配色)")
} else {
  message("   ! 提示: 下调基因未富集到显著的 GO-BP 通路。")
}

message(">>> 正在执行下调基因的 KEGG 通路富集分析...")
ekegg_down <- enrichKEGG(gene          = down_genes,
                         organism      = 'hsa',
                         pvalueCutoff  = 0.05)

if (!is.null(ekegg_down) && nrow(ekegg_down) > 0) {
  ekegg_down <- setReadable(ekegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(ekegg_down), file.path(out_dir, "05_KEGG_Enrichment_C11_Down.csv"), row.names = FALSE)

  p_kegg_down <- dotplot(ekegg_down, showCategory = 15) + 
    scale_color_gradient(low = "#272E6A", high = "#8A9FD1") + 
    nature_theme +
    ggtitle("KEGG Pathways: C11 Down-regulated")

  ggsave(file.path(out_dir, "05_KEGG_Dotplot_C11_Down.pdf"), p_kegg_down, width = 8, height = 7)
  ggsave(file.path(out_dir, "05_KEGG_Dotplot_C11_Down.png"), p_kegg_down, width = 8, height = 7, dpi = 600)

}

# Figure 3A: HALLMARK GSEA

suppressMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell" 
out_dir  <- file.path(work_dir, "03_Endo_Enrichment")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 1. 正在加载全景差异基因数据...")
deg_res <- read.csv("02_Endo_DE_Analysis/DEG_C11_vs_C7_Full.csv")

gsea_df <- deg_res %>%
  filter(!is.na(avg_log2FC) & !is.na(gene) & gene != "") %>%
  arrange(desc(avg_log2FC))

gene_list <- gsea_df$avg_log2FC
names(gene_list) <- gsea_df$gene

message(">>> 2. 正在提取 HALLMARK 数据库并构建安全映射矩阵...")
m_df <- msigdbr(species = "Homo sapiens", category = "H")

m_t2g <- m_df %>% dplyr::select(gs_name, gene_symbol)

m_t2n <- m_df %>% 
  dplyr::select(gs_name) %>% 
  distinct() %>%
  mutate(display_name = gsub("HALLMARK_", "", gs_name))

message(">>> 3. 正在运行 GSEA 分析计算...")
set.seed(42)
gsea_hallmark <- GSEA(geneList      = gene_list, 
                      TERM2GENE     = m_t2g,
                      TERM2NAME     = m_t2n,
                      pvalueCutoff  = 0.05, 
                      pAdjustMethod = "BH",
                      minGSSize     = 10,
                      maxGSSize     = 500,
                      seed          = TRUE,
                      verbose       = FALSE)

if (!is.null(gsea_hallmark) && nrow(gsea_hallmark@result) > 0) {

  write.csv(as.data.frame(gsea_hallmark), file.path(out_dir, "03_GSEA_Hallmark_Results.csv"), row.names = FALSE)
  message(sprintf(">>> 成功！富集到 %d 条 HALLMARK 显著通路。", nrow(gsea_hallmark@result)))

  message(">>> 4. 正在绘制 GSEA 山脊图...")

  p_ridge <- suppressWarnings({
    ridgeplot(gsea_hallmark, showCategory = 15, core_enrichment = TRUE) +
      scale_fill_gradientn(colors = c("#272E6A", "#E0F3F8", "#D51F26")) + 
      theme_classic() +
      labs(title = "HALLMARK GSEA Ridgeplot: C11 vs C7", x = "Log2 Fold Change") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)),
        axis.text.y = element_text(face = "bold", size = 11, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title = element_text(face = "bold", size = 13)
      )
  })

  ggsave(file.path(out_dir, "03_GSEA_Ridgeplot_NatureStyle.pdf"), p_ridge, width = 9, height = 8)
  ggsave(file.path(out_dir, "03_GSEA_Ridgeplot_NatureStyle.png"), p_ridge, width = 9, height = 8, dpi = 600)

  message(">>> 5. 正在提取并绘制 Top 显著通路跑马图...")

  top_n <- min(4, nrow(gsea_hallmark@result))
  hit_pathways <- gsea_hallmark@result$ID[1:top_n]

  plot_title <- ifelse(top_n == 1, 
                       paste0("GSEA: ", hit_pathways), 
                       "Top Enriched HALLMARK Pathways in ACKR1+ ECs")

  nature_colors <- c("#D51F26", "#272E6A", "#00A087", "#F39B7F")

  p_lines <- gseaplot2(gsea_hallmark, 
                       geneSetID = hit_pathways, 
                       title = plot_title,
                       color = nature_colors[1:top_n], 
                       pvalue_table = TRUE,     
                       base_size = 14)          

  ggsave(file.path(out_dir, "03_GSEA_Top_Pathways_RunningScore.pdf"), p_lines, width = 11, height = 6.5)
  ggsave(file.path(out_dir, "03_GSEA_Top_Pathways_RunningScore.png"), p_lines, width = 11, height = 6.5, dpi = 600)

}

# Figure 3B: Transcription-factor activity and volcano plot

suppressMessages({
  library(Seurat)
  library(dorothea)
  library(viper)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "04_Endo_TF_Analysis")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 1. 正在加载内皮细胞亚群数据...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")
Idents(sc_Endo) <- "seurat_clusters"

message(">>> 2. 正在加载 DoRothEA 人类转录因子调控网络...")
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

message(">>> 3. 正在运行 Viper 算法推断单细胞 TF 活性 (可能需要几分钟)...")
sc_Endo <- run_viper(sc_Endo, regulon, 
                     options = list(method = "scale", minsize = 4, 
                                    eset.filter = FALSE, cores = 1, 
                                    verbose = FALSE))

message(">>> 4. 正在计算 C11 vs C7 的差异转录因子活性...")

DefaultAssay(sc_Endo) <- "dorothea"

sc_Endo <- ScaleData(sc_Endo, features = rownames(sc_Endo), verbose = FALSE)

tf_markers <- FindMarkers(sc_Endo, 
                          ident.1 = "11", 
                          ident.2 = "7", 
                          only.pos = FALSE, 
                          min.pct = 0, 
                          logfc.threshold = 0)

tf_markers$TF <- rownames(tf_markers)

tf_markers <- tf_markers %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

write.csv(tf_markers, file.path(out_dir, "01_TF_Activity_C11_vs_C7.csv"), row.names = FALSE)
message("   - 转录因子活性差异表格已保存！")

message(">>> 5. 正在绘制 Top 差异转录因子条形图...")

top_tfs <- tf_markers %>% filter(p_val_adj < 0.05) %>% top_n(15, wt = avg_log2FC)
bottom_tfs <- tf_markers %>% filter(p_val_adj < 0.05) %>% top_n(-15, wt = avg_log2FC)
plot_tfs <- rbind(top_tfs, bottom_tfs)

plot_tfs$TF <- factor(plot_tfs$TF, levels = plot_tfs$TF[order(plot_tfs$avg_log2FC)])
plot_tfs$Status <- ifelse(plot_tfs$avg_log2FC > 0, "Activated in C11", "Repressed in C11")

p_bar <- ggplot(plot_tfs, aes(x = avg_log2FC, y = TF, fill = Status)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Activated in C11" = "#D51F26", "Repressed in C11" = "#272E6A")) +
  theme_classic() +
  labs(title = "Top Regulated Transcription Factors in C11 (ACKR1+ ECs)",
       x = "TF Activity Score Difference (C11 vs C7)", y = "") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 15)),
    axis.text.y = element_text(face = "bold.italic", size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black")
  )

ggsave(file.path(out_dir, "02_TF_Activity_Barplot.pdf"), p_bar, width = 8, height = 7)
ggsave(file.path(out_dir, "02_TF_Activity_Barplot.png"), p_bar, width = 8, height = 7, dpi = 600)

message(">>> 6. 正在绘制 TF 活性热图...")

heatmap_tfs <- tf_markers %>% head(30) %>% pull(TF)

tf_matrix <- GetAssayData(sc_Endo, slot = "scale.data", assay = "dorothea")[heatmap_tfs, ]

cell_info <- data.frame(Cluster = Idents(sc_Endo))
rownames(cell_info) <- colnames(sc_Endo)

anno_colors <- list(Cluster = c("7" = "#8A9FD1", "11" = "#D51F26"))

pdf(file.path(out_dir, "03_TF_Activity_Heatmap.pdf"), width = 7, height = 8)
pheatmap(tf_matrix, 
         annotation_col = cell_info,
         annotation_colors = anno_colors,
         show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         fontsize_row = 10,
         color = colorRampPalette(c("#272E6A", "white", "#D51F26"))(100),
         main = "TF Activity Heatmap (DoRothEA)",
         border_color = NA)
dev.off()

png(file.path(out_dir, "03_TF_Activity_Heatmap.png"), width = 7*600, height = 8*600, res = 600)
pheatmap(tf_matrix, annotation_col = cell_info, annotation_colors = anno_colors,
         show_colnames = FALSE, show_rownames = TRUE, cluster_cols = TRUE, cluster_rows = TRUE,
         fontsize_row = 10, color = colorRampPalette(c("#272E6A", "white", "#D51F26"))(100),
         main = "TF Activity Heatmap (DoRothEA)", border_color = NA)
dev.off()

message("\n>>> 完美！转录因子活性分析计算与绘图全部完成！")


suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "04_Endo_TF_Analysis")
setwd(work_dir)

message(">>> 正在加载转录因子活性差异表格...")
tf_res <- read.csv("04_Endo_TF_Analysis/01_TF_Activity_C11_vs_C7.csv")

padj_cutoff <- 0.05

tf_res <- tf_res %>%
  mutate(
    Significance = case_when(
      p_val_adj < padj_cutoff & avg_log2FC > 0 ~ "Activated",
      p_val_adj < padj_cutoff & avg_log2FC < 0 ~ "Repressed",
      TRUE ~ "Not Sig"
    )
  )

message(">>> 正在筛选需要标注的 Top TF (Activated 15 + Repressed 15)...")

top_tfs_df <- tf_res %>%
  filter(Significance %in% c("Activated", "Repressed")) %>%
  group_by(Significance) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 15) %>%
  ungroup()

top_tfs <- top_tfs_df$TF

tf_res <- tf_res %>%
  mutate(Label = ifelse(TF %in% top_tfs, TF, NA))

message(sprintf("   - 最终准备在 TF 火山图上标注 %d 个转录因子", length(top_tfs)))

message(">>> 正在渲染 TF 活性火山图...")

tf_colors <- c("Activated" = "#D51F26", "Not Sig" = "#D3D3D3", "Repressed" = "#272E6A")

p_tf_volcano <- ggplot(tf_res, aes(x = avg_log2FC, y = -log10(p_val_adj))) +

  geom_point(data = filter(tf_res, Significance == "Not Sig"), 
             aes(color = Significance), size = 2, alpha = 0.5, stroke = 0) +
  geom_point(data = filter(tf_res, Significance != "Not Sig"), 
             aes(color = Significance), size = 3, alpha = 0.85, stroke = 0) +

  scale_color_manual(values = tf_colors) +

  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  coord_cartesian(clip = "off") +

  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black", linewidth = 0.5) +

  geom_text_repel(aes(label = Label),
                  size = 5, 
                  fontface = "bold.italic", 
                  color = "black",
                  bg.color = "white", bg.r = 0.15,
                  segment.color = "black", segment.size = 0.4,
                  min.segment.length = 0,         
                  box.padding = 0.8,              
                  point.padding = 0.4,            
                  force = 4,                      
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
    title = "TF Activity Volcano Plot: C11 vs C7",
    x = "TF Activity Score Difference",
    y = "-Log10(Adjusted P-value)"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

pdf_path <- file.path(out_dir, "04_TF_Activity_Volcano.pdf")
png_path <- file.path(out_dir, "04_TF_Activity_Volcano.png")

ggsave(pdf_path, p_tf_volcano, width = 9, height = 7)
ggsave(png_path, p_tf_volcano, width = 9, height = 7, dpi = 600)

# Figure 3C/D: Glycolysis and OXPHOS module scores

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "05_Glycolysis_Validation")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 正在加载内皮细胞亚群数据...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")

target_genes <- c("PFKFB3", "HK2", "PKM", "PFKM", "PFKL", "PFKP", "ALDOA", "GAPDH", "LDHA", "PDK1")

valid_genes <- target_genes[target_genes %in% rownames(sc_Endo)]
message("   - 实际检测到的糖酵解基因: ", paste(valid_genes, collapse = ", "))

Idents(sc_Endo) <- factor(sc_Endo$seurat_clusters, levels = c("11", "7"))

message(">>> 正在绘制糖酵解核心酶 DotPlot...")

p_dot <- DotPlot(sc_Endo, features = valid_genes, dot.scale = 8) +
  scale_color_gradient2(low = "#272E6A", mid = "#FFFFFF", high = "#D51F26", midpoint = 0) +

  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold.italic", size = 13, color="black"),
    axis.text.y = element_text(face = "bold", size = 13, color="black"),
    axis.title = element_blank(),

    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)),

    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),

    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),

    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  ggtitle("Expression of Core Glycolytic Enzymes (C11 vs C7)")

pdf_path <- file.path(out_dir, "01_Glycolytic_Enzymes_DotPlot.pdf")
png_path <- file.path(out_dir, "01_Glycolytic_Enzymes_DotPlot.png")

ggsave(pdf_path, p_dot, width = 8.5, height = 4.5)
ggsave(png_path, p_dot, width = 8.5, height = 4.5, dpi = 600)

message(">>> 糖酵解气泡图绘制完成！实事求是的红蓝对比图已存入 05_Glycolysis_Validation。")


suppressMessages({
  library(Seurat)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "05_Glycolysis_Validation")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 1. 正在加载内皮细胞亚群数据...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")

message(">>> 2. 正在提取 HALLMARK_OXIDATIVE_PHOSPHORYLATION 基因集...")
m_df <- msigdbr(species = "Homo sapiens", category = "H")

oxphos_genes <- m_df %>% 
  filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>% 
  pull(gene_symbol)

message(">>> 3. 正在计算每个细胞的 OXPHOS 通路打分...")
sc_Endo <- AddModuleScore(sc_Endo, 
                          features = list(oxphos_genes), 
                          name = "OXPHOS_Score")

message(">>> 4. 正在准备小提琴+箱式图数据...")
plot_data <- data.frame(
  Cluster = Idents(sc_Endo),
  OXPHOS_Score = sc_Endo$OXPHOS_Score1
)

plot_data$Cluster <- factor(plot_data$Cluster, levels = c("7", "11"))

message(">>> 5. 正在渲染 OXPHOS 打分统计图...")

max_score <- max(plot_data$OXPHOS_Score, na.rm = TRUE)

p_vln <- ggplot(plot_data, aes(x = Cluster, y = OXPHOS_Score, fill = Cluster)) +

  geom_violin(trim = FALSE, alpha = 0.85, color = "black", linewidth = 0.6) +

  geom_boxplot(width = 0.2, fill = "white", color = "black", 
               outlier.shape = NA, linewidth = 0.6) +

  scale_fill_manual(values = c("7" = "#272E6A", "11" = "#D51F26")) +

  stat_compare_means(comparisons = list(c("7", "11")),
                     method = "wilcox.test",
                     label = "p.format", 
                     label.y = max_score + 0.08,
                     size = 5, fontface = "bold") +

  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +

  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)),
    axis.text.x = element_text(face = "bold", size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r=10)),
    axis.title.x = element_blank(),
    legend.position = "none", 
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    plot.margin = margin(t = 20, r = 20, b = 10, l = 20)
  ) +
  labs(
    title = "OXPHOS Pathway Activity",
    y = "Module Score (MSigDB Hallmark)"
  )

pdf_path <- file.path(out_dir, "02_OXPHOS_Score_VlnBoxPlot.pdf")
png_path <- file.path(out_dir, "02_OXPHOS_Score_VlnBoxPlot.png")

ggsave(pdf_path, p_vln, width = 5.5, height = 6)
ggsave(png_path, p_vln, width = 5.5, height = 6, dpi = 600)

# Figure 3E: Lactylation-potential score

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "05_Glycolysis_Validation")
setwd(work_dir)

message(">>> 1. 正在加载内皮细胞亚群数据...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")

message(">>> 2. 正在构建文献衍生的 Lactylation 核心基因集...")

lactylation_genes <- c(
  "LDHA", "LDHB",
  "SLC16A1", "SLC16A3",
  "EP300", "KAT2A",
  "HDAC1", "HDAC2", "HDAC3",
  "SIRT1", "SIRT2", "SIRT3"
)

valid_lactylation_genes <- lactylation_genes[lactylation_genes %in% rownames(sc_Endo)]
message("   - 实际用于打分的基因: ", paste(valid_lactylation_genes, collapse = ", "))

message(">>> 3. 正在计算每个细胞的乳酸化 (Lactylation) 潜力打分...")
sc_Endo <- AddModuleScore(sc_Endo, 
                          features = list(valid_lactylation_genes), 
                          name = "Lactylation_Score")

plot_data <- data.frame(
  Cluster = Idents(sc_Endo),
  Lactylation_Score = sc_Endo$Lactylation_Score1 
)

plot_data$Cluster <- factor(plot_data$Cluster, levels = c("7", "11"))

message(">>> 4. 正在渲染乳酸化打分统计图...")

max_score <- max(plot_data$Lactylation_Score, na.rm = TRUE)
min_score <- min(plot_data$Lactylation_Score, na.rm = TRUE)
score_range <- max_score - min_score

p_vln <- ggplot(plot_data, aes(x = Cluster, y = Lactylation_Score, fill = Cluster)) +

  geom_violin(trim = FALSE, alpha = 0.85, color = "black", linewidth = 0.6) +
  geom_boxplot(width = 0.2, fill = "white", color = "black", 
               outlier.shape = NA, linewidth = 0.6) +

  scale_fill_manual(values = c("7" = "#272E6A", "11" = "#D51F26")) +

  stat_compare_means(comparisons = list(c("7", "11")),
                     method = "wilcox.test",
                     label = "p.format", 
                     label.y = max_score + (score_range * 0.15), 
                     size = 5, fontface = "bold") +

  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +

  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=15)),
    axis.text.x = element_text(face = "bold", size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r=10)),
    axis.title.x = element_blank(),
    legend.position = "none", 
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    plot.margin = margin(t = 20, r = 20, b = 10, l = 20)
  ) +
  labs(
    title = "Lactylation Potential Score",
    y = "Signature Score (Writers/Transporters)"
  )

pdf_path <- file.path(out_dir, "03_Lactylation_Score_VlnBoxPlot.pdf")
png_path <- file.path(out_dir, "03_Lactylation_Score_VlnBoxPlot.png")

ggsave(pdf_path, p_vln, width = 5.5, height = 6)
ggsave(png_path, p_vln, width = 5.5, height = 6, dpi = 600)

# Figure 3F: Metabolic and toxicity reprogramming (C7 vs C11)

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(msigdbr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "23_Metabolic_Exclusion_Lactate")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 1. 正在加载内皮细胞亚群数据...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")
DefaultAssay(sc_Endo) <- "RNA"

sc_Endo$seurat_clusters <- factor(sc_Endo$seurat_clusters, levels = c("7", "11"))
Idents(sc_Endo) <- "seurat_clusters"

message(">>> 2. 正在提取乳酸与脂质的权威基因集...")

m_df <- msigdbr(species = "Homo sapiens")

go_lactate <- m_df %>% filter(gs_name == "GOBP_LACTATE_METABOLIC_PROCESS") %>% pull(gene_symbol) %>% unique()
core_lactate <- c("SLC16A1", "SLC16A3", "SLC16A4", "LDHA", "LDHB", "PDHK1")
genes_lactate <- unique(c(go_lactate, core_lactate))

genes_lipid_metab <- m_df %>% filter(gs_name == "HALLMARK_FATTY_ACID_METABOLISM") %>% pull(gene_symbol) %>% unique()
tox_terms <- c("GOBP_RESPONSE_TO_OXIDIZED_LIPID", "GOBP_CELLULAR_RESPONSE_TO_LIPID_HYDROPEROXIDE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
genes_lipid_tox <- m_df %>% filter(gs_name %in% tox_terms) %>% pull(gene_symbol) %>% unique()

valid_lactate     <- intersect(genes_lactate, rownames(sc_Endo))
valid_lipid_metab <- intersect(genes_lipid_metab, rownames(sc_Endo))
valid_lipid_tox   <- intersect(genes_lipid_tox, rownames(sc_Endo))

cols_to_remove <- grep("GOLactate|HallmarkLipid|GOTox", colnames(sc_Endo@meta.data), value = TRUE)
if(length(cols_to_remove) > 0) { sc_Endo@meta.data <- sc_Endo@meta.data %>% dplyr::select(-all_of(cols_to_remove)) }

sc_Endo <- AddModuleScore(sc_Endo, features = list(valid_lactate), name = "GOLactate_Score")
sc_Endo <- AddModuleScore(sc_Endo, features = list(valid_lipid_metab), name = "HallmarkLipid_Score")
sc_Endo <- AddModuleScore(sc_Endo, features = list(valid_lipid_tox), name = "GOTox_Score")

col_lactate <- grep("GOLactate_Score", colnames(sc_Endo@meta.data), value = TRUE)[1]
col_lipidm  <- grep("HallmarkLipid_Score", colnames(sc_Endo@meta.data), value = TRUE)[1]
col_lipidt  <- grep("GOTox_Score", colnames(sc_Endo@meta.data), value = TRUE)[1]

message(">>> 3. 正在生成通路富集气泡图...")

meta_df <- sc_Endo@meta.data %>%
  dplyr::select(seurat_clusters, all_of(c(col_lactate, col_lipidm, col_lipidt)))
colnames(meta_df) <- c("Cluster", "GO: Lactate Metabolism", "HALLMARK: Lipid Metabolism", "GO: Response to ox-LDL/ROS")

meta_long <- meta_df %>%
  pivot_longer(cols = -Cluster, names_to = "Pathway", values_to = "Score") %>%
  mutate(Cluster = factor(ifelse(Cluster == "7", "C7 (ACKR1-)", "C11 (ACKR1+)"), 
                          levels = c("C7 (ACKR1-)", "C11 (ACKR1+)")))

bubble_data <- meta_long %>%
  group_by(Cluster, Pathway) %>%
  summarize(
    Mean_Score = mean(Score, na.rm = TRUE),
    Percent_Active = (sum(Score > 0) / n()) * 100,
    .groups = "drop"
  ) %>%
  group_by(Pathway) %>%
  mutate(Scaled_Mean_Score = scale(Mean_Score)[, 1]) %>%
  ungroup() %>%
  mutate(Pathway = factor(Pathway, levels = rev(c("GO: Lactate Metabolism", "HALLMARK: Lipid Metabolism", "GO: Response to ox-LDL/ROS"))))

p_bubble <- ggplot(bubble_data, aes(x = Cluster, y = Pathway)) +
  geom_point(aes(size = Percent_Active, color = Scaled_Mean_Score)) +
  scale_color_gradientn(colors = c("#272E6A", "white", "#D51F26"), 
                        name = "Relative Pathway\nActivity (Z-score)") +
  scale_size_continuous(range = c(4, 12), name = "Percent Active (%)") +
  theme_bw() +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text.x = element_text(size = 14, face = "bold", color = "black", margin = margin(t=5)),
    axis.text.y = element_text(size = 13, face = "bold.italic", color = "black", margin = margin(r=5)),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.5),
    panel.grid.major = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b=15)),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  labs(title = "Metabolic & Toxicity Reprogramming (C7 vs C11)")

ggsave(file.path(out_dir, "01_Lactate_Metabolic_Bubble_Plot.pdf"), p_bubble, width = 7.5, height = 5)
ggsave(file.path(out_dir, "01_Lactate_Metabolic_Bubble_Plot.png"), p_bubble, width = 7.5, height = 5, dpi = 600)

message(">>> 4. 正在执行 乳酸-HIF1A-CXCL12 级联驱动解耦分析...")

sc_c11 <- subset(sc_Endo, idents = "11")

plot_cor_df <- data.frame(
  CXCL12_Exp = GetAssayData(sc_c11, slot = "data")["CXCL12", ],
  HIF1A_Exp  = GetAssayData(sc_c11, slot = "data")["HIF1A", ],
  Lactate_Score = sc_c11@meta.data[[col_lactate]],
  Tox_Score     = sc_c11@meta.data[[col_lipidt]]
)

plot_cor_df <- plot_cor_df %>% filter(CXCL12_Exp > 0)

cor_long <- plot_cor_df %>%
  pivot_longer(cols = c(Lactate_Score, HIF1A_Exp, Tox_Score), 
               names_to = "Driver", values_to = "Module_Score") %>%
  mutate(Driver = factor(Driver, levels = c("Lactate_Score", "HIF1A_Exp", "Tox_Score")))

p_cor <- ggscatter(cor_long, x = "Module_Score", y = "CXCL12_Exp",
                   color = "Driver", size = 2, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "black", fill = "lightgray", size = 1),
                   conf.int = TRUE) +
  facet_wrap(~ Driver, scales = "free_x", ncol = 3,
             labeller = as_labeller(c("Lactate_Score" = "Lactate Metabolism Drive", 
                                      "HIF1A_Exp" = "HIF1A Transcriptional Drive",
                                      "Tox_Score" = "Lipid Toxicity (ox-LDL) Drive"))) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
           size = 4.5, fontface = "bold", color = "black") +
  scale_color_manual(values = c("Lactate_Score" = "#D51F26", "HIF1A_Exp" = "#D51F26", "Tox_Score" = "grey50")) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    strip.background = element_rect(fill = "#F0F0F0", color = "black", linewidth = 1),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 13),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b=15))
  ) +
  labs(title = "Mechanism Decoupling: Lactate-HIF1A Axis Dictates CXCL12 Secretion",
       x = "Metabolic Module Score / Gene Expression", 
       y = "CXCL12 Expression (Log-normalized)")

ggsave(file.path(out_dir, "02_Lactate_HIF1A_Decoupling_Scatter.pdf"), p_cor, width = 11, height = 4.5)
ggsave(file.path(out_dir, "02_Lactate_HIF1A_Decoupling_Scatter.png"), p_cor, width = 11, height = 4.5, dpi = 600)

message(">>> 极其完美的闭环！乳酸-HIF1A 级联解耦图纸已保存在 23_Metabolic_Exclusion_Lactate 中！")


# Figure 2C; Figures 3G/H: Endothelial-subtype and T-cell receptor validation

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "02_Endo_DE_Analysis") 
setwd(work_dir)

message(">>> 1. 正在加载内皮细胞亚群数据...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")

message(">>> 2. 正在提取文献靶基因集...")

target_genes <- c(
  "GJA4", "GJA5", "GATA2", "MECOM",  

  "ACKR1", "NR2F2", "PLVAP",         

  "CXCL12"                           
)

valid_genes <- target_genes[target_genes %in% rownames(sc_Endo)]
message("   - 测序矩阵中成功检测到基因: ", paste(valid_genes, collapse = ", "))

Idents(sc_Endo) <- factor(sc_Endo$seurat_clusters, levels = c("11", "7"))

sc_Endo <- RenameIdents(sc_Endo, 
                        "11" = "C11 (Venular/ACKR1+ ECs)", 
                        "7"  = "C7 (Arterial/ACKR1- ECs)")

message(">>> 3. 正在绘制同款 DotPlot...")

p_dot <- DotPlot(sc_Endo, features = valid_genes, dot.scale = 8) +
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

pdf_path <- file.path(out_dir, "04_Validation_Circulation_Fig1G.pdf")
png_path <- file.path(out_dir, "04_Validation_Circulation_Fig1G.png")

ggsave(pdf_path, p_dot, width = 8, height = 4.5)
ggsave(png_path, p_dot, width = 8, height = 4.5, dpi = 600)

message(">>> 完美！致敬版气泡图已生成！请前往 02_Endo_DE_Analysis 查看结果。")


suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "07_Tcell_Validation")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(">>> 1. 正在加载全局单细胞数据...")
sc_global <- readRDS("Plaque_Integrated_Unannotated.rds") 
Idents(sc_global) <- "seurat_clusters"

message(">>> 2. 正在提取 T 细胞亚群并进行 CD8+ 鉴定...")

sc_Tcells <- subset(sc_global, idents = c("0", "4", "10", "15"))

cd8_markers <- c("CD3D", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7")

p_cd8_id <- DotPlot(sc_Tcells, features = cd8_markers, dot.scale = 8) +
  scale_color_gradient2(low = "#272E6A", mid = "#FFFFFF", high = "#D51F26", midpoint = 0) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 12, color="black"),
    axis.text.y = element_text(face = "bold", size = 12, color="black"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black")
  ) +
  ggtitle("Identification of Cytotoxic CD8+ T cells")

ggsave(file.path(out_dir, "01_CD8_Tcell_Identification.pdf"), p_cd8_id, width = 7, height = 4)
ggsave(file.path(out_dir, "01_CD8_Tcell_Identification.png"), p_cd8_id, width = 7, height = 4, dpi=600)

message(">>> 3. 正在绘制同款趋化因子受体全景图 (升级真实标注)...")

sc_Tcells <- RenameIdents(sc_Tcells, 
                          "4"  = "CD8+ Cytotoxic T (C4)", 
                          "15" = "CD8+ Cytotoxic T (C15)",
                          "0"  = "CD4+/Other T (C0)", 
                          "10" = "CD4+/Other T (C10)")

Idents(sc_Tcells) <- factor(Idents(sc_Tcells), 
                            levels = c("CD8+ Cytotoxic T (C15)", "CD8+ Cytotoxic T (C4)", 
                                       "CD4+/Other T (C10)", "CD4+/Other T (C0)"))

chemokine_receptors <- c("CXCR4", "CXCR3", "CXCR6", "CCR5", "CCR7", "CCR2", "CX3CR1")
valid_receptors <- chemokine_receptors[chemokine_receptors %in% rownames(sc_Tcells)]

p_receptor <- DotPlot(sc_Tcells, features = valid_receptors, dot.scale = 9) +
  scale_color_gradient2(low = "#272E6A", mid = "#FFFFFF", high = "#D51F26", midpoint = 0) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 14, color="black"),
    axis.text.y = element_text(face = "bold", size = 13, color="black"),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 15))
  ) +
  ggtitle("Chemokine Receptors in T cells (Ref: Circulation Fig 1K)")

ggsave(file.path(out_dir, "02_Validation_Circulation_Fig1K_Annotated.pdf"), p_receptor, width = 8.5, height = 4.5)
ggsave(file.path(out_dir, "02_Validation_Circulation_Fig1K_Annotated.png"), p_receptor, width = 8.5, height = 4.5, dpi=600)

# Figure 3I/J: Patient-level clinical association analysis

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "08_Clinical_Relevance")
setwd(work_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

Get_Patient_Mean_Expression <- function(sobj, genes) {
  DefaultAssay(sobj) <- "RNA"
  if ("JoinLayers" %in% getNamespaceExports("Seurat")) sobj <- JoinLayers(sobj)

  valid_g <- intersect(genes, rownames(sobj))
  if(length(valid_g) == 0) return(NULL)

  mat <- GetAssayData(sobj, slot = "data", assay = "RNA")[valid_g, , drop = FALSE]
  df <- as.data.frame(t(as.matrix(mat)))
  df$PatientID <- sobj$PatientID
  df$Stability <- sobj$Stability

  patient_df <- df %>%
    group_by(PatientID, Stability) %>%
    summarise(across(all_of(valid_g), mean, na.rm = TRUE), .groups = 'drop')

  return(patient_df)
}

message(">>> 1. 正在加载内皮细胞并合并计算患者均值...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")
DefaultAssay(sc_Endo) <- "RNA"
sc_Endo$Clinical_Status <- factor(sc_Endo$Stability, levels = c("Stable", "Unstable"))
Idents(sc_Endo) <- "Clinical_Status"

endo_genes <- c("HIF1A", "CXCL12", "CXCR4", "PFKFB3", "HK2", "ALDOA", "LDHA")
valid_endo_genes <- intersect(endo_genes, rownames(sc_Endo))

if (length(valid_endo_genes) > 0) {
  p_endo_dot <- DotPlot(sc_Endo, features = valid_endo_genes, dot.scale = 10) +
    scale_color_gradient2(low = "#272E6A", mid = "#FFFFFF", high = "#D51F26", midpoint = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 13, color="black"),
          axis.text.y = element_text(face = "bold", size = 14, color="black"),
          axis.title = element_blank(), axis.line = element_line(linewidth = 1, color = "black"),
          axis.ticks = element_line(linewidth = 1, color = "black")) +
    ggtitle("Endothelial Markers: Stable vs Unstable Plaques")

  ggsave(file.path(out_dir, "02_Endo_Markers_Clinical_DotPlot.pdf"), p_endo_dot, width = 7.5, height = 4)
  ggsave(file.path(out_dir, "02_Endo_Markers_Clinical_DotPlot.png"), p_endo_dot, width = 7.5, height = 4, dpi = 600)

  patient_endo_df <- Get_Patient_Mean_Expression(sc_Endo, valid_endo_genes)
  patient_endo_long <- pivot_longer(patient_endo_df, cols = all_of(valid_endo_genes), names_to = "Gene", values_to = "Expression")
  patient_endo_long$Gene <- factor(patient_endo_long$Gene, levels = valid_endo_genes)
  patient_endo_long$Clinical_Status <- factor(patient_endo_long$Stability, levels = c("Stable", "Unstable"))

  p_endo_box <- ggplot(patient_endo_long, aes(x = Clinical_Status, y = Expression, fill = Clinical_Status)) +
    geom_boxplot(alpha = 0.8, color = "black", outlier.shape = NA, linewidth=0.8) +
    geom_jitter(width = 0.2, size = 3, color = "black", shape = 21, stroke = 1, fill = "white") +
    scale_fill_manual(values = c("Stable" = "#272E6A", "Unstable" = "#D51F26")) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
    stat_compare_means(comparisons = list(c("Stable", "Unstable")), method = "wilcox.test", label = "p.format", size=4, fontface="bold") +
    theme_classic() +
    theme(strip.text = element_text(face = "bold.italic", size = 14),
          axis.text.x = element_text(face = "bold", size = 14, color="black"),
          axis.text.y = element_text(size = 12, color="black"),
          axis.title = element_blank(), legend.position = "none",
          axis.line = element_line(linewidth = 1, color = "black"))

  ggsave(file.path(out_dir, "03_Endo_Markers_PatientLevel_Boxplot.pdf"), p_endo_box, width = 10, height = 8.5)
  ggsave(file.path(out_dir, "03_Endo_Markers_PatientLevel_Boxplot.png"), p_endo_box, width = 10, height = 8.5, dpi = 600)
}

message(">>> 运行完毕！CXCR4 已成功加入，画布尺寸已智能适配。")

message("\n>>> 2. 正在合并计算患者级别的 C4 T细胞比例...")
sc_global <- readRDS("Plaque_Integrated_Unannotated.rds") 
meta_global <- sc_global@meta.data

sample_prop <- meta_global %>%
  group_by(orig.ident, PatientID, Stability) %>%
  summarise(Total_Cells = n(), C4_Cells = sum(seurat_clusters == "4"), Prop = (C4_Cells / Total_Cells) * 100, .groups = 'drop')

patient_prop <- sample_prop %>%
  group_by(PatientID, Stability) %>%
  summarise(Patient_Mean_Prop = mean(Prop), .groups = 'drop')

patient_prop$Clinical_Status <- factor(patient_prop$Stability, levels = c("Stable", "Unstable"))

p_c4_prop <- ggplot(patient_prop, aes(x = Clinical_Status, y = Patient_Mean_Prop, fill = Clinical_Status)) +
  geom_boxplot(alpha = 0.8, color = "black", outlier.shape = NA, linewidth=0.8) +
  geom_jitter(width = 0.2, size = 3, color = "black", shape = 21, stroke = 1, fill = "white") +
  scale_fill_manual(values = c("Stable" = "#272E6A", "Unstable" = "#D51F26")) +
  stat_compare_means(comparisons = list(c("Stable", "Unstable")), method = "wilcox.test", label = "p.format", size=5, fontface="bold") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 14),
        legend.position = "none", axis.line = element_line(linewidth = 1, color = "black")) +
  labs(title = "C4 Proportion per Patient", y = "Mean C4 % of Total Cells")

ggsave(file.path(out_dir, "04_C4_Tcell_Proportion_PatientLevel.pdf"), p_c4_prop, width = 4.5, height = 5)
ggsave(file.path(out_dir, "04_C4_Tcell_Proportion_PatientLevel.png"), p_c4_prop, width = 4.5, height = 5, dpi = 600)

message(">>> 3. 正在提取并合并计算 C4 细胞的受体表达...")
sc_c4 <- subset(sc_global, subset = seurat_clusters == "4")
DefaultAssay(sc_c4) <- "RNA"
sc_c4$Clinical_Status <- factor(sc_c4$Stability, levels = c("Stable", "Unstable"))
Idents(sc_c4) <- "Clinical_Status"

c4_genes <- c("CXCR4", "CXCR3", "CCR5")
valid_c4_genes <- intersect(c4_genes, rownames(sc_c4))

if (length(valid_c4_genes) > 0) {
  p_c4_dot <- DotPlot(sc_c4, features = valid_c4_genes, dot.scale = 12) +
    scale_color_gradient2(low = "#272E6A", mid = "#FFFFFF", high = "#D51F26", midpoint = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 14, color="black"),
          axis.text.y = element_text(face = "bold", size = 14, color="black"),
          axis.title = element_blank(), axis.line = element_line(linewidth=1, color="black")) +
    ggtitle("Receptors in CD8+ T cells (C4)")
  ggsave(file.path(out_dir, "05_C4_Receptors_Clinical_DotPlot.pdf"), p_c4_dot, width = 5.5, height = 4)
  ggsave(file.path(out_dir, "05_C4_Receptors_Clinical_DotPlot.png"), p_c4_dot, width = 5.5, height = 4, dpi = 600)

  patient_c4_df <- Get_Patient_Mean_Expression(sc_c4, valid_c4_genes)
  patient_c4_long <- pivot_longer(patient_c4_df, cols = all_of(valid_c4_genes), names_to = "Gene", values_to = "Expression")
  patient_c4_long$Clinical_Status <- factor(patient_c4_long$Stability, levels = c("Stable", "Unstable"))

  p_c4_box <- ggplot(patient_c4_long, aes(x = Clinical_Status, y = Expression, fill = Clinical_Status)) +
    geom_boxplot(alpha = 0.8, color = "black", outlier.shape = NA, linewidth=0.8) +
    geom_jitter(width = 0.2, size = 3, color = "black", shape = 21, stroke = 1, fill = "white") +
    scale_fill_manual(values = c("Stable" = "#272E6A", "Unstable" = "#D51F26")) +
    facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
    stat_compare_means(comparisons = list(c("Stable", "Unstable")), method = "wilcox.test", label = "p.format", size=4, fontface="bold") +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", size = 14, color="black"),
          axis.text.y = element_text(size = 12, color="black"),
          strip.text = element_text(face = "bold.italic", size = 14),
          axis.title = element_blank(), legend.position = "none",
          axis.line = element_line(linewidth = 1, color = "black"))
  ggsave(file.path(out_dir, "06_C4_Receptors_PatientLevel_Boxplot.pdf"), p_c4_box, width = 3 * length(valid_c4_genes), height = 5)
  ggsave(file.path(out_dir, "06_C4_Receptors_PatientLevel_Boxplot.png"), p_c4_box, width = 3 * length(valid_c4_genes), height = 5, dpi = 600)

  message("\n>>> 完美通关！以患者(N=10)为单位的顶级严密统计图表已全部生成！")
}


suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "08_Clinical_Relevance")
setwd(work_dir)

sc_global <- readRDS("Plaque_Integrated_Unannotated.rds")

meta_global <- sc_global@meta.data
meta_global$Clinical_Status <- factor(meta_global$Stability, levels = c("Stable", "Unstable"))

sample_prop <- meta_global %>%
  group_by(orig.ident, PatientID, Clinical_Status) %>%
  summarise(Total_Cells = n(), 
            C4_Cells = sum(seurat_clusters == "4"), 
            Prop = (C4_Cells / Total_Cells) * 100, 
            .groups = 'drop')

patient_prop <- sample_prop %>%
  group_by(PatientID, Clinical_Status) %>%
  summarise(Patient_Mean_Prop = mean(Prop), .groups = 'drop')

p_c4_patient <- ggplot(patient_prop, aes(x = Clinical_Status, y = Patient_Mean_Prop, fill = Clinical_Status)) +
  geom_boxplot(alpha = 0.8, color = "black", outlier.shape = NA, linewidth=0.8) +
  geom_jitter(width = 0.2, size = 4, color = "black", shape = 21, stroke = 1, fill = "white") +
  scale_fill_manual(values = c("Stable" = "#272E6A", "Unstable" = "#D51F26")) +

  stat_compare_means(comparisons = list(c("Stable", "Unstable")), 
                     method = "wilcox.test", label = "p.format", size=5, fontface="bold") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14, color="black"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold", size = 14),
        legend.position = "none", 
        axis.line = element_line(linewidth = 1, color = "black")) +
  labs(title = "CD8+ T cell (C4) Proportion per Patient", y = "Mean C4 % of Total Cells")

ggsave(file.path(out_dir, "04_C4_Tcell_Proportion_PatientLevel.pdf"), p_c4_patient, width = 4.5, height = 5)
