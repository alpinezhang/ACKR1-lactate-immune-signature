# Supplementary Figures 5-6: ACKR1 revision supplementary analyses (English; current revision order)
# Corresponds to current Supplementary Figure 5 (ACKR1/glycolysis/internal single-cell communication) and Supplementary Figure 6 (revision model and supplementary validation).
# Supplementary Figure 5A-D: Final ACKR1/glycolysis redraw

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "02_Endo_DE_Analysis")
setwd(work_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

col_C7  <- "#272E6A"
col_C11 <- "#D51F26"

message(">>> 正在加载内皮细胞亚群 (C7 & C11)...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")
DefaultAssay(sc_Endo) <- "RNA"

umap_name <- grep("umap", Reductions(sc_Endo), ignore.case = TRUE, value = TRUE)[1]
if (is.na(umap_name)) stop("对象里找不到 UMAP 降维，请先 RunUMAP()。")

sc_Endo$seurat_clusters <- factor(as.character(sc_Endo$seurat_clusters),
                                  levels = c("7", "11"))
Idents(sc_Endo) <- "seurat_clusters"

message(">>> 正在绘制小提琴图...")

vln_df <- data.frame(
  cluster = sc_Endo$seurat_clusters,
  expr    = FetchData(sc_Endo, vars = "ACKR1")[, 1]
)

x_labels <- c("7" = "Cluster 7\n(ACKR1\u2212)",
              "11" = "Cluster 11\n(ACKR1+)")

p_vln <- ggplot(vln_df, aes(x = cluster, y = expr, fill = cluster)) +
  geom_jitter(width = 0.32, size = 0.35, alpha = 0.35,
              color = "grey35", show.legend = FALSE) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.9,
              color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("7" = col_C7, "11" = col_C11)) +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = expression(italic("ACKR1") ~ "Expression Distribution"),
       y = "Expression Level", x = NULL) +
  theme_classic(base_size = 13) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold.italic", size = 15,
                                margin = ggplot2::margin(b = 12)),
    axis.text.x  = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, margin = ggplot2::margin(r = 8)),
    axis.line    = element_line(linewidth = 0.9, color = "black"),
    axis.ticks   = element_line(linewidth = 0.9, color = "black"),
    legend.position = "none",
    plot.margin  = ggplot2::margin(t = 15, r = 12, b = 8, l = 15)
  )


message(">>> 正在绘制 UMAP FeaturePlot...")

emb <- as.data.frame(Embeddings(sc_Endo, umap_name))[, 1:2]
colnames(emb) <- c("UMAP_1", "UMAP_2")
emb$cluster <- sc_Endo$seurat_clusters
centroids <- emb %>%
  group_by(cluster) %>%
  summarise(x = median(UMAP_1), y = median(UMAP_2), .groups = "drop")
c11 <- centroids[centroids$cluster == "11", ]
c7  <- centroids[centroids$cluster == "7",  ]
yspan <- diff(range(emb$UMAP_2))

reds_ramp <- c("#F2F2F2", "#FDD8CC", "#FCA98C", "#F4694C",
               "#D42B20", "#A50F15")

p_umap <- FeaturePlot(sc_Endo, features = "ACKR1", reduction = umap_name,
                      order = TRUE, pt.size = 0.5) +
  scale_color_gradientn(colours = reds_ramp, name = expression(italic("ACKR1"))) +
  labs(title = expression(italic("ACKR1") ~ "Expression on Endothelial UMAP")) +
  annotate("text", x = c11$x, y = c11$y, label = "11", size = 6, fontface = "bold") +
  annotate("text", x = c7$x,  y = c7$y,  label = "7",  size = 6, fontface = "bold") +
  annotate("text", x = c11$x, y = c11$y + yspan * 0.18,
           label = "(ACKR1+)", size = 4.2, fontface = "bold") +
  annotate("text", x = c7$x,  y = c7$y  - yspan * 0.18,
           label = "(ACKR1\u2212)", size = 4.2, fontface = "bold") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 15,
                              margin = ggplot2::margin(b = 12)),
    axis.title = element_text(face = "italic", size = 12),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold.italic", size = 11)
  )

p_final <- p_vln + p_umap + plot_layout(widths = c(1, 1.15))

pdf_path <- file.path(out_dir, "ACKR1_Violin_UMAP.pdf")
png_path <- file.path(out_dir, "ACKR1_Violin_UMAP.png")
ggsave(pdf_path, p_final, width = 12, height = 5.2)
ggsave(png_path, p_final, width = 12, height = 5.2, dpi = 600)

message("\n>>> 完成！已保存：\n    ", pdf_path, "\n    ", png_path)


suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(msigdbr)
  library(patchwork)
})

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "05_Glycolysis_Validation")
setwd(work_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

col_C7  <- "#272E6A"
col_C11 <- "#D51F26"

message(">>> 正在加载内皮细胞亚群 (C7 & C11)...")
sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")
DefaultAssay(sc_Endo) <- "RNA"
Idents(sc_Endo) <- "seurat_clusters"
obj_genes <- rownames(sc_Endo)

message(">>> 正在提取糖酵解基因集...")

m_H <- msigdbr(species = "Homo sapiens", category = "H")
glyco_hallmark <- m_H %>%
  filter(gs_name == "HALLMARK_GLYCOLYSIS") %>% pull(gene_symbol) %>% unique()

glyco_go <- tryCatch({
  m_GO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  g <- m_GO %>% filter(gs_name == "GOBP_GLYCOLYTIC_PROCESS") %>%
    pull(gene_symbol) %>% unique()
  g
}, error = function(e) character(0))

core_glyco <- c("HK1","HK2","GPI","PFKL","PFKM","PFKP","PFKFB3",
                "ALDOA","ALDOB","ALDOC","TPI1","GAPDH","PGK1",
                "PGAM1","ENO1","ENO2","PKM","LDHA","SLC2A1","SLC2A3")

present <- function(genes) intersect(unique(genes), obj_genes)

set_hallmark <- present(glyco_hallmark)
set_go       <- present(glyco_go)
set_core     <- present(core_glyco)

if (length(set_go) < 5) {
  message(sprintf("   ! GO 集合交集仅 %d 个基因，第二面板改用手工核心糖酵解集。", length(set_go)))
  set2       <- set_core
  set2_label <- "Glycolysis (Core Enzymes)"
} else {
  set2       <- set_go
  set2_label <- "Glycolysis (GO:0006096)"
}

message(sprintf("   - HALLMARK 交集: %d 基因；第二集 (%s) 交集: %d 基因",
                length(set_hallmark), set2_label, length(set2)))

score_set <- function(obj, genes, name) {
  if (length(genes) < 3) {
    message("   ! ", name, " 可用基因不足 3 个，已跳过。")
    return(NULL)
  }
  obj <- AddModuleScore(obj, features = list(genes), name = name,
                        ctrl = min(100, length(obj_genes) %/% 24))
  obj[[paste0(name, "1")]][, 1]
}

message(">>> 正在计算糖酵解模块打分...")
sc_HALLMARK <- score_set(sc_Endo, set_hallmark, "Glyco_HALLMARK")
sc_SET2     <- score_set(sc_Endo, set2,         "Glyco_SET2")

score_df <- data.frame(Cluster = factor(Idents(sc_Endo), levels = c("7", "11")))
if (!is.null(sc_HALLMARK)) score_df$Glyco_HALLMARK <- sc_HALLMARK
if (!is.null(sc_SET2))     score_df$Glyco_SET2     <- sc_SET2

cat("\n================ 糖酵解模块打分统计 (C11 vs C7) ================\n")
for (col in intersect(c("Glyco_HALLMARK", "Glyco_SET2"), colnames(score_df))) {
  p   <- wilcox.test(score_df[[col]] ~ score_df$Cluster)$p.value
  m7  <- median(score_df[[col]][score_df$Cluster == "7"],  na.rm = TRUE)
  m11 <- median(score_df[[col]][score_df$Cluster == "11"], na.rm = TRUE)
  cat(sprintf("%-16s | 中位数 C7 = %.4f, C11 = %.4f | Wilcoxon P = %.3g\n",
              col, m7, m11, p))
}
cat("================================================================\n\n")

make_vln <- function(df, yvar, title) {
  max_score <- max(df[[yvar]], na.rm = TRUE)
  ggplot(df, aes(x = Cluster, y = .data[[yvar]], fill = Cluster)) +
    geom_violin(trim = FALSE, alpha = 0.85, color = "black", linewidth = 0.6) +
    geom_boxplot(width = 0.2, fill = "white", color = "black",
                 outlier.shape = NA, linewidth = 0.6) +
    scale_fill_manual(values = c("7" = col_C7, "11" = col_C11)) +
    scale_x_discrete(labels = c("7" = "C7\n(ACKR1\u2212)", "11" = "C11\n(ACKR1+)")) +
    stat_compare_means(comparisons = list(c("7", "11")),
                       method = "wilcox.test", label = "p.format",
                       label.y = max_score + 0.08, size = 5, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_classic() +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 15,
                                  margin = ggplot2::margin(b = 12)),
      axis.text.x  = element_text(face = "bold", size = 13, color = "black"),
      axis.text.y  = element_text(size = 11, color = "black"),
      axis.title.y = element_text(face = "bold", size = 13,
                                  margin = ggplot2::margin(r = 8)),
      axis.title.x = element_blank(),
      legend.position = "none",
      axis.line  = element_line(linewidth = 1, color = "black"),
      axis.ticks = element_line(linewidth = 1, color = "black"),
      plot.margin = ggplot2::margin(t = 15, r = 15, b = 8, l = 15)
    ) +
    labs(title = title, y = "Module Score")
}

plots <- list()
if ("Glyco_HALLMARK" %in% colnames(score_df))
  plots[["h"]] <- make_vln(score_df, "Glyco_HALLMARK", "Glycolysis (HALLMARK)")
if ("Glyco_SET2" %in% colnames(score_df))
  plots[["s"]] <- make_vln(score_df, "Glyco_SET2", set2_label)

p_final <- if (length(plots) == 2) plots[["h"]] + plots[["s"]] else plots[[1]]

w <- if (length(plots) == 2) 9 else 5
pdf_path <- file.path(out_dir, "SuppFig_Glycolysis_ModuleScore.pdf")
png_path <- file.path(out_dir, "SuppFig_Glycolysis_ModuleScore.png")
ggsave(pdf_path, p_final, width = w, height = 5.2)
ggsave(png_path, p_final, width = w, height = 5.2, dpi = 600)

# Supplementary Figure 5E: Internal single-cell global cell-cell communication (ACKR1+ ECs -> microenvironment)

suppressMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
})

options(future.globals.maxSize = 10000 * 1024^2)

work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "06_CellChat_Communication")
setwd(work_dir)

message(">>> 1. 正在加载全局单细胞数据...")
sc_global <- readRDS("Plaque_Integrated_Unannotated.rds")

message(">>> 2. 正在映射真实的细胞生物学名称并合并同类项...")
cell_type_mapping <- c(
  "0"  = "T cells",
  "1"  = "NKT cells",
  "2"  = "SMCs",
  "3"  = "Macrophages",
  "4"  = "T cells",
  "5"  = "Macrophages",
  "6"  = "Monocytes",
  "7"  = "ACKR1- ECs",
  "8"  = "B cells",
  "9"  = "cDC2",
  "10" = "T cells",
  "11" = "ACKR1+ ECs",
  "12" = "Proliferating",
  "13" = "Mast cells",
  "14" = "Plasma cells",
  "15" = "T cells",
  "16" = "SMCs",
  "17" = "pDCs"
)

sc_global$CellType <- as.character(cell_type_mapping[as.character(sc_global$seurat_clusters)])

Idents(sc_global) <- "CellType"

data.input <- GetAssayData(sc_global, assay = "RNA", slot = "data")
meta <- sc_global@meta.data

message(">>> 3. 正在构建 CellChat 对象...")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellType")

message(">>> 4. 正在加载 CellChat 分泌型数据库...")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

message(">>> 5. 正在进行细胞通讯概率计算 (可能需要几分钟)...")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

message("   - CellChat 基础网络构建完成！")

message(">>> 6. 正在绘制 ACKR1+ ECs 专属招募配受体气泡图...")

all_celltypes <- levels(cellchat@idents)
target_celltypes <- all_celltypes[all_celltypes != "ACKR1+ ECs"]

p_bubble <- netVisual_bubble(cellchat, 
                             sources.use = "ACKR1+ ECs", 
                             targets.use = target_celltypes, 
                             remove.isolate = FALSE,
                             title.name = "Secreted Signaling from ACKR1+ ECs to Microenvironment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face="bold"),
        axis.text.y = element_text(color = "black", face="bold.italic"))

ggsave(file.path(out_dir, "01_Sender_BubblePlot.pdf"), p_bubble, width = 11, height = 8)
ggsave(file.path(out_dir, "01_Sender_BubblePlot.png"), p_bubble, width = 11, height = 8, dpi = 600)

saveRDS(cellchat, file.path(out_dir, "CellChat_Object_Secreted_Annotated.rds"))

# Figure 7G; Supplementary Figure 6A: Fully nested LOOCV and feature stability

suppressMessages({
  library(glmnet)
  library(pROC)
})

set.seed(2024)

y_num <- as.integer(y_label == "IPH")
X     <- as.matrix(X_mat_scaled)
n     <- nrow(X)
p     <- ncol(X)

cat(sprintf(">>> 样本数 n=%d, 候选特征数 p=%d\n", n, p))
cat(sprintf(">>> 事件数: IPH=%d, non-IPH=%d\n", sum(y_num==1), sum(y_num==0)))

fit_lasso_inner <- function(Xtr, ytr) {
  cvf <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1,
                   type.measure = "deviance",
                   nfolds = min(10, length(ytr)))
  lam <- cvf$lambda.min
  fit <- glmnet(Xtr, ytr, family = "binomial", alpha = 1, lambda = lam)
  list(
    predict_fun = function(Xnew) {
      as.numeric(predict(fit, newx = Xnew, type = "response", s = lam))
    },
    selected = rownames(coef(fit))[which(as.numeric(coef(fit)) != 0)]
  )
}

cat("\n>>> [1] 正在运行完全嵌套 LOOCV ...\n")

nested_pred <- numeric(n)
feature_freq <- setNames(integer(p), colnames(X))

for (i in seq_len(n)) {
  Xtr <- X[-i, , drop = FALSE]
  ytr <- y_num[-i]
  Xte <- X[i, , drop = FALSE]

  inner <- fit_lasso_inner(Xtr, ytr)
  nested_pred[i] <- inner$predict_fun(Xte)

  sel <- setdiff(inner$selected, "(Intercept)")
  feature_freq[sel] <- feature_freq[sel] + 1L
}

roc_nested <- roc(y_num, nested_pred, quiet = TRUE, direction = "auto")
ci_nested  <- ci.auc(roc_nested)

cat(sprintf("   - 嵌套 LOOCV AUC = %.3f (95%% CI: %.3f - %.3f)\n",
            as.numeric(roc_nested$auc), ci_nested[1], ci_nested[3]))

cat("\n   - 各候选特征在 LOOCV 各折中被 LASSO 选中的频率（稳定性）:\n")
ff <- sort(feature_freq, decreasing = TRUE)
for (nm in names(ff)) {
  cat(sprintf("       %-16s : %d / %d 折 (%.0f%%)\n",
              nm, ff[nm], n, 100 * ff[nm] / n))
}

cat("\n>>> [2] 正在运行 Bootstrap 乐观度校正 (B=1000) ...\n")

B <- 1000

app_model <- fit_lasso_inner(X, y_num)
app_pred  <- app_model$predict_fun(X)
apparent_auc <- as.numeric(roc(y_num, app_pred, quiet = TRUE, direction = "auto")$auc)
cat(sprintf("   - Apparent AUC (乐观, 未校正) = %.3f\n", apparent_auc))

optimism_vec <- numeric(B)
b_done <- 0
for (b in seq_len(B)) {
  idx <- sample(seq_len(n), n, replace = TRUE)
  if (length(unique(y_num[idx])) < 2) next
  Xb <- X[idx, , drop = FALSE]
  yb <- y_num[idx]

  boot_model <- tryCatch(fit_lasso_inner(Xb, yb), error = function(e) NULL)
  if (is.null(boot_model)) next

  pred_boot_on_boot <- boot_model$predict_fun(Xb)
  auc_boot <- tryCatch(as.numeric(roc(yb, pred_boot_on_boot, quiet = TRUE, direction = "auto")$auc),
                       error = function(e) NA)

  pred_boot_on_orig <- boot_model$predict_fun(X)
  auc_orig <- tryCatch(as.numeric(roc(y_num, pred_boot_on_orig, quiet = TRUE, direction = "auto")$auc),
                       error = function(e) NA)

  if (is.na(auc_boot) || is.na(auc_orig)) next
  b_done <- b_done + 1
  optimism_vec[b_done] <- auc_boot - auc_orig
}
optimism_vec <- optimism_vec[seq_len(b_done)]

mean_optimism  <- mean(optimism_vec, na.rm = TRUE)
corrected_auc  <- apparent_auc - mean_optimism

cat(sprintf("   - 有效 Bootstrap 次数 = %d\n", b_done))
cat(sprintf("   - 平均乐观度 (optimism) = %.3f\n", mean_optimism))
cat(sprintf("   - 乐观度校正后 AUC (optimism-corrected) = %.3f\n", corrected_auc))

cat("\n============================================================\n")
cat(" 汇总 (可直接填入回复信 XX 处)\n")
cat("------------------------------------------------------------\n")
cat(sprintf(" Apparent AUC (未校正, 乐观)        : %.3f\n", apparent_auc))
cat(sprintf(" Fully nested LOOCV AUC             : %.3f (95%% CI %.3f-%.3f)\n",
            as.numeric(roc_nested$auc), ci_nested[1], ci_nested[3]))
cat(sprintf(" Bootstrap optimism                 : %.3f\n", mean_optimism))
cat(sprintf(" Optimism-corrected AUC             : %.3f\n", corrected_auc))
cat("============================================================\n")

res_summary <- data.frame(
  Metric = c("Apparent_AUC", "Nested_LOOCV_AUC", "Nested_LOOCV_CI_low",
             "Nested_LOOCV_CI_high", "Bootstrap_optimism", "Optimism_corrected_AUC"),
  Value  = c(apparent_auc, as.numeric(roc_nested$auc), ci_nested[1],
             ci_nested[3], mean_optimism, corrected_auc)
)
print(res_summary)


suppressMessages({
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(cowplot)
})

if(!exists("out_dir")) out_dir <- getwd()

theme_nature <- function() {
  theme_classic() +
    theme(text = element_text(color = "black", family = "sans"),
          axis.text = element_text(size = 13, color = "black", face = "bold"),
          axis.title = element_text(size = 14, face = "bold", color = "black"),
          axis.line = element_line(linewidth = 1.2, color = "black"),
          axis.ticks = element_line(linewidth = 1.2, color = "black"),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5,
                                    margin = ggplot2::margin(b = 5)),
          plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0.5,
                                       margin = ggplot2::margin(b = 12)))
}

auc_nested <- as.numeric(roc_nested$auc)

roc_df <- data.frame(
  spec = rev(roc_nested$specificities),
  sens = rev(roc_nested$sensitivities)
)

label_txt <- sprintf(
  "Apparent AUC = %.3f\nNested LOOCV AUC = %.3f (95%% CI: %.3f-%.3f)\nOptimism-corrected AUC = %.3f",
  apparent_auc, auc_nested, ci_nested[1], ci_nested[3], corrected_auc
)

p_roc_nested <- ggplot(roc_df, aes(x = spec, y = sens)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 1) +
  geom_line(color = "#E64B35", linewidth = 2.5) +
  scale_x_reverse(limits = c(1, 0), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  annotate("text", x = 0.62, y = 0.18, label = label_txt,
           hjust = 0, size = 4.2, fontface = "bold", color = "black") +
  labs(title = "Internal Validation: Fully Nested LOOCV",
       subtitle = sprintf("Dataset: GSE163154 | N = %d", n),
       x = "Specificity", y = "Sensitivity") +
  theme_nature()

print(p_roc_nested)
ggsave(file.path(out_dir, "Figure9G_Nested_LOOCV_ROC.pdf"), p_roc_nested, width = 7, height = 7)
ggsave(file.path(out_dir, "Figure9G_Nested_LOOCV_ROC.png"), p_roc_nested, width = 7, height = 7, dpi = 600)

if (exists("feature_freq")) {
  stab_df <- data.frame(
    Feature = names(feature_freq),
    Pct = 100 * as.numeric(feature_freq) / n
  ) %>%
    filter(Feature != "(Intercept)") %>%
    arrange(desc(Pct))

  stab_df$Feature <- fct_reorder(stab_df$Feature, stab_df$Pct)

  core_feats <- c("ACKR1", "SLC16A3", "SLC16A4", "CD8_Tox_Score")
  stab_df$Core <- ifelse(stab_df$Feature %in% core_feats, "Core", "Other")

  p_stab <- ggplot(stab_df, aes(x = Pct, y = Feature, fill = Core)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f%%", Pct)),
              hjust = -0.15, size = 4.2, fontface = "bold") +
    scale_fill_manual(values = c("Core" = "#E64B35", "Other" = "grey70"), guide = "none") +
    scale_x_continuous(limits = c(0, 110), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    labs(title = "Feature Selection Stability across Nested LOOCV",
         subtitle = sprintf("Frequency selected by LASSO over %d folds", n),
         x = "Selection Frequency (%)", y = "") +
    theme_nature()

  print(p_stab)
  ggsave(file.path(out_dir, "Figure9H_Feature_Stability.pdf"), p_stab, width = 7.5, height = 6)
  ggsave(file.path(out_dir, "Figure9H_Feature_Stability.png"), p_stab, width = 7.5, height = 6, dpi = 600)

  p_merge <- plot_grid(p_roc_nested, p_stab, labels = c("A", "B"),
                       label_size = 18, ncol = 2, align = "h", rel_widths = c(1, 1))
  ggsave(file.path(out_dir, "Figure9_Nested_Combined.pdf"), p_merge, width = 15, height = 6.5)
  ggsave(file.path(out_dir, "Figure9_Nested_Combined.png"), p_merge, width = 15, height = 6.5, dpi = 600)

  cat(">>> 已输出: Figure9G_Nested_LOOCV_ROC / Figure9H_Feature_Stability / Figure9_Nested_Combined\n")
} else {
  cat(">>> 未找到 feature_freq 对象，仅输出 ROC 图。请重跑嵌套 LOOCV 那段以生成 feature_freq。\n")
}

# Supplementary Figure 6B: Cell-composition adjustment and partial correlations

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


library(GEOquery)
library(illuminaHumanv2.db)
library(limma)

gset <- getGEO(filename = "GSE163154_series_matrix.txt.gz", getGPL = FALSE)
expr_gse <- exprs(gset)
pheno_gse <- pData(gset)

group_gse <- ifelse(grepl("non-IPH", pheno_gse$title, ignore.case = TRUE), "non_IPH", "IPH")
group_factor <- factor(group_gse, levels = c("non_IPH", "IPH"))
cat(">>> 分组:", table(group_factor), "\n")

sym <- mapIds(illuminaHumanv2.db,
              keys = rownames(expr_gse),
              column = "SYMBOL", keytype = "PROBEID",
              multiVals = "first")
keep <- !is.na(sym)
expr_sym <- expr_gse[keep, ]
sym <- sym[keep]

expr_bulk <- limma::avereps(expr_sym, ID = sym)
expr_bulk <- as.matrix(expr_bulk)

cat(">>> expr_bulk 就绪:", nrow(expr_bulk), "基因 x", ncol(expr_bulk), "样本\n")


exists("sc_filtered"); exists("expr_bulk"); exists("group_factor")
identical(ncol(expr_bulk), length(group_factor))


suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ppcor)
})

set.seed(2024)
if(!exists("out_dir")) out_dir <- getwd()

glyco_genes   <- c("HK2","PFKL","PFKM","PKM","ALDOA","GAPDH","LDHA","PFKFB3")
lactate_genes <- c("SLC16A3","SLC16A4","SLC16A7")
immune_genes  <- c("CD8A","GZMA","GZMB","PRF1")
key_genes     <- c("ACKR1","CXCL12","HIF1A", glyco_genes, lactate_genes, immune_genes)

message(">>> [分析一] 构建 pseudo-bulk 并在 ACKR1+ EC 内部验证 ...")

DefaultAssay(sc_filtered) <- "RNA"
counts <- GetAssayData(sc_filtered, layer = "counts")

meta <- sc_filtered@meta.data
meta$cell_type <- as.character(meta$cell_type_annotated)
meta$sample    <- as.character(meta$Orig_Sample)

make_pseudobulk <- function(cell_idx, genes) {
  sub_counts <- counts[genes, cell_idx, drop = FALSE]
  sm <- meta$sample[cell_idx]
  agg <- t(apply(sub_counts, 1, function(g) tapply(g, sm, sum)))
  lib <- tapply(Matrix::colSums(counts[, cell_idx, drop = FALSE]), sm, sum)
  cpm <- sweep(agg, 2, lib[colnames(agg)], "/") * 1e6
  log2(cpm + 1)
}

idx_ackr1ec <- which(meta$seurat_clusters == "11")
n_samp_ec <- length(unique(meta$sample[idx_ackr1ec]))
message(sprintf("    ACKR1+ EC 覆盖 %d 个样本", n_samp_ec))

pb_ec <- make_pseudobulk(idx_ackr1ec, intersect(key_genes, rownames(counts)))
pb_ec <- pb_ec[, colSums(is.na(pb_ec)) == 0, drop = FALSE]

cor_within_ec <- data.frame()
target_genes <- intersect(c(glyco_genes, lactate_genes, immune_genes, "CXCL12","HIF1A"),
                          rownames(pb_ec))
for (g in target_genes) {
  x <- pb_ec["ACKR1", ]; y <- pb_ec[g, ]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) >= 4 && sd(x[ok])>0 && sd(y[ok])>0) {
    ct <- cor.test(x[ok], y[ok], method = "spearman")
    cor_within_ec <- rbind(cor_within_ec,
                           data.frame(Gene=g, R=unname(ct$estimate), P=ct$p.value, n=sum(ok)))
  }
}
message("    ACKR1+ EC 内部 ACKR1 vs 机制分子 相关性 (pseudo-bulk, 单一细胞类型内):")
print(cor_within_ec)
write.csv(cor_within_ec, file.path(out_dir, "PseudoBulk_within_ACKR1EC_correlation.csv"), row.names=FALSE)

message("\n>>> [分析二] 解卷积估计细胞比例 + 偏相关校正 ...")

big_map <- c(
  "0 T cells"="T","4 T cells"="T","10 T cells"="T","15 T cells"="T","1 NKT cells"="T",
  "3 Macrophages"="Mac","5 Macrophages"="Mac","6 Monocytes"="Mac","9 cDC2"="Mac","17 pDCs"="Mac",
  "2 SMCs"="SMC","16 SMCs"="SMC",
  "7 ACKR1- ECs"="EC","11 ACKR1+ ECs"="EC",
  "8 B cells"="B","14 Plasma cells"="B","13 Mast cells"="Mast","12 Proliferating"="Prolif"
)
meta$big <- big_map[as.character(meta$cell_type_annotated)]

common_genes <- intersect(rownames(counts), rownames(expr_bulk))
sc_logcpm <- log2(sweep(counts[common_genes,], 2, Matrix::colSums(counts), "/")*1e6 + 1)

sig <- sapply(unique(na.omit(meta$big)), function(bt){
  cells <- which(meta$big == bt)
  Matrix::rowMeans(sc_logcpm[, cells, drop=FALSE])
})
sig_var <- apply(sig, 1, function(x) max(x) - mean(x))
top_sig_genes <- names(sort(sig_var, decreasing = TRUE))[1:min(1500, length(sig_var))]
sig_mat <- sig[top_sig_genes, ]

library(nnls)
bulk_sub <- expr_bulk[top_sig_genes, , drop=FALSE]
props <- apply(bulk_sub, 2, function(y){
  fit <- nnls(sig_mat, y)
  p <- fit$x; p/sum(p)
})
props <- t(props)
colnames(props) <- colnames(sig_mat)
message("    估计的平均细胞比例:")
print(round(colMeans(props), 3))
write.csv(props, file.path(out_dir, "Deconvolution_cellproportions.csv"))

message("\n    偏相关 (控制 EC/T/Mac 细胞比例) vs 普通相关:")
ctrl <- as.data.frame(props[, intersect(c("EC","T","Mac","SMC"), colnames(props))])

pairs_to_test <- list(
  c("ACKR1","HK2"), c("ACKR1","SLC16A3"), c("ACKR1","SLC16A4"),
  c("ACKR1","CXCL12"), c("ACKR1","CD8A"), c("ACKR1","GZMB"),
  c("HIF1A","CXCL12")
)
res_pcor <- data.frame()
for (pr in pairs_to_test) {
  g1 <- pr[1]; g2 <- pr[2]
  if (!(g1 %in% rownames(expr_bulk) && g2 %in% rownames(expr_bulk))) next
  x <- as.numeric(expr_bulk[g1, rownames(ctrl)])
  y <- as.numeric(expr_bulk[g2, rownames(ctrl)])
  raw <- cor.test(x, y, method = "pearson")
  pc  <- pcor.test(x, y, ctrl, method = "pearson")
  res_pcor <- rbind(res_pcor, data.frame(
    Pair = paste(g1, g2, sep = " vs "),
    Raw_R = unname(raw$estimate), Raw_P = raw$p.value,
    Partial_R = pc$estimate, Partial_P = pc$p.value
  ))
}
print(res_pcor)
write.csv(res_pcor, file.path(out_dir, "PartialCorrelation_adjusted_composition.csv"), row.names=FALSE)

cat("\n============================================================\n")
cat(" 汇总：细胞组成校正结果 (可填入回复 XX 处)\n")
cat("------------------------------------------------------------\n")
cat(" [A] Pseudo-bulk (仅 ACKR1+ EC 内部, 单细胞类型内相关):\n")
for(i in seq_len(nrow(cor_within_ec))){
  cat(sprintf("     ACKR1 vs %-8s : R=%.2f, P=%.2g (n=%d)\n",
              cor_within_ec$Gene[i], cor_within_ec$R[i], cor_within_ec$P[i], cor_within_ec$n[i]))
}
cat("\n [B] Bulk 偏相关 (校正 EC/T/Mac 比例后):\n")
for(i in seq_len(nrow(res_pcor))){
  cat(sprintf("     %-16s : 原始 R=%.2f(P=%.2g) -> 偏相关 R=%.2f(P=%.2g)\n",
              res_pcor$Pair[i], res_pcor$Raw_R[i], res_pcor$Raw_P[i],
              res_pcor$Partial_R[i], res_pcor$Partial_P[i]))
}
cat("============================================================\n")


suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
})
select <- dplyr::select

if(!exists("out_dir")) out_dir <- getwd()

theme_nature <- function() {
  theme_classic() +
    theme(text = element_text(color="black", family="sans"),
          axis.text = element_text(size=13, color="black", face="bold"),
          axis.title = element_text(size=14, face="bold", color="black"),
          axis.line = element_line(linewidth=1.2, color="black"),
          axis.ticks = element_line(linewidth=1.2, color="black"),
          plot.title = element_text(size=15, face="bold", hjust=0.5,
                                    margin = ggplot2::margin(b=4)),
          plot.subtitle = element_text(size=11, color="grey30", hjust=0.5,
                                       margin = ggplot2::margin(b=10)),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size=12, face="bold"),
          panel.grid.major.y = element_blank())
}

if(!exists("res_pcor")) {
  res_pcor <- read.csv(file.path(out_dir, "PartialCorrelation_adjusted_composition.csv"))
}

res_pcor$Sig_after <- ifelse(res_pcor$Partial_P < 0.05, "*", "")

plotA_df <- res_pcor %>%
  dplyr::select(Pair, Raw_R, Partial_R) %>%
  pivot_longer(cols = c(Raw_R, Partial_R),
               names_to = "Type", values_to = "R") %>%
  mutate(Type = recode(Type,
                       Raw_R = "Unadjusted (raw)",
                       Partial_R = "Adjusted (partial)"),
         Type = factor(Type, levels = c("Unadjusted (raw)", "Adjusted (partial)")))

order_pairs <- res_pcor %>% arrange(desc(abs(Raw_R))) %>% dplyr::pull(Pair)
plotA_df$Pair <- factor(plotA_df$Pair, levels = rev(order_pairs))

lab_map <- setNames(
  paste0(res_pcor$Pair, ifelse(res_pcor$Sig_after=="*", "  *", "")),
  res_pcor$Pair)

pA <- ggplot(plotA_df, aes(x = R, y = Pair, fill = Type)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66,
           color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linewidth = 0.9, color = "grey40") +
  scale_fill_manual(values = c("Unadjusted (raw)" = "#4C72B0",
                               "Adjusted (partial)" = "#C44E52")) +
  scale_y_discrete(labels = lab_map) +
  scale_x_continuous(limits = c(-0.9, 0.9), breaks = seq(-0.8, 0.8, 0.4)) +
  labs(title = "Correlations before vs after adjusting for cell composition",
       subtitle = "GSE163154 | Pearson partial correlation (* = still significant after adjustment)",
       x = "Correlation coefficient (R)", y = "") +
  theme_nature()

print(pA)
ggsave(file.path(out_dir, "FigureSX_Composition_Adjustment.pdf"), pA, width = 8.5, height = 6)
ggsave(file.path(out_dir, "FigureSX_Composition_Adjustment.png"), pA, width = 8.5, height = 6, dpi = 600)

# Supplementary Figure 6H: Metabolic module scores across cell types

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
})

if(!exists("out_dir")) out_dir <- getwd()
Idents(sc_filtered) <- "cell_type_annotated"
DefaultAssay(sc_filtered) <- "RNA"

glycolysis_genes <- list(c("HK2","GPI","PFKL","PFKM","PFKP","ALDOA","ALDOC",
                           "GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA","PFKFB3",
                           "SLC2A1","PGM1"))
hypoxia_genes    <- list(c("HIF1A","VEGFA","SLC2A1","CA9","PGK1","LDHA","BNIP3",
                           "NDRG1","ADM","PDK1","ANKRD37","EGLN3","P4HA1"))
lactate_genes    <- list(c("LDHA","LDHB","SLC16A1","SLC16A3","SLC16A4","PDK1"))

sc_filtered <- AddModuleScore(sc_filtered, features = glycolysis_genes,
                              name = "Glycolysis_Score", seed = 2024)
sc_filtered <- AddModuleScore(sc_filtered, features = hypoxia_genes,
                              name = "Hypoxia_Score", seed = 2024)
sc_filtered <- AddModuleScore(sc_filtered, features = lactate_genes,
                              name = "Lactate_Score", seed = 2024)

score_df <- FetchData(sc_filtered,
                      vars = c("cell_type_annotated",
                               "Glycolysis_Score1","Hypoxia_Score1","Lactate_Score1"))
colnames(score_df) <- c("CellType","Glycolysis","Hypoxia","Lactate")

mean_tab <- score_df %>%
  group_by(CellType) %>%
  summarise(across(c(Glycolysis,Hypoxia,Lactate), mean), .groups="drop")
cat(">>> 各细胞类型平均模块评分:\n")
print(as.data.frame(mean_tab %>% mutate(across(where(is.numeric), ~round(.,3)))))
write.csv(mean_tab, file.path(out_dir, "CellType_Metabolic_ModuleScores.csv"), row.names=FALSE)

long_df <- score_df %>%
  pivot_longer(cols = c(Glycolysis,Hypoxia,Lactate),
               names_to = "Module", values_to = "Score")
long_df$Module <- factor(long_df$Module, levels = c("Glycolysis","Hypoxia","Lactate"))

long_df$Highlight <- case_when(
  grepl("ACKR1\\+", long_df$CellType) ~ "ACKR1+ ECs",
  grepl("Macrophage", long_df$CellType) ~ "Macrophages",
  TRUE ~ "Other"
)

order_ct <- mean_tab %>% arrange(desc(Glycolysis)) %>% pull(CellType)
long_df$CellType <- factor(long_df$CellType, levels = order_ct)

theme_nature <- function() {
  theme_classic() +
    theme(text=element_text(color="black"),
          axis.text.y=element_text(size=10,color="black",face="bold"),
          axis.text.x=element_text(size=10,color="black",face="bold"),
          axis.title=element_text(size=12,face="bold"),
          strip.text=element_text(size=12,face="bold"),
          strip.background=element_rect(fill="grey92",color=NA),
          plot.title=element_text(size=14,face="bold",hjust=0.5),
          legend.position="top", legend.title=element_blank(),
          legend.text=element_text(size=11,face="bold"))
}

p <- ggplot(long_df, aes(x = Score, y = CellType, fill = Highlight)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.35) +
  facet_wrap(~Module, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c("ACKR1+ ECs" = "#C44E52",
                               "Macrophages" = "#4C72B0",
                               "Other" = "grey80")) +
  labs(title = "Metabolic module scores across all cell types",
       x = "Module score", y = "") +
  theme_nature()

print(p)
ggsave(file.path(out_dir, "FigureSX_ModuleScore_byCellType.pdf"), p, width = 12, height = 6)
ggsave(file.path(out_dir, "FigureSX_ModuleScore_byCellType.png"), p, width = 12, height = 6, dpi = 600)

# Supplementary Figure 6C: GSE253903 CD8+ T-cell proportion

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

sc_filtered <- readRDS(file.path("PATH_TO_DATA/GSE253903", "GSE253903_Integrated_Unannotated.rds"))

if(!"cell_type_annotated" %in% colnames(sc_filtered@meta.data)) {
  Idents(sc_filtered) <- "seurat_clusters"
  cluster2celltype <- c(
    "0"="0 T cells", "1"="1 Monocyte", "2"="2 CD8+T cells", "3"="3 pDCs",
    "4"="4 SMCs", "5"="5 Macrophages", "6"="6 Monocyte", "7"="7 ACKR1+ ECs",
    "8"="8 T cells", "9"="9 Proliferating", "10"="10 ACKR1- ECs",
    "11"="11 Plasma cells", "12"="12 Mast cells", "13"="13 SMCs",
    "14"="14 SMCs", "15"="15 B cells", "16"="16 T cells"
  )
  sc_filtered <- RenameIdents(sc_filtered, cluster2celltype)
  sc_filtered$cell_type_annotated <- Idents(sc_filtered)
}

meta <- sc_filtered@meta.data
cat(">>> 症状分组样本数:\n"); print(table(meta$Condition))
cat("\n>>> 细胞类型:\n"); print(table(meta$cell_type_annotated))

cd8_label <- "2 CD8+T cells"
meta$is_cd8 <- as.character(meta$cell_type_annotated) == cd8_label

prop_tab <- meta %>%
  filter(!is.na(Condition)) %>%
  group_by(Orig_Sample, Condition) %>%
  summarise(total = n(),
            cd8_n = sum(is_cd8),
            cd8_pct = 100 * mean(is_cd8),
            .groups = "drop")

cat("\n>>> 各样本 CD8+ 毒性 T 占比:\n")
print(as.data.frame(prop_tab))

wt <- wilcox.test(cd8_pct ~ Condition, data = prop_tab)
summ <- prop_tab %>% group_by(Condition) %>%
  summarise(mean_pct = round(mean(cd8_pct),2),
            median_pct = round(median(cd8_pct),2),
            n = n(), .groups="drop")
cat("\n>>> Symptomatic vs Asymptomatic:\n"); print(as.data.frame(summ))
cat(sprintf("    Wilcoxon P = %.3f\n", wt$p.value))

sym_mean   <- summ$mean_pct[summ$Condition=="Symptomatic"]
asym_mean  <- summ$mean_pct[summ$Condition=="Asymptomatic"]
direction  <- ifelse(sym_mean > asym_mean, "有症状组更高(方向一致)", "有症状组更低(方向不一致)")
cat(sprintf("    >>> 结论方向: %s (Sym=%.2f%% vs Asym=%.2f%%)\n", direction, sym_mean, asym_mean))

prop_tab$Condition <- factor(prop_tab$Condition, levels=c("Asymptomatic","Symptomatic"))
p <- ggplot(prop_tab, aes(x=Condition, y=cd8_pct, fill=Condition)) +
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75) +
  geom_jitter(width=0.12, size=3, shape=21, color="black") +
  scale_fill_manual(values=c("Asymptomatic"="#4C72B0","Symptomatic"="#C44E52")) +
  labs(title="CD8+ cytotoxic T cell proportion (GSE253903)",
       subtitle=sprintf("Symptomatic vs Asymptomatic | Wilcoxon P = %.3f", wt$p.value),
       x="", y="CD8+ cytotoxic T cells (% of total)") +
  theme_classic() +
  theme(legend.position="none",
        axis.text=element_text(size=12,face="bold",color="black"),
        axis.title=element_text(size=13,face="bold"),
        plot.title=element_text(size=14,face="bold",hjust=0.5))
print(p)
ggsave("FigureSX_GSE253903_CD8_by_symptom.png", p, width=5, height=6, dpi=600)
ggsave("FigureSX_GSE253903_CD8_by_symptom.pdf", p, width=5, height=6)

# Supplementary Figure 6D-F: Multiple-testing and expression validation of key genes

suppressMessages(library(dplyr))

stopifnot(exists("deg_clean"))

key_genes <- c(
  "ACKR1", "HIF1A", "CXCL12", "CXCR4",
  "HK2", "PFKL", "PFKM", "PKM", "ALDOA", "ALDOC", "GAPDH", "PFKFB3",
  "LDHA", "LDHB", "SLC16A3", "SLC16A4", "SLC16A7",
  "CD8A", "GZMA", "GZMB", "PRF1", "IFNG"
)

check_tab <- deg_clean %>%
  filter(Gene %in% key_genes) %>%
  dplyr::select(Gene, logFC, P.Value, adj.P.Val) %>%
  mutate(
    raw_sig = ifelse(P.Value < 0.05, "*", ""),
    fdr_sig = ifelse(adj.P.Val < 0.05, "*", ""),
    P.Value   = signif(P.Value, 3),
    adj.P.Val = signif(adj.P.Val, 3),
    logFC     = round(logFC, 3)
  ) %>%
  arrange(adj.P.Val)

cat("========================================================\n")
cat(" 关键基因：原始 P 值 vs FDR 校正后 adj.P.Val\n")
cat("--------------------------------------------------------\n")
print(as.data.frame(check_tab), row.names = FALSE)
cat("========================================================\n")

n_total   <- nrow(check_tab)
n_fdr_sig <- sum(check_tab$adj.P.Val < 0.05)
n_raw_sig <- sum(check_tab$P.Value < 0.05)
cat(sprintf("\n找到 %d 个关键基因：\n", n_total))
cat(sprintf("  原始 P<0.05 显著: %d 个\n", n_raw_sig))
cat(sprintf("  FDR<0.05  显著: %d 个\n", n_fdr_sig))

lost <- check_tab %>% filter(P.Value < 0.05 & adj.P.Val >= 0.05)
if(nrow(lost) > 0){
  cat("\n>>> 以下基因原始P显著、但FDR校正后不显著（回复需如实标注）:\n")
  print(as.data.frame(lost %>% dplyr::select(Gene, P.Value, adj.P.Val)), row.names = FALSE)
} else {
  cat("\n>>> 所有原始P显著的关键基因，FDR校正后依然显著。\n")
}

write.csv(check_tab, "KeyGenes_FDR_check.csv", row.names = FALSE)
cat("\n>>> 已保存 KeyGenes_FDR_check.csv\n")


suppressMessages({
  library(dplyr)
  library(openxlsx)
})

stopifnot(exists("deg_clean"))

gene_module <- tribble(
  ~Gene,      ~Category,
  "ACKR1",    "Core target",
  "HIF1A",    "Core target",
  "CXCL12",   "Core target",
  "CXCR4",    "Chemokine receptor",
  "HK2",      "Glycolysis",
  "PFKL",     "Glycolysis",
  "PFKM",     "Glycolysis",
  "PKM",      "Glycolysis",
  "ALDOA",    "Glycolysis",
  "ALDOC",    "Glycolysis",
  "GAPDH",    "Glycolysis",
  "PFKFB3",   "Glycolysis",
  "LDHA",     "Lactate metabolism",
  "LDHB",     "Lactate metabolism",
  "SLC16A3",  "Lactate transport",
  "SLC16A4",  "Lactate transport",
  "SLC16A7",  "Lactate transport",
  "CD8A",     "CD8+ T cytotoxicity",
  "GZMA",     "CD8+ T cytotoxicity",
  "GZMB",     "CD8+ T cytotoxicity",
  "PRF1",     "CD8+ T cytotoxicity",
  "IFNG",     "CD8+ T cytotoxicity"
)

fmt_p <- function(p){
  ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
         formatC(p, format = "f", digits = 3))
}

sup_tab <- deg_clean %>%
  inner_join(gene_module, by = "Gene") %>%
  transmute(
    Category,
    Gene,
    log2FC        = round(logFC, 3),
    `Raw P value` = fmt_p(P.Value),
    `FDR (adj.P)` = fmt_p(adj.P.Val),
    `Significant (FDR<0.05)` = ifelse(adj.P.Val < 0.05, "Yes", "No")
  ) %>%
  arrange(factor(Category, levels = c(
    "Core target","Chemokine receptor","Glycolysis",
    "Lactate metabolism","Lactate transport","CD8+ T cytotoxicity")),
    `FDR (adj.P)`)

cat(">>> 补充表 SX 预览：\n")
print(as.data.frame(sup_tab), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Key genes FDR")
writeData(wb, "Key genes FDR", sup_tab, headerStyle =
            createStyle(textDecoration = "bold", halign = "center",
                        fgFill = "#D9E1F2", border = "TopBottom"))
setColWidths(wb, "Key genes FDR", cols = 1:6, widths = c(20,10,10,14,14,20))
freezePane(wb, "Key genes FDR", firstRow = TRUE)
saveWorkbook(wb, "Supplementary_Table_SX_KeyGenes_FDR.xlsx", overwrite = TRUE)

write.csv(sup_tab, "Supplementary_Table_SX_KeyGenes_FDR.csv", row.names = FALSE)

cat("\n>>> 已输出 Supplementary_Table_SX_KeyGenes_FDR.xlsx / .csv\n")


suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

if(!exists("sc_Endo")) sc_Endo <- readRDS("Endothelial_Subpopulation_C7_C11.rds")
Idents(sc_Endo) <- "seurat_clusters"
DefaultAssay(sc_Endo) <- "RNA"

iph_genes <- c("BLVRB", "HMOX1", "CD164")
iph_genes <- iph_genes[iph_genes %in% rownames(sc_Endo)]
cat(">>> 可用 IPH 标志物:", paste(iph_genes, collapse=", "), "\n\n")

expr <- FetchData(sc_Endo, vars = c(iph_genes, "seurat_clusters"))
expr$Group <- ifelse(expr$seurat_clusters == "11", "ACKR1+ EC (C11)",
                     ifelse(expr$seurat_clusters == "7", "ACKR1- EC (C7)", NA))
expr <- expr[!is.na(expr$Group), ]

summ <- expr %>%
  group_by(Group) %>%
  summarise(across(all_of(iph_genes),
                   list(mean = ~round(mean(.),3),
                        pct  = ~round(mean(.>0)*100,1)),
                   .names = "{.col}_{.fn}"), .groups="drop")
cat(">>> ACKR1+ EC vs ACKR1- EC 的 IPH 标志物表达:\n")
print(as.data.frame(summ), row.names = FALSE)

cat("\n>>> FindMarkers (C11 vs C7) 差异检验:\n")
mk <- FindMarkers(sc_Endo, ident.1 = "11", ident.2 = "7",
                  features = iph_genes, logfc.threshold = 0, min.pct = 0)
mk$gene <- rownames(mk)
print(mk[, c("gene","avg_log2FC","p_val","p_val_adj","pct.1","pct.2")])
cat("\n   注: avg_log2FC>0 表示在 ACKR1+ EC(C11)中更高; pct.1=C11阳性率, pct.2=C7阳性率\n")

cat("\n>>> 方向小结:\n")
for(g in iph_genes){
  row <- mk[mk$gene==g, ]
  if(nrow(row)==0){ cat(sprintf("   %s: 未检出\n", g)); next }
  dir <- ifelse(row$avg_log2FC > 0, "ACKR1+ EC 更高", "ACKR1- EC 更高")
  sig <- ifelse(row$p_val_adj < 0.05, "显著", "不显著")
  cat(sprintf("   %-8s: %s (log2FC=%.2f, adj.P=%.2g, %s)\n",
              g, dir, row$avg_log2FC, row$p_val_adj, sig))
}

library(tidyr)
long <- expr %>% pivot_longer(all_of(iph_genes), names_to="Gene", values_to="Expr")
p <- ggplot(long, aes(Group, Expr, fill=Group)) +
  geom_violin(alpha=0.5, trim=FALSE) +
  geom_boxplot(width=0.2, outlier.shape=NA) +
  facet_wrap(~Gene, scales="free_y") +
  scale_fill_manual(values=c("ACKR1- EC (C7)"="#4C72B0","ACKR1+ EC (C11)"="#C44E52")) +
  labs(title="Intraplaque hemorrhage markers: ACKR1+ vs ACKR1- ECs", x="", y="Expression") +
  theme_classic() +
  theme(legend.position="none", axis.text.x=element_text(angle=20, hjust=1, face="bold"))
ggsave("IPH_markers_ACKR1EC_vs_ACKR1negEC.png", p, width=8, height=4.5, dpi=600)
ggsave("IPH_markers_ACKR1EC_vs_ACKR1negEC.pdf", p, width=8, height=4.5)
cat("\n>>> 已输出图 IPH_markers_ACKR1EC_vs_ACKR1negEC\n")


# Supplementary Figure 6G: ACKR1+ endothelial signature validation in bulk

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

if(!exists("deg_results_full")){
  deg_results_full <- read.csv("02_Endo_DE_Analysis/DEG_C11_vs_C7_Full.csv")
}
if(!"gene" %in% colnames(deg_results_full)) deg_results_full$gene <- rownames(deg_results_full)

ackr1ec_sig <- deg_results_full %>%
  filter(avg_log2FC > 0.5, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

cat(sprintf(">>> 单细胞 ACKR1+ EC 上调特征基因数: %d\n", length(ackr1ec_sig)))
cat("    (前20个):", paste(head(ackr1ec_sig,20), collapse=", "), "\n\n")

sig_in_bulk <- intersect(ackr1ec_sig, rownames(expr_bulk))
cat(sprintf(">>> 其中在 bulk 中可检测到的: %d 个\n", length(sig_in_bulk)))

mat <- expr_bulk[sig_in_bulk, , drop=FALSE]
mat_z <- t(scale(t(mat)))
sig_score <- colMeans(mat_z, na.rm=TRUE)

score_df <- data.frame(Sample = names(sig_score),
                       Score = sig_score,
                       Group = group_factor)

wt <- wilcox.test(Score ~ Group, data = score_df)
cat("\n>>> ACKR1+ EC signature 评分 IPH vs non-IPH:\n")
print(score_df %>% group_by(Group) %>%
        summarise(mean=round(mean(Score),3), median=round(median(Score),3), n=n()))
cat(sprintf("    Wilcoxon P = %.3g\n", wt$p.value))

iph_idx  <- which(group_factor == "IPH")
niph_idx <- which(group_factor == "non_IPH")
up_count <- 0; test_tab <- data.frame()
for(g in sig_in_bulk){
  v_iph  <- as.numeric(expr_bulk[g, iph_idx])
  v_niph <- as.numeric(expr_bulk[g, niph_idx])
  lfc <- mean(v_iph) - mean(v_niph)
  p   <- tryCatch(wilcox.test(v_iph, v_niph)$p.value, error=function(e) NA)
  if(!is.na(lfc) && lfc > 0) up_count <- up_count + 1
  test_tab <- rbind(test_tab, data.frame(gene=g, bulk_diff=round(lfc,3), p=signif(p,3)))
}
cat(sprintf("\n>>> %d 个可检测特征基因中, 在 IPH 组表达更高(上调)的有 %d 个 (%.0f%%)\n",
            length(sig_in_bulk), up_count, 100*up_count/length(sig_in_bulk)))
n_sig_up <- sum(test_tab$bulk_diff>0 & test_tab$p<0.05, na.rm=TRUE)
cat(sprintf("    其中显著上调(P<0.05)的: %d 个\n", n_sig_up))

cat("\n>>> 在 bulk IPH 组上调最明显的特征基因(前15):\n")
print(head(test_tab %>% arrange(desc(bulk_diff)), 15), row.names=FALSE)

write.csv(test_tab, "ACKR1EC_signature_in_bulk_IPH.csv", row.names=FALSE)

score_df$Group <- factor(score_df$Group, levels=c("non_IPH","IPH"))
p <- ggplot(score_df, aes(Group, Score, fill=Group)) +
  geom_violin(alpha=0.5, trim=FALSE) +
  geom_boxplot(width=0.2, outlier.shape=NA) +
  geom_jitter(width=0.1, size=2, alpha=0.6) +
  scale_fill_manual(values=c("non_IPH"="#4C72B0","IPH"="#C44E52")) +
  labs(title="ACKR1+ EC signature score in bulk (GSE163154)",
       subtitle=sprintf("IPH vs non-IPH | Wilcoxon P = %.3g", wt$p.value),
       x="", y="ACKR1+ EC signature score") +
  theme_classic() +
  theme(legend.position="none", axis.text=element_text(face="bold",size=12))
ggsave("ACKR1EC_signature_bulk_IPH.png", p, width=5, height=6, dpi=600)
ggsave("ACKR1EC_signature_bulk_IPH.pdf", p, width=5, height=6)

# Supplementary Figure 5F: Targeted recruitment of CD8+ T cells (C4) by ACKR1+ ECs

suppressMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
})

options(future.globals.maxSize = 10000 * 1024^2)
work_dir <- "PATH_TO_DATA/ACKR1_single_cell"
out_dir  <- file.path(work_dir, "06_CellChat_Communication")
setwd(work_dir)

message(">>> 1. 正在加载全局单细胞数据以分离 C4...")
sc_global <- readRDS("Plaque_Integrated_Unannotated.rds")

cell_type_mapping <- c(
  "0"  = "Other T cells",
  "1"  = "NKT cells",
  "2"  = "SMCs",
  "3"  = "Macrophages",
  "4"  = "CD8+ T cells (C4)",
  "5"  = "Macrophages",
  "6"  = "Monocytes",
  "7"  = "ACKR1- ECs",
  "8"  = "B cells",
  "9"  = "cDC2",
  "10" = "Other T cells",
  "11" = "ACKR1+ ECs",
  "12" = "Proliferating",
  "13" = "Mast cells",
  "14" = "Plasma cells",
  "15" = "Other CD8+ T",
  "16" = "SMCs",
  "17" = "pDCs"
)

sc_global$FineCellType <- as.character(cell_type_mapping[as.character(sc_global$seurat_clusters)])
Idents(sc_global) <- "FineCellType"

message(">>> 2. 正在开启高灵敏度模式计算细胞通讯网络...")
data.input <- GetAssayData(sc_global, assay = "RNA", slot = "data")
meta <- sc_global@meta.data

cc_targeted <- createCellChat(object = data.input, meta = meta, group.by = "FineCellType")
cc_targeted@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cc_targeted <- subsetData(cc_targeted)

cc_targeted <- identifyOverExpressedGenes(cc_targeted, thresh.pc = 0.01, thresh.fc = 1.05, thresh.p = 0.1)
cc_targeted <- identifyOverExpressedInteractions(cc_targeted)

cc_targeted <- computeCommunProb(cc_targeted, type = "truncatedMean", trim = 0.05)

cc_targeted <- filterCommunication(cc_targeted, min.cells = 5)

cc_targeted <- computeCommunProbPathway(cc_targeted)
cc_targeted <- aggregateNet(cc_targeted)

message(">>> 当前高灵敏度网络中捕捉到的通路有：")
print(cc_targeted@netP$pathways)
message(">>> 3. 正在绘制 1V1 靶向气泡图...")

p_bubble_1v1 <- netVisual_bubble(cc_targeted, 
                                 sources.use = "ACKR1+ ECs", 
                                 targets.use = "CD8+ T cells (C4)", 
                                 remove.isolate = FALSE,
                                 title.name = "Targeted Signaling: ACKR1+ ECs -> CD8+ T cells (C4)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face="bold", size = 12),
        axis.text.y = element_text(color = "black", face="bold.italic", size = 12))

ggsave(file.path(out_dir, "04_Targeted_Bubble_C11_to_C4.pdf"), p_bubble_1v1, width = 5, height = 7)
ggsave(file.path(out_dir, "04_Targeted_Bubble_C11_to_C4.png"), p_bubble_1v1, width = 5, height = 7, dpi = 600)
