# Figure 9: ACKR1 spatial transcriptomics analysis (English; current revision order)
# Corresponds to current main Figure 9 (spatial transcriptomics).
# GSE100927 external validation is in script 03; this file retains only the spatial-transcriptomics code for current Figure 9.


suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

work_dir <- "PATH_TO_DATA/ACKR1_spatial"
setwd(work_dir)
out_dir  <- file.path(work_dir, "02_Spatial_Analysis")
proc_dir <- file.path(out_dir, "processed")
for (d in c(out_dir, proc_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

sample_group <- c(
  Sample_XXXX01 = "Unknown", Sample_XXXX02 = "Unknown", Sample_XXXX03 = "Unknown",
  Sample_XXXX04 = "Unknown", Sample_XXXX05 = "Unknown", Sample_XXXX06 = "Unknown",
  Sample_XXXX07 = "Unknown", Sample_XXXX08 = "Unknown", Sample_XXXX09 = "Unknown",
  Sample_XXXX10 = "Unknown"
)

ec_markers <- c("PECAM1", "CDH5", "VWF", "CLDN5", "FLT1", "EGFL7", "RAMP2", "PLVAP")
hb_genes   <- c("HBA1","HBA2","HBB","HBD","HBM","HBQ1","HBG1","HBG2","HBZ","HBE1")

get_counts <- function(obj, assay) {
  tryCatch(GetAssayData(obj, assay = assay, layer = "counts"),
           error = function(e) GetAssayData(obj, assay = assay, slot = "counts"))
}

filter_files <- list.files(work_dir, pattern = "_filter\\.rds$")
filter_files <- filter_files[!grepl("^\\._", filter_files)]
sample_ids   <- gsub("_filter\\.rds$", "", filter_files)
message(">>> 待处理样本: ", paste(sample_ids, collapse = ", "))

all_summary <- list()

for (sid in sample_ids) {
  message("\n========== 处理样本: ", sid, " ==========")
  obj  <- readRDS(file.path(work_dir, paste0(sid, "_filter.rds")))
  meta <- readRDS(file.path(work_dir, paste0(sid, "_meta.rds")))

  main_assay <- Assays(obj)[which.max(sapply(Assays(obj), function(a) nrow(obj[[a]])))]
  cnt <- get_counts(obj, main_assay)
  std <- CreateSeuratObject(counts = cnt, assay = "Spatial", project = sid)
  rm(obj); gc()

  if (all(colnames(std) %in% rownames(meta))) {
    meta <- meta[colnames(std), , drop = FALSE]
  } else if (nrow(meta) == ncol(std)) {
    message("   ! barcode 名对不上，按相同行序兜底对齐")
    rownames(meta) <- colnames(std)
  } else {
    stop("meta 行数与 spot 数不一致，无法对齐，请检查 ", sid)
  }
  add_cols <- setdiff(colnames(meta), colnames(std@meta.data))
  std <- AddMetaData(std, metadata = meta[, add_cols, drop = FALSE])

  std$Sample <- sid
  std$Group  <- ifelse(sid %in% names(sample_group), sample_group[[sid]], "Unknown")

  if (all(c("x","y") %in% colnames(std@meta.data))) {
    std$coord_x <- as.numeric(std$x); std$coord_y <- as.numeric(std$y)
  } else stop("样本 ", sid, " 缺少 x/y 坐标列")

  std[["percent.mt"]] <- PercentageFeatureSet(std, pattern = "^MT-")
  hb_in <- intersect(hb_genes, rownames(std))
  std[["percent.hb"]] <- if (length(hb_in) > 0) PercentageFeatureSet(std, features = hb_in) else 0

  std <- NormalizeData(std, verbose = FALSE)

  ec_in <- intersect(ec_markers, rownames(std))
  if (length(ec_in) >= 2) {
    std <- tryCatch(AddModuleScore(std, features = list(ec_in), name = "EC_score", ctrl = 50, assay = "Spatial"),
                    error = function(e) { std$EC_score1 <- NA_real_; std })
    ec_score <- std$EC_score1
  } else ec_score <- NA_real_

  dat <- GetAssayData(std, assay = "Spatial", layer = "data")
  ackr1_expr <- if ("ACKR1" %in% rownames(dat)) as.numeric(dat["ACKR1", ]) else NA_real_

  saveRDS(std, file.path(proc_dir, paste0(sid, "_std.rds")))

  m <- std@meta.data
  all_summary[[sid]] <- data.frame(
    Sample      = sid,
    Group       = m$Group,
    barcode     = colnames(std),
    x           = m$coord_x,
    y           = m$coord_y,
    nCount      = m$nCount_Spatial,
    nFeature    = m$nFeature_Spatial,
    percent.mt  = m$percent.mt,
    percent.hb  = m$percent.hb,
    predicted.id = if ("predicted.id" %in% colnames(m)) as.character(m$predicted.id) else NA_character_,
    EC_score    = ec_score,
    ACKR1_expr  = ackr1_expr,
    stringsAsFactors = FALSE
  )
  message(sprintf("   - %s 完成：%d spots，ACKR1 %s", sid, ncol(std),
                  ifelse(all(is.na(ackr1_expr)), "未检出(C2)", "已记录")))
  rm(std, dat); gc()
}

summary_df <- do.call(rbind, all_summary)
saveRDS(summary_df, file.path(out_dir, "AllSpots_Summary.rds"))
write.csv(summary_df, file.path(out_dir, "AllSpots_Summary.csv"), row.names = FALSE)
message("\n>>> 总表已保存：", nrow(summary_df), " 个 spot")

nature_pal <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F",
                "#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")

qc_long <- summary_df %>%
  dplyr::select(Sample, nFeature, nCount, percent.mt, percent.hb) %>%
  tidyr::pivot_longer(-Sample, names_to = "Feature", values_to = "Value")
qc_long$Feature <- factor(qc_long$Feature,
                          levels = c("nFeature","nCount","percent.mt","percent.hb"))
qc_long$Sample  <- factor(qc_long$Sample, levels = sample_ids)

p_qc <- ggplot(qc_long, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.4, color = "black") +
  facet_wrap(~ Feature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = nature_pal) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y  = element_text(color = "black"),
    strip.text   = element_text(face = "bold", size = 13),
    strip.background = element_blank(),
    legend.position = "none",
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 15)
  ) +
  labs(x = NULL, y = NULL, title = "Spatial QC Metrics across 10 Samples (post-filter)")

ggsave(file.path(out_dir, "01_QC_AllSamples_Violin.pdf"), p_qc, width = 15, height = 4.5)
ggsave(file.path(out_dir, "01_QC_AllSamples_Violin.png"), p_qc, width = 15, height = 4.5, dpi = 300)

message("\n>>> 第 1 步完成！")
message(">>> QC 图：", file.path(out_dir, "01_QC_AllSamples_Violin.pdf"))
message(">>> 请打开 QC 图看一眼 percent.mt / percent.hb 分布，告诉我是否需要补做过滤；")
message(">>> 然后我们进第 2 步：定义 ACKR1+ 内皮 spot 并画 10 个患者的空间定位图。")


# Figure 9: Publication-quality colocalization and proximity analysis
suppressMessages({ library(Seurat); library(ggplot2); library(dplyr) })
use_rast <- requireNamespace("ggrastr", quietly = TRUE)
use_fnn  <- requireNamespace("FNN", quietly = TRUE)

work_dir <- "PATH_TO_DATA/ACKR1_spatial"
out_dir  <- file.path(work_dir, "02_Spatial_Analysis")
proc_dir <- file.path(out_dir, "processed")
pp_dir   <- file.path(out_dir, "per_patient"); if(!dir.exists(pp_dir)) dir.create(pp_dir, recursive=TRUE)
hb_cut   <- 20; flip_y <- FALSE; min_spots <- 5

get_data  <- function(o) tryCatch(GetAssayData(o,"Spatial",layer="data"),
                                  error=function(e) GetAssayData(o,"Spatial",slot="data"))
safe_expr <- function(dat,g) if (g %in% rownames(dat)) as.numeric(dat[g,]) else rep(0, ncol(dat))
nn_dist   <- function(q,r){ if(use_fnn) FNN::get.knnx(r,q,k=1)$nn.dist[,1]
  else apply(q,1,function(p) sqrt(min((r[,1]-p[1])^2+(r[,2]-p[2])^2))) }

COL <- c("Other"="grey90", "CD8 T cells"="#2C6FA6", "ACKR1+ EC"="#C0392B")
SZ  <- c("Other"=0.28, "CD8 T cells"=0.85, "ACKR1+ EC"=0.85)

theme_nat <- function(base=12){
  theme_classic(base_size=base) +
    theme(plot.title=element_text(hjust=.5, face="bold", size=base+2),
          plot.subtitle=element_text(hjust=.5, color="grey35", size=base-1),
          axis.text=element_text(color="black"),
          axis.line=element_line(linewidth=.5, color="black"),
          axis.ticks=element_line(color="black"),
          legend.title=element_text(face="bold"),
          strip.text=element_text(face="bold", size=base),
          strip.background=element_blank())
}
theme_sp <- theme_void(base_size=12) +
  theme(plot.title=element_text(hjust=.5,face="bold",size=14),
        plot.subtitle=element_text(hjust=.5,size=10,color="grey35"),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=11), strip.text=element_text(face="bold",size=12))
rast <- function(g) if(use_rast) ggrastr::rasterise(g, dpi=400) else g

std_files <- list.files(proc_dir, pattern="_std\\.rds$"); sample_ids <- gsub("_std\\.rds$","",std_files)
enr <- list()
for (sid in sample_ids) {
  std <- readRDS(file.path(proc_dir, paste0(sid,"_std.rds"))); m <- std@meta.data; dat <- get_data(std)
  pid <- if ("predicted.id" %in% colnames(m)) as.character(m$predicted.id) else NA_character_
  is_EC <- grepl("endothel", pid, ignore.case=TRUE); is_EC[is.na(is_EC)] <- FALSE
  ackr1 <- safe_expr(dat,"ACKR1")
  cd3 <- pmax(safe_expr(dat,"CD3D"), safe_expr(dat,"CD3E"))
  cd8 <- pmax(safe_expr(dat,"CD8A"), safe_expr(dat,"CD8B"))
  hb  <- m$percent.hb
  enr[[sid]] <- data.frame(Sample=sid, x=m$coord_x, y=m$coord_y,
                           is_ackr1_EC = is_EC & ackr1>0 & !(hb>=hb_cut),
                           is_cd8t     = cd3>0 & cd8>0, stringsAsFactors=FALSE)
  message("   - 提取: ", sid); rm(std,dat); gc()
}
df <- do.call(rbind, enr)
if (flip_y) df$y <- -df$y
df$Sample <- factor(df$Sample, levels=sort(unique(df$Sample)))
df$Cat <- "Other"; df$Cat[df$is_cd8t] <- "CD8 T cells"; df$Cat[df$is_ackr1_EC] <- "ACKR1+ EC"
df$Cat <- factor(df$Cat, levels=c("Other","CD8 T cells","ACKR1+ EC"))
saveRDS(df, file.path(out_dir,"Coloc_markerbased.rds"))

d_ord <- df %>% arrange(Cat)
p_co <- ggplot(d_ord, aes(x,y,color=Cat,size=Cat)) + rast(geom_point(stroke=0)) +
  scale_color_manual(values=COL) + scale_size_manual(values=SZ, guide="none") +
  facet_wrap(~Sample, scales="free", nrow=2) +
  guides(color=guide_legend(override.aes=list(size=3.2))) +
  labs(title="Spatial co-localization of ACKR1+ EC and CD8 T cells") + theme_sp
ggsave(file.path(out_dir,"06b_Coloc_AllPatients.pdf"), p_co, width=16, height=7)
ggsave(file.path(out_dir,"06b_Coloc_AllPatients.png"), p_co, width=16, height=7, dpi=300)

for (sid in levels(df$Sample)) {
  d <- df %>% filter(Sample==sid) %>% arrange(Cat)
  p <- ggplot(d, aes(x,y,color=Cat,size=Cat)) + rast(geom_point(stroke=0)) +
    scale_color_manual(values=COL) + scale_size_manual(values=SZ, guide="none") +
    coord_fixed() + guides(color=guide_legend(override.aes=list(size=3.2))) +
    labs(title=paste0("Patient ",sid),
         subtitle=sprintf("ACKR1+ EC: %d  |  CD8 T: %d", sum(d$is_ackr1_EC), sum(d$is_cd8t))) +
    theme_sp
  ggsave(file.path(pp_dir,paste0(sid,"_Coloc_v2.pdf")), p, width=7, height=6)
  ggsave(file.path(pp_dir,paste0(sid,"_Coloc_v2.png")), p, width=7, height=6, dpi=300)
}
message(">>> 共定位图 06b 完成")

dist_all <- list(); stat_rows <- list()
for (sid in levels(df$Sample)) {
  d <- df %>% filter(Sample==sid)
  ec <- as.matrix(d[d$is_ackr1_EC,c("x","y")]); td <- as.matrix(d[d$is_cd8t,c("x","y")])
  if (nrow(ec)<min_spots || nrow(td)<min_spots) next
  set.seed(1); rnd <- as.matrix(d[sample(nrow(d), nrow(td)),c("x","y")])
  dc <- nn_dist(td,ec); dr <- nn_dist(rnd,ec); nm <- median(dr)
  dist_all[[sid]] <- rbind(data.frame(Sample=sid,Group="CD8 T cells",dist=dc,rel=dc/nm),
                           data.frame(Sample=sid,Group="Random spot",dist=dr,rel=dr/nm))
  stat_rows[[sid]] <- data.frame(Sample=sid,n_ACKR1EC=nrow(ec),n_CD8T=nrow(td),
                                 med_CD8T=round(median(dc),1),med_Random=round(median(dr),1),
                                 p_value=signif(wilcox.test(dc,dr)$p.value,3))
}
dist_df <- do.call(rbind, dist_all); stat_df <- do.call(rbind, stat_rows)
write.csv(stat_df, file.path(out_dir,"07_Proximity_stats.csv"), row.names=FALSE)
dist_df$Group <- factor(dist_df$Group, levels=c("Random spot","CD8 T cells"))
gcol <- c("CD8 T cells"="#2C6FA6","Random spot"="grey72")

p_overall <- wilcox.test(rel~Group, data=dist_df)$p.value
ymax <- as.numeric(quantile(dist_df$rel,0.97))
med_pt <- dist_df %>% group_by(Sample,Group) %>% summarise(m=median(rel), .groups="drop")
p_v <- ggplot(dist_df, aes(Group, rel, fill=Group)) +
  geom_violin(scale="width", width=.8, linewidth=.3, color="grey35", alpha=.5, trim=TRUE) +
  geom_boxplot(width=.14, outlier.shape=NA, linewidth=.45, fatten=1.4) +
  geom_point(data=med_pt, aes(Group,m), inherit.aes=FALSE,
             position=position_jitter(width=.06,height=0), size=2.4, shape=21,
             fill="white", color="black", stroke=.7) +
  geom_hline(yintercept=1, linetype="dashed", color="grey55", linewidth=.4) +
  scale_fill_manual(values=gcol) + coord_cartesian(ylim=c(0,ymax)) +
  annotate("text", x=1.5, y=ymax*.95, label=sprintf("p = %s",
                                                    format(p_overall,scientific=TRUE,digits=3)), fontface="bold", size=4.3) +
  labs(x=NULL, y="Relative distance to nearest ACKR1+ EC\n(1 = random median)",
       title="CD8 T cells are closer to ACKR1+ EC",
       subtitle=sprintf("Pooled %d patients · white dots = per-patient medians · lower = closer",
                        length(unique(dist_df$Sample)))) +
  theme_nat(14) + theme(legend.position="none",
                        axis.text.x=element_text(face="bold", size=13))
ggsave(file.path(out_dir,"07b_Proximity_violin.pdf"), p_v, width=6.8, height=6)
ggsave(file.path(out_dir,"07b_Proximity_violin.png"), p_v, width=6.8, height=6, dpi=300)

p_ec1 <- ggplot(dist_df, aes(rel, color=Group)) + stat_ecdf(linewidth=1.1) +
  scale_color_manual(values=gcol, name=NULL) + coord_cartesian(xlim=c(0, ymax)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60", linewidth=.4) +
  labs(x="Relative distance to nearest ACKR1+ EC", y="Cumulative fraction",
       title="CD8 T shift toward ACKR1+ EC (pooled)") +
  theme_nat(13) + theme(legend.position=c(.8,.3))
ggsave(file.path(out_dir,"07c_Proximity_ECDF_pooled.pdf"), p_ec1, width=6.5, height=5.5)
ggsave(file.path(out_dir,"07c_Proximity_ECDF_pooled.png"), p_ec1, width=6.5, height=5.5, dpi=300)

p_ec2 <- ggplot(dist_df, aes(dist, color=Group)) + stat_ecdf(linewidth=1) +
  scale_color_manual(values=gcol, name=NULL) + facet_wrap(~Sample, scales="free_x", nrow=2) +
  labs(x="Distance to nearest ACKR1+ EC spot", y="Cumulative fraction",
       title="Distance distribution per patient") +
  theme_nat(12) + theme(legend.position="top", panel.grid=element_blank())
ggsave(file.path(out_dir,"07c_Proximity_ECDF_facet.pdf"), p_ec2, width=14, height=6.5)
ggsave(file.path(out_dir,"07c_Proximity_ECDF_facet.png"), p_ec2, width=14, height=6.5, dpi=300)

sl <- med_pt %>% tidyr::pivot_wider(names_from=Group, values_from=m)
sl$dir <- ifelse(sl$`CD8 T cells` < sl$`Random spot`, "Closer", "Farther")
p_sl <- ggplot(med_pt, aes(Group, m, group=Sample)) +
  geom_line(data=med_pt %>% left_join(sl[,c("Sample","dir")],by="Sample"),
            aes(color=dir), linewidth=.9) +
  geom_point(size=2.6, shape=21, fill="white", color="black", stroke=.7) +
  geom_text(data=sl, aes(x=2, y=`CD8 T cells`, label=Sample), inherit.aes=FALSE,
            hjust=-.25, size=3.4) +
  geom_hline(yintercept=1, linetype="dashed", color="grey60", linewidth=.4) +
  scale_color_manual(values=c("Closer"="#2C6FA6","Farther"="grey60"), name=NULL) +
  labs(x=NULL, y="Median relative distance to ACKR1+ EC",
       title="Per-patient median distance: CD8 T vs random",
       subtitle="Line sloping down = CD8 T closer to ACKR1+ EC") +
  theme_nat(13) + theme(axis.text.x=element_text(face="bold",size=12))
ggsave(file.path(out_dir,"07d_Proximity_slopegraph.pdf"), p_sl, width=6.5, height=6)
ggsave(file.path(out_dir,"07d_Proximity_slopegraph.png"), p_sl, width=6.5, height=6, dpi=300)

message("\n>>> 全部出版级图重画完成！统计："); print(stat_df, row.names=FALSE)
message(sprintf(">>> 合并 Wilcoxon p = %.2e", p_overall))

