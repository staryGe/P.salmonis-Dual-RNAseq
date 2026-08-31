# ---------------------------------------------------------------
# 04_host_deseq2.R: 大西洋鲑鱼 (SHK-1 巨噬细胞) 宿主侧差异表达分析 + 火山图 + Top10热图
# ---------------------------------------------------------------

# 0. 环境准备：检测并自动安装绘图所需依赖包
needed_packages <- c("DESeq2", "ggplot2", "dplyr", "pheatmap", "RColorBrewer", "STRINGdb")
for (pkg in needed_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("pak", quietly = TRUE)) {
      install.packages("pak", repos = "https://cran.rstudio.com/")
    }
    cat(sprintf("正在使用 pak 安装缺少包: %s ...\n", pkg))
    pak::pkg_install(pkg)
  }
}

library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(STRINGdb)
options(timeout = 600)

# ===============================================================
# 模块一：DESeq2 宿主差异表达分析
# ===============================================================

# 1. 读取宿主侧 Counts 矩阵
host_file <- file.path("raw_data", "GSE254974_SHK_counts_infected-control.txt.gz")
host_counts <- read.table(host_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 2. 构建宿主侧 Metadata (元数据)
host_coldata <- data.frame(
  row.names = colnames(host_counts),
  group = factor(c(rep("Infected", 3), rep("Control", 3)), 
                 levels = c("Control", "Infected")) # Control 作为对照基准
)

cat("================ 宿主侧 Metadata ================\n")
print(host_coldata)

# 3. 构建 DESeq2 对象
dds_host <- DESeqDataSetFromMatrix(countData = host_counts,
                                   colData = host_coldata,
                                   design = ~ group)

# 4. 过滤低表达基因 (至少在 3 个样本中 count >= 10)
keep_host <- rowSums(counts(dds_host) >= 10) >= 3
dds_host <- dds_host[keep_host, ]
cat("\n过滤低表达后剩余宿主基因数:", sum(keep_host), "\n")

# 5. 运行 DESeq2 差异分析
cat("\n=== 正在运行宿主侧 DESeq2 差异分析... ===\n")
dds_host <- DESeq(dds_host)

# 6. 提取结果 (Infected vs Control)
res_host <- results(dds_host, contrast = c("group", "Infected", "Control"))
res_host_df <- as.data.frame(res_host)
res_host_df <- res_host_df[order(res_host_df$padj), ]

# 7. 统计显著差异基因 (padj < 0.05 且 |log2FC| > 1)
host_up <- sum(res_host_df$padj < 0.05 & res_host_df$log2FoldChange > 1, na.rm = TRUE)
host_down <- sum(res_host_df$padj < 0.05 & res_host_df$log2FoldChange < -1, na.rm = TRUE)

cat("\n================ 宿主侧差异分析统计 ================\n")
cat("宿主显著上调基因数 (Up in Infected):", host_up, "\n")
cat("宿主显著下调基因数 (Down in Infected):", host_down, "\n")
cat("宿主差异基因总数 (Total DEG):", host_up + host_down, "\n\n")

cat("最显著的前 10 个宿主差异基因:\n")
print(head(res_host_df, 10))

# 8. 保存宿主侧差异分析结果
if(!dir.exists("results")) dir.create("results", recursive = TRUE)
write.csv(res_host_df, file = "results/host_DEG_results.csv")
cat("\n宿主侧差异表达结果已保存至: results/host_DEG_results.csv\n")


# ===============================================================
# 模块二：宿主侧火山图绘制 (Volcano Plot)
# ===============================================================
cat("\n=== 正在生成宿主侧差异火山图... ===\n")

# 1. 添加显著性分类标签
res_host_df <- res_host_df %>%
  mutate(change = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up (Infected)",
    padj < 0.05 & log2FoldChange < -1 ~ "Down (Infected)",
    TRUE ~ "Not Significant"
  ))

# 2. 筛选 Top 10 上调基因
top10_host_up <- res_host_df %>%
  filter(change == "Up (Infected)") %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(10)

# 3. 绘制火山图
p_host_volcano <- ggplot(res_host_df, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Up (Infected)" = "#e41a1c", 
                                "Down (Infected)" = "#377eb8", 
                                "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", size = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.6) +
  labs(title = "Salmo salar (Macrophage) Host Immune Response to P. salmonis",
       subtitle = "Infected vs. Uninfected Control (GSE254974)",
       x = "Log2 Fold Change (Infected / Control)",
       y = "-Log10 (Adjusted P-value)",
       color = "Expression State") +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# 4. 保存火山图
if(!dir.exists("plots")) dir.create("plots", recursive = TRUE)
ggsave("plots/host_volcano_plot.png", plot = p_host_volcano, width = 8, height = 6, dpi = 300)
cat("宿主火山图已成功保存至: plots/host_volcano_plot.png\n")


# ===============================================================
# 模块三：宿主上调基因 Top 10 热图制作 (修复分组颜色对齐版)
# ===============================================================
cat("\n=== 正在生成宿主上调基因 Top 10 热图... ===\n")

top10_host_gene_ids <- rownames(top10_host_up)

# 1. 提取 VST 表达矩阵 (确保列顺序与 host_coldata 一致)
vsd_host <- vst(dds_host, blind = FALSE)
mat_host_top10 <- assay(vsd_host)[top10_host_gene_ids, colnames(dds_host)]

# 2. 构建顶部样本分组注释表 (严格与表达矩阵列名一对应)
annotation_col_host <- data.frame(
  Group = host_coldata$group,
  row.names = rownames(host_coldata)
)

# 3. 自定义色彩：Infected 对应橙红 (#d95f02)，Control 对应深蓝 (#2b5c8f)
ann_colors_host <- list(
  Group = c("Infected" = "#d95f02", "Control" = "#2b5c8f")
)

# 4. 重新绘制并保存
color_palette <- colorRampPalette(c("#313695", "#ffffbf", "#a50026"))(100)

pheatmap(
  mat_host_top10,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  cluster_cols = FALSE,              # 保持表达矩阵原有列顺序
  show_rownames = TRUE,             # 显示基因名 (LOCxx)
  show_colnames = TRUE,             # 显示样本名 (CGR02_x / CTROL_x)
  annotation_col = annotation_col_host,
  annotation_colors = ann_colors_host,
  color = color_palette,
  main = "Salmo salar Top 10 Upregulated Genes in Response to Infection",
  filename = "plots/host_top10_up_heatmap.png",
  width = 8.5,
  height = 6.5,
  dpi = 300
)

cat("修正后的宿主 Top 10 热图已重新保存至: plots/host_top10_up_heatmap.png\n")

# ===============================================================
# 模块四：宿主下调基因 Top 10 热图制作 (Down-regulated Top 10)
# ===============================================================
cat("\n=== 正在生成宿主下调基因 Top 10 热图... ===\n")

# 1. 筛选 Top 10 下调基因 (log2FC < -1 且按 padj 升序、log2FC 降序/绝对值最大排序)
top10_host_down <- res_host_df %>%
  filter(change == "Down (Infected)") %>%
  arrange(padj, log2FoldChange) %>%
  head(10)

top10_host_down_ids <- rownames(top10_host_down)

# 2. 提取 VST 表达矩阵
mat_host_down_top10 <- assay(vsd_host)[top10_host_down_ids, colnames(dds_host)]

# 3. 绘制下调基因 Top 10 热图
pheatmap(
  mat_host_down_top10,
  scale = "row",                      # Z-score 按行标准化
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  cluster_cols = FALSE,               # 保持样本按组顺序排列 (Infected 在前, Control 在后)
  show_rownames = TRUE,              # 显示下调基因 ID
  show_colnames = TRUE,              # 显示样本名
  annotation_col = annotation_col_host,  # 复用上文的分组注释
  annotation_colors = ann_colors_host,   # 复用上文的颜色绑定 (Infected=橙, Control=蓝)
  color = color_palette,             # 蓝白红渐变
  main = "Salmo salar Top 10 Downregulated Genes in Response to Infection",
  filename = "plots/host_top10_down_heatmap.png",
  width = 8.5,
  height = 6.5,
  dpi = 300
)

cat("宿主 Top 10 下调基因热图已保存至: plots/host_top10_down_heatmap.png\n")
cat("\n================ 宿主侧全套热图生成完毕！ ================\n")