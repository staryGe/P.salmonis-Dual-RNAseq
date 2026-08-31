# ---------------------------------------------------------------
# 01_bact_deseq2.R: Piscirickettsia salmonis 细菌侧差异表达分析
# ---------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)


# 1. 读取细菌侧 Counts 矩阵
bact_file <- file.path("raw_data", "GSE254974_CGR02_counts_extra-intracellular.txt.gz")
bact_counts <- read.table(bact_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 2. 构建细菌侧的 Metadata (元数据)
# 对照组基准: ExtraCellular (体外培养)
# 实验组: Intracellular (胞内感染)
bact_coldata <- data.frame(
  row.names = colnames(bact_counts),
  group = factor(c(rep("ExtraCellular", 3), rep("Intracellular", 3)), 
                 levels = c("ExtraCellular", "Intracellular"))
)

cat("================ 细菌侧 Metadata ================\n")
print(bact_coldata)

# 3. 构建 DESeq2 对象
dds_bact <- DESeqDataSetFromMatrix(countData = bact_counts,
                                   colData = bact_coldata,
                                   design = ~ group)

# 4. 过滤低表达基因 (至少在 3 个样本中 count >= 10)
keep_bact <- rowSums(counts(dds_bact) >= 10) >= 3
dds_bact <- dds_bact[keep_bact, ]
cat("\n过滤低表达后剩余细菌基因数:", sum(keep_bact), "\n")

# 5. 运行 DESeq2 差异分析
cat("\n=== 正在运行细菌侧 DESeq2 差异分析... ===\n")
dds_bact <- DESeq(dds_bact)

# 6. 提取结果 (Intracellular vs ExtraCellular)
res_bact <- results(dds_bact, contrast = c("group", "Intracellular", "ExtraCellular"))
res_bact_df <- as.data.frame(res_bact)
res_bact_df <- res_bact_df[order(res_bact_df$padj), ]

# 7. 统计显著差异基因 (padj < 0.05 且 |log2FC| > 1)
bact_up <- sum(res_bact_df$padj < 0.05 & res_bact_df$log2FoldChange > 1, na.rm = TRUE)
bact_down <- sum(res_bact_df$padj < 0.05 & res_bact_df$log2FoldChange < -1, na.rm = TRUE)

cat("\n================ 细菌侧差异分析统计 ================\n")
cat("胞内显著上调基因数 (Up in Intracellular):", bact_up, "\n")
cat("胞内显著下调基因数 (Down in Intracellular):", bact_down, "\n")
cat("差异基因总数 (Total DEG):", bact_up + bact_down, "\n\n")

cat("最显著的前 10 个细菌差异基因:\n")
print(head(res_bact_df, 10))

# 8. 保存细菌侧差异分析结果到本地
if(!dir.exists("results")) dir.create("results")
write.csv(res_bact_df, file = "results/bact_DEG_results.csv")
cat("\n细菌侧差异表达结果已保存至: results/bact_DEG_results.csv\n")

#########################火山图绘制##########################

# ---------------------------------------------------------------
# 02_bact_volcano.R: 绘制 P. salmonis 细菌侧差异火山图
# ---------------------------------------------------------------

library(ggplot2)
library(dplyr)

# 1. 读取上一步保存的差异结果
res_bact_df <- read.csv("results/bact_DEG_results.csv", row.names = 1)

# 2. 添加显著性分类标签
res_bact_df <- res_bact_df %>%
  mutate(change = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up (Intracellular)",
    padj < 0.05 & log2FoldChange < -1 ~ "Down (Intracellular)",
    TRUE ~ "Not Significant"
  ))

# 3. 筛选 Top 10 上调基因进行标记
top10_up <- res_bact_df %>%
  filter(change == "Up (Intracellular)") %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(10)

# 4. 绘制火山图
p_volcano <- ggplot(res_bact_df, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Up (Intracellular)" = "#d95f02", 
                                "Down (Intracellular)" = "#7570b3", 
                                "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", size = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.6) +
  labs(title = "Piscirickettsia salmonis Transcriptional Response in Host Cells",
       subtitle = "Intracellular Infection vs. Axenic Growth (GSE254974)",
       x = "Log2 Fold Change (Intracellular / Axenic)",
       y = "-Log10 (Adjusted P-value)",
       color = "Expression State") +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# 5. 保存图件
if(!dir.exists("plots")) dir.create("plots")
ggsave("plots/bact_volcano_plot.png", plot = p_volcano, width = 8, height = 6, dpi = 300)

cat("火山图已成功生成并保存至: plots/bact_volcano_plot.png\n")

# ===============================================================
# 模块三：上调基因 Top 10 热图制作 (自动提取 STRING 蛋白描述)
# ===============================================================
cat("\n=== 正在生成上调基因 Top 10 热图 (含 STRING 蛋白自动注释)... ===\n")

# 1. 获取 Top 10 上调基因的 Gene ID 列表
top10_gene_ids <- rownames(top10_up)

# 2. 提取 DESeq2 标准化表达量 (VST 变换)
vsd_bact <- vst(dds_bact, blind = FALSE)
mat_top10 <- assay(vsd_bact)[top10_gene_ids, ]

# 3. 利用 STRINGdb 自动获取这 10 个基因的官方描述 / Symbol
cat("=== 正在通过 STRINGdb 获取 Top 10 基因的可读注释... ===\n")
string_db <- STRINGdb$new(
  version = "12.0",
  species = 1238,
  score_threshold = 400,
  input_directory = ""
)

# 映射这 Top 10 个基因
top10_df <- data.frame(gene = top10_gene_ids)
top10_mapped <- string_db$map(top10_df, "gene")
top10_info <- string_db$get_proteins() %>% 
  filter(protein_external_id %in% top10_mapped$STRING_id)

# 匹配并生成可读标签 "preferred_name (AWJ11_xxxx)"
# 如果没有匹配到通用名，则保留原始 Tag
readable_labels <- sapply(top10_gene_ids, function(gid) {
  sid <- top10_mapped$STRING_id[top10_mapped$gene == gid]
  p_name <- top10_info$preferred_name[top10_info$protein_external_id == sid]
  
  if (length(p_name) > 0 && !is.na(p_name) && p_name != "") {
    return(paste0(p_name, " (", gid, ")"))
  } else {
    return(gid) # 兜底逻辑：若无通用名则保留原始 ID
  }
})

# 4. 将可读标签赋予表达矩阵的行名
rownames(mat_top10) <- readable_labels

# 5. 构建顶部样本分组颜色条
annotation_col <- data.frame(
  Group = bact_coldata$group,
  row.names = rownames(bact_coldata)
)

ann_colors <- list(
  Group = c("ExtraCellular" = "#2b5c8f", "Intracellular" = "#d95f02")
)

# 6. 绘制热图并保存
color_palette <- colorRampPalette(c("#313695", "#ffffbf", "#a50026"))(100)

pheatmap(
  mat_top10,
  scale = "row",                     # Z-score 按行标准化
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  cluster_cols = FALSE,              # 保持样本按组顺序排列
  show_rownames = TRUE,             # 显示转换后的可读基因名
  show_colnames = TRUE,             # 显示样本名
  annotation_col = annotation_col,  # 样本分组顶部注释
  annotation_colors = ann_colors,   # 分组颜色匹配
  color = color_palette,
  main = "P. salmonis Top 10 Upregulated Genes (Annotated)",
  filename = "plots/bact_top10_up_heatmap.png",
  width = 8.5,
  height = 6.5,
  dpi = 300
)

cat("带蛋白注释的 Top 10 热图已更新保存至: plots/bact_top10_up_heatmap.png\n")