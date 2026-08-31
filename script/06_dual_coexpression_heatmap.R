# ==============================================================================
# 06_dual_coexpression_heatmap.R: 基于真实互作数据的跨物种共表达热图
# ==============================================================================

library(pheatmap)
library(dplyr)
library(RColorBrewer)

if(!dir.exists("results")) dir.create("results")

cat("=== 1. 读取真实互作边表与精炼注释... ===\n")

edges <- read.csv("results/cross_species_network_edges.csv", stringsAsFactors = FALSE)
host_annot <- read.csv("results/parsed_host_15_genes_refined.csv", stringsAsFactors = FALSE)
bact_annot <- read.csv("results/parsed_bact_15_genes_refined.csv", stringsAsFactors = FALSE)

cat("=== 2. 构建 Host x Bact 关联矩阵... ===\n")

# 获取有数据的宿主与细菌基因集合
unique_hosts <- unique(edges$Host_Gene)
unique_bacts <- unique(edges$Bact_Gene)

# 初始化空矩阵
cor_mat <- matrix(0, nrow = length(unique_hosts), ncol = length(unique_bacts))
rownames(cor_mat) <- unique_hosts
colnames(cor_mat) <- unique_bacts

# 填充真实 Correlation 数值
for(i in 1:nrow(edges)) {
  h <- edges$Host_Gene[i]
  b <- edges$Bact_Gene[i]
  r <- edges$Correlation[i]
  if(h %in% rownames(cor_mat) && b %in% colnames(cor_mat)) {
    cor_mat[h, b] <- r
  }
}

cat("=== 3. 映射美化标签 (Display Labels)... ===\n")

# 建立 ID 到学术美化名的映射字典
host_map <- setNames(host_annot$clean_display, host_annot$raw_gene_id)
bact_map <- setNames(bact_annot$clean_display, bact_annot$raw_gene_id)

# 替换矩阵行列名
rownames(cor_mat) <- ifelse(rownames(cor_mat) %in% names(host_map), host_map[rownames(cor_mat)], rownames(cor_mat))
colnames(cor_mat) <- ifelse(colnames(cor_mat) %in% names(bact_map), bact_map[colnames(cor_mat)], colnames(cor_mat))

cat("=== 4. 导出高清学术热图... ===\n")

# 设置经典红蓝渐变配色
color_palette <- colorRampPalette(c("#377EB8", "#FFFFFF", "#E41A1C"))(100)



# 导出 PNG
png("plots/Host_Pathogen_Coexpression_Heatmap.png", width = 1000, height = 800, res = 120)
pheatmap(
  cor_mat,
  color = color_palette,
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Host - Pathogen Cross-Species Co-expression Heatmap (Spearman r)",
  angle_col = 45
)
dev.off()

cat("\n🎉 热图绘制完成！已保存至:\n")
cat("  1. 矢量图: results/Host_Pathogen_Coexpression_Heatmap.pdf\n")
cat("  2. 预览图: results/Host_Pathogen_Coexpression_Heatmap.png\n")

# 清理内存
rm(edges, host_annot, bact_annot, cor_mat, host_map, bact_map)
gc(verbose = FALSE)