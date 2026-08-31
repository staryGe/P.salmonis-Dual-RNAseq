# ---------------------------------------------------------------
# 02_bact_GO_analysis.R: P. salmonis 细菌侧上调基因全维度功能/结构域富集分析
# ---------------------------------------------------------------

# 1. 环境准备：优先使用 pak 安装缺少包
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://cran.rstudio.com/")
}

needed_packages <- c("STRINGdb", "dplyr", "ggplot2")
for (pkg in needed_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    pak::pkg_install(pkg)
  }
}

library(STRINGdb)
library(dplyr)
library(ggplot2)

# 2. 读取 01 脚本输出的差异结果
deg_file <- "results/bact_DEG_results.csv"
if (!file.exists(deg_file)) {
  stop("错误: 未找到 results/bact_DEG_results.csv，请先运行 01_bact_deseq2.R！")
}

res_bact_df <- read.csv(deg_file, row.names = 1)

# 3. 提取显著上调基因
up_genes_df <- res_bact_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  mutate(gene = rownames(.))

cat("================ 细菌上调基因统计 ================\n")
cat("用于富集分析的胞内上调基因数量:", nrow(up_genes_df), "\n\n")

# 4. 初始化 STRINGdb (P. salmonis Taxonomy ID: 1238)
cat("=== 正在连接 STRING 数据库 (Taxonomy ID: 1238) ... ===\n")
string_db <- STRINGdb$new(
  version = "12.0",
  species = 1238,
  score_threshold = 400,
  input_directory = ""
)

# 5. 映射基因 ID
mapped_genes <- string_db$map(up_genes_df, "gene")
mapped_genes <- mapped_genes %>% filter(!is.na(STRING_id))

cat("成功映射到 STRING 数据库的基因数:", nrow(mapped_genes), "\n\n")

# ---------------------------------------------------------------
# 6. 运行全量富集分析 (直接提取全量 categories)
# ---------------------------------------------------------------
cat("=== 正在计算全维度功能与结构域富集... ===\n")

# 直接获取全量富集结果
enrichment <- string_db$get_enrichment(mapped_genes$STRING_id)

# 导出原始全量表格到 CSV
if (!dir.exists("results/enrichment")) dir.create("results/enrichment", recursive = TRUE)
write.csv(enrichment, "results/enrichment/bact_up_GO_enrichment.csv", row.names = FALSE)
cat("原始富集结果已保存至: results/enrichment/bact_up_GO_enrichment.csv\n")
cat("提取到的富集总条目数:", nrow(enrichment), "\n\n")

# ---------------------------------------------------------------
# 7. 按基因数/显著性选取 Top 15 进行绘图
# ---------------------------------------------------------------
# 优先选按 FDR 排序的前 15 个最显著条目
top_terms <- enrichment %>%
  arrange(fdr, desc(number_of_genes)) %>%
  head(15)

cat("================ 绘图选择的 Top 条目 ================\n")
print(top_terms %>% dplyr::select(term, description, fdr, number_of_genes))

# ---------------------------------------------------------------
# 8. 绘制富集气泡图 (Dotplot)
# ---------------------------------------------------------------
p_enrich <- ggplot(top_terms, aes(x = number_of_genes, y = reorder(description, number_of_genes))) +
  geom_point(aes(size = number_of_genes, color = fdr)) +
  scale_color_gradient(low = "#d95f02", high = "#2b5c8f", name = "FDR (p.adj)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "P. salmonis Intracellular Functional Enrichment",
    subtitle = "Top Enriched Domains & Pathways in Host Macrophage (STRING v12.0)",
    x = "Number of Genes",
    y = "Functional Category / Domain / Pathway"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.minor = element_blank()
  )

# 9. 保存图片
if (!dir.exists("plots")) dir.create("plots")
ggsave("plots/bact_GO_dotplot.png", plot = p_enrich, width = 10, height = 6.5, dpi = 300)

cat("\n================ 运行完成！ ================\n")
cat("气泡图已成功更新并保存至: plots/bact_GO_dotplot.png\n")