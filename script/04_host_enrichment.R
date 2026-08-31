# ---------------------------------------------------------------
# 05_host_enrichment.R: 大西洋鲑鱼 (SHK-1 巨噬细胞) 宿主侧 GO/KEGG 功能富集分析
# ---------------------------------------------------------------

# 0. 环境准备：检测并自动安装依赖包
needed_packages <- c("ggplot2", "dplyr", "STRINGdb")
for (pkg in needed_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("pak", quietly = TRUE)) {
      install.packages("pak", repos = "https://cran.rstudio.com/")
    }
    cat(sprintf("正在使用 pak 安装缺少包: %s ...\n", pkg))
    pak::pkg_install(pkg)
  }
}

library(ggplot2)
library(dplyr)
library(STRINGdb)


# 1. 设置超时时间与路径准备
options(timeout = 600)
if(!dir.exists("results/enrichment")) dir.create("results/enrichment", recursive = TRUE)
if(!dir.exists("plots")) dir.create("plots", recursive = TRUE)

cat("================ 正在读取宿主差异分析结果... ================\n")
host_deg_file <- "results/host_DEG_results.csv"

if (!file.exists(host_deg_file)) {
  stop("未找到宿主差异表 results/host_DEG_results.csv，请先运行 04_host_deseq2.R！")
}

res_host_df <- read.csv(host_deg_file, row.names = 1)

# 2. 提取宿主显著上调基因 (Up in Infected)
host_up_df <- res_host_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  arrange(padj)

host_up_genes <- rownames(host_up_df)
cat("宿主显著上调基因总数:", length(host_up_genes), "\n")

# 导出上调基因列表文件 (供备用)
write.table(host_up_genes, file = "results/enrichment/host_up_gene_list.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# 3. 初始化 STRINGdb (Salmo salar, TaxID: 8030)
cat("\n=== 正在连接 STRING 数据库获取大西洋鲑鱼功能富集分析结果... ===\n")

string_db_host <- STRINGdb$new(
  version = "12.0",
  species = 8030,           # 大西洋鲑鱼 NCBI Taxonomy ID
  score_threshold = 400,
  input_directory = ""
)

# 4. 映射基因 ID 并计算富集分析
mapped_host <- string_db_host$map(
  my_data_frame = data.frame(gene = host_up_genes),
  my_data_frame_id_col_names = "gene",
  removeUnmapped = TRUE
)

cat("成功映射到 STRING 数据库的宿主基因数:", nrow(mapped_host), "\n")

# 获取 GO / KEGG 富集条目
enrichment_host <- string_db_host$get_enrichment(mapped_host$STRING_id)

# 保存完整富集结果表
write.csv(enrichment_host, file = "results/enrichment/host_up_STRING_enrichment.csv", row.names = FALSE)
cat("宿主全量富集分析结果已保存至: results/enrichment/host_up_STRING_enrichment.csv\n")

# 5. 筛选并清洗数据：严格拆分 BP, CC, MF 并过滤 BP 顶层宽泛词条
cat("\n=== 正在按 BP, CC, MF 拆分 GO 富集条目并过滤宽泛大类... ===\n")

# 定义需要剔除的 GO-BP 顶层/超宽泛父节点列表
bp_parent_terms <- c(
  "Biological regulation", 
  "Cellular process", 
  "Regulation of biological process", 
  "Negative regulation of biological process", 
  "Regulation of cellular process",
  "Positive regulation of biological process",
  "Positive regulation of cellular process",
  "Negative regulation of cellular process",
  "Process", 
  "Biological process",
  "Cellular process regulation"
)

plot_data_go <- enrichment_host %>%
  # 1. 剔除 FDR 为空或不显著的行
  filter(!is.na(fdr) & as.numeric(fdr) < 0.05) %>%
  
  # 2. 映射与识别 GO 三大类 (Process = BP, Component = CC, Function = MF)
  mutate(
    GO_Class = case_when(
      category %in% c("Process", "Biological Process (GO)") | grepl("Process", category) ~ "Biological Process (BP)",
      category %in% c("Component", "Cellular Component (GO)") | grepl("Component", category) ~ "Cellular Component (CC)",
      category %in% c("Function", "Molecular Function (GO)") | grepl("Function", category) ~ "Molecular Function (MF)",
      TRUE ~ NA_character_
    )
  ) %>%
  
  # 3. 只保留能明确归类为 BP, CC, MF 的 GO 条目
  filter(!is.na(GO_Class)) %>%
  
  # 4. 精细过滤：剔除全局背景词条 以及 BP 中的顶层父节点 (基因数过大或属于通用词)
  filter(!description %in% c("Animal", "Whole body", "Tissues, cell types and enzyme sources")) %>%
  filter(!(GO_Class == "Biological Process (BP)" & (description %in% bp_parent_terms | as.numeric(number_of_genes) > 2500))) %>%
  
  # 5. 在每个类别 (BP/CC/MF) 内部按 FDR 升序，各取 Top 5 最显著的精细条目
  group_by(GO_Class) %>%
  arrange(as.numeric(fdr)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  
  # 6. 构造绘图所需变量
  mutate(
    term_desc = factor(description, levels = rev(unique(description))),
    fdr_val = as.numeric(fdr),
    gene_count = as.numeric(number_of_genes),
    neg_log10_fdr = -log10(fdr_val)
  )

cat("=== 精细化分类后的 Top 15 (BP/CC/MF 各 Top 5) 预览 ===\n")
print(plot_data_go %>% dplyr::select(GO_Class, term_desc, gene_count, fdr_val, neg_log10_fdr))

# 6. 绘制按 GO 三大类分面板 (Facet) 展示的气泡图
p_host_go_split <- ggplot(plot_data_go, aes(x = neg_log10_fdr, y = term_desc)) +
  geom_point(aes(size = gene_count, color = neg_log10_fdr)) +
  scale_color_gradient(low = "#377EB8", high = "#E41A1C", name = expression(-log[10](FDR))) +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  # 按 GO_Class 分面板，保持各自纵轴独立
  facet_grid(GO_Class ~ ., scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Salmo salar Macrophage GO Enrichment Analysis",
    subtitle = "Categorized by Biological Process (BP), Cellular Component (CC), and Molecular Function (MF)",
    x = expression(-log[10]("False Discovery Rate")),
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    strip.background = element_rect(fill = "#e0e0e0", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(color = "black", size = 10, face = "bold"),
    axis.text.x = element_text(color = "black", size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# 7. 保存分类 GO 富集图
ggsave("plots/host_GO_split_dotplot.png", plot = p_host_go_split, width = 10, height = 8, dpi = 300)


cat("\n严格区分 BP/CC/MF 的规范 GO 富集气泡图已保存至: plots/host_GO_split_dotplot.png\n")