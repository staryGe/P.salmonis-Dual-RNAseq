# ==============================================================================
# 05_calc_cross_species_edges.R: 细菌 Locus Tag 提取与全量跨物种共表达计算
# ==============================================================================

library(dplyr)
library(stringr)

if(!dir.exists("results")) dir.create("results")

cat("=== Step 1. 读取 15 个宿主与 15 个细菌精炼注释表... ===\n")

host_annot <- read.csv("results/parsed_host_15_genes_refined.csv", stringsAsFactors = FALSE)
bact_annot <- read.csv("results/parsed_bact_15_genes_refined.csv", stringsAsFactors = FALSE)

cat("=== Step 2. 精准提取 clean_display 中的 AWJ11_ Tag 并对接 bact_count... ===\n")

bact_all_rownames <- rownames(bact_count)

# 智能提取与映射函数
bact_annot <- bact_annot %>%
  mutate(
    # 提取 clean_display 中括号内的 AWJ11_XXXXX 编号
    extracted_tag = stringr::str_extract(clean_display, "AWJ11_\\d+"),
    
    # 确定最终矩阵匹配 ID 优先逻辑：
    # 1. raw_gene_id 本身就是 AWJ11_/Symbol 且存在于矩阵中
    # 2. 从 clean_display 提取出的 AWJ11_ Tag 存在于矩阵中
    matrix_row_id = case_when(
      raw_gene_id %in% bact_all_rownames ~ raw_gene_id,
      !is.na(extracted_tag) & extracted_tag %in% bact_all_rownames ~ extracted_tag,
      TRUE ~ NA_character_
    )
  )

cat("细菌 15 个基因的矩阵 ID 提取与映射结果预览：\n")
print(bact_annot %>% dplyr::select(raw_gene_id, clean_display, matrix_row_id))

# 筛选出有效匹配的 ID
valid_bact_df <- bact_annot %>% filter(!is.na(matrix_row_id))
valid_host_ids <- intersect(host_annot$raw_gene_id, rownames(host_count))

cat(sprintf("\n✅ 宿主基因匹配成功: %d / 15 个\n", length(valid_host_ids)))
cat(sprintf("✅ 细菌基因匹配成功: %d / 15 个\n\n", nrow(valid_bact_df)))

cat("=== Step 3. 提取表达量并做 Log2-CPM 标准化... ===\n")

# 1. 宿主 Log2-CPM
host_sub <- host_count[valid_host_ids, , drop = FALSE]
host_cpm <- log2(t(t(host_sub) / colSums(host_count)) * 1e6 + 1)

# 2. 细菌 Log2-CPM
bact_sub <- bact_count[valid_bact_df$matrix_row_id, , drop = FALSE]
bact_cpm <- log2(t(t(bact_sub) / colSums(bact_count)) * 1e6 + 1)

cat("=== Step 4. 计算 15x15 跨物种 Spearman 相关性... ===\n")

host_t <- t(host_cpm)
bact_t <- t(bact_cpm)

# 建立了映射向量： matrix_row_id -> raw_gene_id / clean_display
bact_id_map <- setNames(valid_bact_df$raw_gene_id, valid_bact_df$matrix_row_id)

edges_df <- data.frame()

for (h in colnames(host_t)) {
  for (b_mat_id in colnames(bact_t)) {
    test_res <- cor.test(host_t[, h], bact_t[, b_mat_id], method = "spearman", exact = FALSE)
    r_val <- as.numeric(test_res$estimate)
    p_val <- as.numeric(test_res$p.value)
    
    # 还原为原始细菌 ID
    b_raw_id <- bact_id_map[b_mat_id]
    
    # 保留显著调控对 (|r| >= 0.7 且 p < 0.05)
    if (!is.na(p_val) && abs(r_val) >= 0.7 && p_val < 0.05) {
      edges_df <- rbind(edges_df, data.frame(
        Host_Gene    = h,
        Bact_Gene    = b_raw_id,
        Bact_Mat_ID  = b_mat_id,
        Correlation  = r_val,
        P_value      = p_val,
        Interaction  = ifelse(r_val > 0, "Positive", "Negative"),
        stringsAsFactors = FALSE
      ))
    }
  }
}

write.csv(edges_df, "results/cross_species_network_edges.csv", row.names = FALSE)

cat(sprintf("\n🎉 跨物种共表达计算完美完成！\n"))
cat(sprintf("📊 在 |r| >= 0.7, p < 0.05 阈值下，筛选出 %d 条显著跨物种互作边。\n", nrow(edges_df)))
cat("📄 结果已保存至: results/cross_species_network_edges.csv\n\n")

print(head(edges_df))

# 内存安全清理
rm(host_annot, bact_annot, bact_all_rownames, valid_bact_df, valid_host_ids, 
   host_sub, bact_sub, host_cpm, bact_cpm, host_t, bact_t, bact_id_map, 
   test_res, r_val, p_val, h, b_mat_id, b_raw_id)
gc(verbose = FALSE)