# ==============================================================================
# Dual RNA-seq 一键全流程主控运行脚本 (run_pipeline.R)
# ==============================================================================
# 本脚本按顺序依次执行依赖安装、数据校验、病原菌分析、宿主分析及跨物种共表达分析
# ==============================================================================

# 定义脚本严格按顺序执行的清单
pipeline_scripts <- c(
  "install_dependencies.R",
  "00_1_data_check.R",
  "01_bact_deseq2.R",
  "02_bact_GO_analysis.R",
  "03_host_deseq2.R",
  "04_host_enrichment.R",
  "05_calc_cross_species_edges.R",
  "06_dual_coexpression_heatmap.R"
)

cat("\n==================================================================\n")
message(" 启动 Dual RNA-seq 全流程主控分析 (Master Pipeline)")
cat("==================================================================\n\n")

start_time_total <- Sys.time()

for (i in seq_along(pipeline_scripts)) {
  script <- pipeline_scripts[i]
  
  cat("\n------------------------------------------------------------------\n")
  message(sprintf("[%d/%d] 正在启动脚本: %s", i, length(pipeline_scripts), script))
  cat("------------------------------------------------------------------\n")
  
  if (!file.exists(script)) {
    stop(sprintf("[错误] 在当前目录下找不到必填脚本 '%s'！请检查文件名或路径。", script))
  }
  
  step_start_time <- Sys.time()
  
  # 执行子脚本并捕获潜在异常
  tryCatch({
    source(script, local = FALSE)
    step_duration <- round(difftime(Sys.time(), step_start_time, units = "secs"), 2)
    message(sprintf("[完成] %s 顺利运行结束，耗时 %s 秒。\n", script, step_duration))
  }, error = function(e) {
    cat("\n==================================================================\n")
    message(sprintf("[中断] 流程在第 %d 步 (%s) 发生错误！", i, script))
    message("错误详细信息如下:")
    print(e)
    cat("==================================================================\n\n")
    stop(sprintf("由于 '%s' 运行报错，主控流程已安全终止。", script))
  })
}

total_duration <- round(difftime(Sys.time(), start_time_total, units = "mins"), 2)

cat("\n==================================================================\n")
message(sprintf("[恭喜] 全部 %d 个脚本成功运行完毕！总计耗时 %s 分钟。", length(pipeline_scripts), total_duration))
cat("==================================================================\n\n")