# ==============================================================================
# Environment Setup & Dependency Installation Script
# ==============================================================================
# 本脚本用于快速环境搭建：
# 1. 自动递归扫描当前文件夹及子文件夹中的所有 .R 脚本；
# 2. 正则提取所有 library() 和 require() 中的 R 包依赖；
# 3. 自动排除 R 内置基础包，并对比本地已安装包；
# 4. 优先并自动安装 `pak`，利用 `pak::pkg_install()` 极速并行安装缺失依赖。
# ==============================================================================

# 1. 检查并确保安装了 pak（优先使用 CRAN 镜像）
if (!requireNamespace("pak", quietly = TRUE)) {
  message("--> 未检测到 `pak`，正在安装以实现高速并行安装...")
  install.packages("pak", repos = "https://cloud.r-project.org")
}

# 2. 递归获取项目根目录下所有的 .R 脚本文件
r_files <- list.files(path = ".", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)

# 过滤掉本安装脚本自身，防止自我扫描
current_script <- suppressWarnings(normalizePath(sys.frames()[[1]]$ofile, winslash = "/"))
if (!is.null(current_script)) {
  r_files <- r_files[normalizePath(r_files, winslash = "/") != current_script]
}

message(sprintf("--> 在项目中共找到 %d 个 R 脚本，开始解析依赖包...", length(r_files)))

# 3. 正则表达式提取 library() 与 require() 中的包名
pkg_pattern <- "(?:library|require)\\s*\\(\\s*[\"']?([a-zA-Z0-9.]+)(?:[\"']?|.*)?\\)"
extracted_pkgs <- character(0)

for (file in r_files) {
  lines <- readLines(file, warn = FALSE)
  # 忽略被注释掉的行 (#)
  clean_lines <- lines[!grepl("^\\s*#", lines)]
  matches <- regmatches(clean_lines, gregexpr(pkg_pattern, clean_lines, perl = TRUE))
  for (m in unlist(matches)) {
    pkg <- gsub(pkg_pattern, "\\1", m, perl = TRUE)
    extracted_pkgs <- c(extracted_pkgs, pkg)
  }
}

# 提取唯一的包列表
needed_pkgs <- unique(extracted_pkgs)

# 排除 R 内置的标准基础包（Base/Recommended packages，无需额外安装）
base_pkgs <- c("base", "compiler", "datasets", "graphics", "grDevices", "grid", 
               "methods", "parallel", "splines", "stats", "stats4", "tcltk", 
               "tools", "utils")
needed_pkgs <- setdiff(needed_pkgs, base_pkgs)

cat("\n==================================================\n")
message("检测到的项目所需 R 包清单：")
print(needed_pkgs)
cat("==================================================\n\n")

# 4. 对比本地环境，筛选缺失的包
installed_pkgs <- rownames(installed.packages())
missing_pkgs <- setdiff(needed_pkgs, installed_pkgs)

# 5. 使用 pak 进行一键极速安装
if (length(missing_pkgs) > 0) {
  message(sprintf("--> 发现 %d 个缺失的依赖包，正在使用 pak 自动安装...", length(missing_pkgs)))
  # pak 会自动区分该包属于 CRAN 还是 Bioconductor 并解决依赖关系
  pak::pkg_install(missing_pkgs)
  message("\n [成功] 所有依赖包已安装完毕！即可开始运行分析。")
} else {
  message("\n [提示] 当前环境中所有必需包均已安装，无需重复安装！")
}