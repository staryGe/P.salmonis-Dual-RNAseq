# ---------------------------------------------------------------
# 检查 GSE254974 数据集结构与基因名称
# ---------------------------------------------------------------

# 指定子文件夹 raw_data 下的文件路径
data_file <- file.path("raw_data", "GSE254974_CGR02_counts_extra-intracellular.txt.gz")

# 1. 读取数据
bact_count <- read.table(data_file, 
                 header = TRUE, 
                 row.names = 1, 
                 sep = "\t", 
                 check.names = FALSE)

# 2. 查看矩阵维度与所有样本列名
cat("=== 矩阵维度 (基因数 x 样本数) ===\n")
print(dim(bact_count))

cat("\n=== 所有样本列名 (Sample IDs) ===\n")
print(colnames(bact_count))

# 3. 查看前 10 个行名 (Gene IDs) 确认是否包含宿主/细菌
cat("\n=== 前 10 个基因名预览 ===\n")
print(head(rownames(bact_count), 10))

# 指定子文件夹 raw_data 下的文件路径
data_file <- file.path("raw_data", "GSE254974_SHK_counts_infected-control.txt.gz")

# 1. 读取数据
host_count <- read.table(data_file, 
                   header = TRUE, 
                   row.names = 1, 
                   sep = "\t", 
                   check.names = FALSE)

# 2. 查看矩阵维度与所有样本列名
cat("=== 矩阵维度 (基因数 x 样本数) ===\n")
print(dim(host_count))

cat("\n=== 所有样本列名 (Sample IDs) ===\n")
print(colnames(host_count))

# 3. 查看前 10 个行名 (Gene IDs) 确认是否包含宿主/细菌
cat("\n=== 前 10 个基因名预览 ===\n")
print(head(rownames(host_count), 10))