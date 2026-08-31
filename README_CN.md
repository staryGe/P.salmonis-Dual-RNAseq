# *P. salmonis* & Host Dual RNA-seq & Cross-Species Co-expression Analysis

![R](https://img.shields.io/badge/R-4.6%2B-blue)
![DESeq2](https://img.shields.io/badge/Bioconductor-DESeq2-green)
![clusterProfiler](https://img.shields.io/badge/Bioconductor-clusterProfiler-orange)
![pak](https://img.shields.io/badge/R_pkg-pak-purple)

## 项目简介

本项目基于发表文献 *"Dual transcriptional analysis provides insights into the replicative niche of P. salmonis and the host response during infection"* 的公共双转录组数据，在大西洋鲑（*Salmo salar*，细胞系 SHK-1）与宿主胞内寄生菌鲑鱼立克次体（*Piscirickettsia salmonis*）的相互作用系统模型中，开展了系统的 **Dual RNA-seq 表达谱分析与跨物种共表达网络（Cross-Species Co-expression Analysis）** 练习。

通过同步解析感染过程中宿主与病原体的转录组动态变化，结合 Spearman 相关性分析与跨物种基因表达热图（Co-expression Heatmap），深度还原病原菌在宿主酸性液泡微环境（Acidified Vacuole Niche）中的适应机制与宿主防御响应。



## 核心生物学发现

1. **宿主溶酶体与液泡酸化响应（Host Lysosome & Acidification Response）**：
   - 宿主细胞大幅激活了**溶酶体发生（Lysosomal biogenesis）\**与\**蛋白水解**相关通路。
   - 跨物种共表达分析显示，驱动液泡酸化的核心质子泵亚基 **`ATP6V1E1` (LOC106609471)**、内体/溶酶体离子通道 **`CLCN5` (LOC106578043)** 以及溶酶体硫醇还原酶 **`GILT` (LOC106568950)** 均呈现显著的共表达调控，证实细菌处于酸化的晚期内体/溶酶体样结构中。
2. **病原菌胞内适应与关键毒力表达（Bacterial Intracellular Adaptation）**：
   - *P. salmonis* 配合宿主微环境显著上调了 **Dot/Icm IVB 型分泌系统（T4SS）**、**应激响应系统**（如传感激酶 **`ntrY` / AWJ11_01955**）以及**铁获取系统（Iron-acquisition systems）**。
3. **跨物种“铁争夺”与协同博弈（Host-Pathogen Cross-Species Interaction）**：
   - 细菌铁吸收关键基因 **`FedB` (AWJ11_07035)** 与宿主免疫受体 **`TNFRSF9`** 及氨基酸/离子转运蛋白 **`LAT4`** 呈现极高的负相关（$r \le -0.83$），从转录组水平印证了“**宿主铁可用性直接决定胞内细菌复制**”以及跨物种代谢博弈的分子假说。





---

## 分析成果展示

| 分析模块                   |                           结果展示                           | 说明                                                         |
| :------------------------- | :----------------------------------------------------------: | :----------------------------------------------------------- |
| **病原菌差异基因筛选**     |     <img src="plots/bact_volcano_plot.png" width="300"/>     | 火山图展现 *P. salmonis* 在宿主胞内与纯培养（Axenic）状态下的差异表达基因分布 |
| **病原菌 Top10 上调基因**  |   <img src="plots/bact_top10_up_heatmap.png" width="300"/>   | 热图展示细菌在胞内高度特异性激活的 Top 10 关键基因（含 Dot/Icm IVB 型分泌系统等） |
| **病原菌功能与结构域富集** |      <img src="plots/bact_GO_dotplot.png" width="300"/>      | 气泡图揭示细菌胞内上调基因显著富集于环境应激、铁获取与毒力相关蛋白结构域 |
| **宿主免疫响应火山图**     |     <img src="plots/host_volcano_plot.png" width="300"/>     | 火山图展现大西洋鲑（SHK-1 巨噬细胞）受感染后的全局转录组重塑 |
| **宿主 Top10 上调基因**    |   <img src="plots/host_top10_up_heatmap.png" width="300"/>   | 感染组中显著上调的 Top 10 宿主免疫与细胞应答基因表达热图     |
| **宿主 Top10 下调基因**    |  <img src="plots/host_top10_down_heatmap.png" width="300"/>  | 感染后表达显著受到抑制或关停的 Top 10 宿主基因表达热图       |
| **宿主 GO 多维度富集分析** |   <img src="plots/host_GO_split_dotplot.png" width="300"/>   | 拆分展示宿主在 BP（生物过程）、CC（细胞组件/溶酶体及囊泡）与 MF（分子功能）维度的功能富集 |
| **跨物种共表达热图**       | <img src="plots/Host_Pathogen_Coexpression_Heatmap.png" width="350"/> | **核心结果**：Spearman 相关性热图，捕获宿主液泡酸化/溶酶体基因（如 `ATP6V1E1`）与细菌毒力/铁吸收基因（如 `FedB`, `ntrY`）的跨物种协同调控网络 |



# 目录结构



```markdown
📁 BioPratice09_Dual_RNA/
├── 📁 plots/
│   ├── 📄 bact_GO_dotplot.png (0.11 MB)
│   ├── 📄 bact_top10_up_heatmap.png (<0.1 MB)
│   ├── 📄 bact_volcano_plot.png (0.50 MB)
│   ├── 📄 host_GO_split_dotplot.png (0.23 MB)
│   ├── 📄 Host_Pathogen_Coexpression_Heatmap.png (<0.1 MB)
│   ├── 📄 host_top10_down_heatmap.png (<0.1 MB)
│   ├── 📄 host_top10_up_heatmap.png (<0.1 MB)
│   └── 📄 host_volcano_plot.png (0.53 MB)
├── 📁 raw_data/
│   ├── 📄 GSE254974_CGR02_counts_extra-intracellular.txt.gz (<0.1 MB)
│   └── 📄 GSE254974_SHK_counts_infected-control.txt.gz (0.52 MB)
├── 📁 results/
│   ├── 📁 enrichment/
│   │   ├── 📄 bact_up_GO_enrichment.csv (<0.1 MB)
│   │   ├── 📄 host_up_gene_list.txt (<0.1 MB)
│   │   └── 📄 host_up_STRING_enrichment.csv (6.43 MB)
│   ├── 📄 bact_DEG_results.csv (0.34 MB)
│   ├── 📄 cross_species_network_edges.csv (<0.1 MB)
│   ├── 📄 heatmap_bact_genes_annotated.csv (<0.1 MB)
│   ├── 📄 heatmap_host_genes_annotated.csv (<0.1 MB)
│   ├── 📄 heatmap_host_genes_ortholog_annotated.csv (<0.1 MB)
│   ├── 📄 host_DEG_results.csv (3.41 MB)
│   ├── 📄 master_host_annotation.csv (6.25 MB)
│   ├── 📄 parsed_bact_15_genes_refined.csv (<0.1 MB)
│   └── 📄 parsed_host_15_genes_refined.csv (<0.1 MB)
├── 📁 script/
│   ├── 📄 00_1_data_check.R (<0.1 MB)
│   ├── 📄 01_bact_deseq2.R (<0.1 MB)
│   ├── 📄 02_bact_GO_analysis.R (<0.1 MB)
│   ├── 📄 03_host_deseq2.R (<0.1 MB)
│   ├── 📄 04_host_enrichment.R (<0.1 MB)
│   ├── 📄 05_calc_cross_species_edges.R (<0.1 MB)
│   └── 📄 06_dual_coexpression_heatmap.R (<0.1 MB)
├── 📄 BioPratice09_Dual_RNA.Rproj (<0.1 MB)
└── 📄 README_CN.md (<0.1 MB)
└── 📄 README.md (<0.1 MB)
```



# 环境依赖与快速复现



在项目根目录下，直接在 RStudio 控制台中运行： 

``` 
R source("run_pipeline.R")
```





# 调试与技术避坑记录



---

### 1. 跨物种基因 ID 规范统一与命名冲突
* **遇到问题**：宿主基因通常采用 NCBI Symbol 或 LOC ID（如 `LOC106609471`），而病原菌采用特定基因组的 Locus Tag（如 `AWJ11_07035`）。在构建跨物种相关性矩阵或绘图时，容易因为列名冲突或缺少注释导致热图标签不可读。
* **排查与解决**：
  * 在脚本 `05_calc_cross_species_edges.R` 中引入自动化命名空间前缀（`Host_` / `Pathogen_`）；
  * 手动构建独立的注释映射表（Metadata），在绘制共表达热图前将细菌 Locus Tag 动态关联至其已知基因名（如 `FedB`, `ntrY`），极大提升了成果图的可读性。

---

### 2. 复杂跨物种热图的标签重叠与聚类树枝剪裁
* **遇到问题**：使用 `pheatmap` 绘制跨物种共表达热图时，横纵轴基因数目较多，字体容易发生挤压和文本重叠，且缺少明确的数值标示。
* **排查与解决**：
  * 在 `06_dual_coexpression_heatmap.R` 中显式设置 `cellwidth` 与 `cellheight` 固定格子尺寸；
  * 使用 `display_numbers = TRUE` 将 Spearman $r$ 值保留两位小数打在热图网格中，并配置 `number_color`；
  * 使用自定义配色方案（蓝-白-红，中心点设为 0），更直观地对比正负相关强度。

---

### 3. 多脚本依赖冲突与环境复现
* **遇到问题**：复现生信项目时，使用者常因缺少 `DESeq2`、`clusterProfiler` 或 `pheatmap` 等依赖（分别分布于 CRAN 和 Bioconductor）而频繁报错中断。
* **排查与解决**：
  * 编写了 `install_dependencies.R` 脚本，引入 **`pak` 极速包管理器**，实现了全自动的依赖扫描、区分源（CRAN/Bioc）与并行安装；
  * 编写了总控脚本 `run_pipeline.R`，加入 `tryCatch()` 异常处理与步骤计时，实现了“一键无人值守”复现。