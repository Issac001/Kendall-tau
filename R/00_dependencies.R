# 00_dependencies.R: 依赖包安装与加载
required_pkgs <- c("gtools", "furrr", "future", "tidyverse", "pheatmap")

# 安装缺失包
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# 加载包
library(gtools)     # 生成全排列
library(furrr)      # 并行计算
library(future)     # 并行计划
library(tidyverse)  # 数据处理
library(pheatmap)   # 距离矩阵热图

# 配置并行环境（保守使用：CPU核心数-1）
set.seed(123)  # 并行随机种子（保证可复现）
workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)

message(sprintf("✅ 依赖包加载完成，并行 workers = %d", workers))