# 05_neighborhood.R: 邻域操作（SA 候选解生成）
source("R/02_utils.R")  # 依赖工具函数

# ---------------------- 1. 小跳跃：交换某一行的两个位置 ----------------------
neighbor_swap_row <- function(H) {
  K <- length(H)
  idx <- sample(K, 1)  # 随机选一行
  perm <- H[[idx]]
  pos <- sample.int(length(perm), 2)  # 随机选两个位置
  
  # 交换位置
  newperm <- perm
  newperm[pos] <- rev(newperm[pos])
  
  list(idx = idx, newperm = newperm)
}

# ---------------------- 2. 大跳跃：随机替换某一行（从全排列采样） ----------------------
neighbor_replace_row <- function(H, m) {
  K <- length(H)
  idx <- sample(K, 1)  # 随机选一行
  newperm <- sample.int(m, m)  # 随机生成一个排列
  
  list(idx = idx, newperm = newperm)
}

message("✅ 邻域操作函数加载完成")