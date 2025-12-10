# 03_distance_matrix.R: 距离矩阵计算与增量更新
source("R/02_utils.R")  # 依赖工具函数

# ---------------------- 1. 并行计算 H 的距离矩阵（K × K） ----------------------
compute_dmat_from_H_parallel <- function(H, m) {
  K <- length(H)
  if (K <= 1) return(matrix(0L, K, K))
  
  # 生成所有 i<j 的索引对
  idx_pairs <- expand.grid(i = 1:(K-1), j = 2:K) %>%
    filter(i < j) %>%
    pmap(function(i, j) list(i = i, j = j))  # 转为列表便于并行
  
  # 并行计算每对的 Kendall 距离
  dvals <- future_map_dbl(
    .x = idx_pairs,
    .f = ~kendall_tau_fast(H[[.x$i]], H[[.x$j]]),
    .options = furrr_options(seed = TRUE)  # 并行种子（可复现）
  )
  
  # 填充对称距离矩阵
  dmat <- matrix(0L, nrow = K, ncol = K)
  for (k in seq_along(idx_pairs)) {
    i <- idx_pairs[[k]]$i; j <- idx_pairs[[k]]$j
    dmat[i, j] <- as.integer(dvals[k])
    dmat[j, i] <- as.integer(dvals[k])
  }
  dmat
}

# ---------------------- 2. 增量更新距离矩阵（替换 H 的某一行） ----------------------
update_dmat_single_row_fast <- function(dmat, H, idx, newperm, m) {
  K <- length(H)
  
  # 计算 newperm 与 H 中所有排列的距离
  newd <- map_dbl(seq_len(K), function(j) {
    if (j == idx) 0L else kendall_tau_fast(newperm, H[[j]])
  }) %>% as.integer()
  
  # 更新对称矩阵
  dmat_new <- dmat
  dmat_new[idx, ] <- newd
  dmat_new[, idx] <- newd
  
  # 更新 H
  H_new <- H
  H_new[[idx]] <- newperm
  
  list(H = H_new, dmat = dmat_new)
}

message("✅ 距离矩阵处理函数加载完成")