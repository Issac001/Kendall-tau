# 02_utils.R: 通用工具函数
source("R/01_config.R")  # 依赖全局参数

# ---------------------- 1. 快速 Kendall-tau 距离（归并排序算逆序） ----------------------
#merge sort kendall——tau
# number of invertions（O(n log n)）
count_inversions <- function(arr) {
  merge_sort <- function(x) {
    n <- length(x)
    if (n <= 1) return(list(arr = x, inv = 0))
    
    mid <- n %/% 2
    left <- merge_sort(x[1:mid])
    right <- merge_sort(x[(mid+1):n])
    
    merge_count <- function(left, right) {
      i <- j <- inv_count <- 0
      merged <- integer(length(left$arr) + length(right$arr))
      li <- length(left$arr)
      ri <- length(right$arr)
      
      while (i < li && j < ri) {
        if (left$arr[i + 1] <= right$arr[j + 1]) {
          merged[i + j + 1] <- left$arr[i + 1]
          i <- i + 1
        } else {
          merged[i + j + 1] <- right$arr[j + 1]
          j <- j + 1
          inv_count <- inv_count + (li - i) # all left elements left are invetions
        }
      }
      
      if (i < li) merged[(i + j + 1):(li + ri)] <- left$arr[(i + 1):li]
      if (j < ri) merged[(i + j + 1):(li + ri)] <- right$arr[(j + 1):ri]
      
      list(arr = merged, inv = left$inv + right$inv + inv_count)
    }
    
    merge_count(left, right)
  }
  
  merge_sort(arr)$inv
}

# Kendall tau dist
kendall_tau_fast <- function(itemsA, itemsB) {
  stopifnot(length(itemsA) == length(itemsB))
  
  n <- length(itemsA)
  
  #  
  indexesB <- integer(n)
  indexesB[itemsB] <- 1:n  # suppose itemsB is the permutaion of  1..n 
  
  #  transform itemsA into  itemsB 's postion index
  indexesA <- indexesB[itemsA]
  
  # 
  count_inversions(indexesA)
}

# ---------------------- 2. Foldover 映射（值反转：x → m+1-x） ----------------------
foldover_perm <- function(perm, m) {
  as.integer(m + 1 - perm)
}
# 若需反序映射，替换为：foldover_perm <- function(perm, m) rev(perm)

# ---------------------- 3. 计算 S_full（精确全谱值，仅 m ≤ 8） ----------------------
compute_S_full <- function(m) {
  if (m > 8) stop("❌ compute_S_full: m > 8 无法精确枚举全排列")
  
  # 生成所有排列
  perms_all <- permutations(m, m, 1:m)
  nfull <- nrow(perms_all)
  p <- choose(m, 2)
  
  # 构建 PWO 向量矩阵 Fmat（nfull × p）
  pairs <- combn(1:m, 2)  # 所有 (i,j) 对（i<j）
  Fmat <- matrix(0L, nrow = nfull, ncol = p)
  
  for (r in seq_len(nfull)) {
    perm <- perms_all[r, ]
    for (k in seq_len(ncol(pairs))) {
      i <- pairs[1, k]; j <- pairs[2, k]
      posi <- which(perm == i); posj <- which(perm == j)
      Fmat[r, k] <- ifelse(posi < posj, 1L, -1L)
    }
  }
  
  # 计算 S_full
  Sfull <- sum(((Fmat %*% t(Fmat))^2)) / (nfull^4)
  list(p = p, n_full = nfull, S_full = Sfull)
}

message("✅ 通用工具函数加载完成")