# 04_metrics.R: 评估指标与目标函数
source("R/03_distance_matrix.R")  # 依赖距离矩阵函数

# ---------------------- 1. 从距离矩阵计算评估指标 ----------------------
evaluate_from_dmat_fast <- function(dmat, m) {
  K <- nrow(dmat)
  p <- choose(m, 2)
  
  # 边缘情况：K=1（仅一个代表排列）
  if (K <= 1) {
    return(list(
      minK = ifelse(K == 1, p, NA_integer_),
      ave_K_all = 0,
      ave_sqK_all = 0,
      var_K_all = 0
    ))
  }
  
  # 1. 最小距离 minK（i<j 中 min(d_ij, p-d_ij) 的最小值）
  dvec <- dmat[upper.tri(dmat)]  # 上三角距离向量
  minK <- pmin(dvec, p - dvec) %>% min() %>% as.integer()
  
  # 2. 平均距离（无序对 → 有序对）
  sumd_unordered <- 2L * p * choose(K, 2) + K * p  # 闭合公式（无需遍历）
  total_unordered <- choose(2 * K, 2)
  ave_K_unordered <- sumd_unordered / total_unordered
  
  # 3. 平均平方距离（无序对 → 有序对）
  sumsq_pairs <- sum(2L * (dvec^2) + 2L * ((p - dvec)^2))
  sumsq_unordered <- sumsq_pairs + K * p^2
  ave_sqK_unordered <- sumsq_unordered / total_unordered
  
  # 4. 有序对指标（所有 n=2K 个排列的有序对）
  nfull <- 2L * K
  ave_K_all <- (2L * sumd_unordered) / (nfull^2)
  ave_sqK_all <- (2L * sumsq_unordered) / (nfull^2)
  var_K_all <- ave_sqK_all - (ave_K_all)^2  # 方差
  
  list(
    minK = minK,
    ave_K_all = as.numeric(ave_K_all),
    ave_sqK_all = as.numeric(ave_sqK_all),
    var_K_all = as.numeric(var_K_all)
  )
}

# ---------------------- 2. 目标函数（转指标为标量，用于 SA 接受准则） ----------------------
objective_scalar_from_dmat_fast <- function(dmat, m) {
  met <- evaluate_from_dmat_fast(dmat, m)
  # 主目标：max minK；次目标：min ave_sqK（用小系数破 ties）
  obj_val <- met$minK - SA_PARAMS$eps_secondary * met$ave_sqK_all
  list(val = as.numeric(obj_val), metrics = met)
}

message("✅ 评估指标与目标函数加载完成")