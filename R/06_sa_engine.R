# 06_sa_engine.R: 模拟退火核心引擎
source("R/05_neighborhood.R")  # 依赖邻域函数

sa_optimize_H_fast <- function(m, n, verbose = TRUE) {
  # 输入检查：n 必须为偶数（foldover 设计）
  if (n %% 2 != 0) stop("❌ SA 引擎：n 必须为偶数")
  K <- n / 2  # 代表集大小
  allp <- factorial(m)  # 全排列总数（用于初始化 H）
  
  # 全局最优初始化
  best_global <- list(obj = -Inf, H = NULL, dmat = NULL, metrics = NULL)
  
  # 多重启策略（避免局部最优）
  for (r in seq_len(SA_PARAMS$restarts)) {
    if (verbose) cat(sprintf("\n[SA 重启 %d/%d] 初始化代表集 H...\n", r, SA_PARAMS$restarts))
    
    # ---------------------- 1. 初始化代表集 H ----------------------
    H <- if (allp <= 20000) {
      # 小 m：生成全排列后随机采样 K 个（保证不重复）
      perms_all <- permutations(m, m, 1:m)
      rows <- sample(nrow(perms_all), K)
      lapply(rows, function(i) as.integer(perms_all[i, ]))
    } else {
      # 大 m：直接随机生成 K 个排列（允许重复，效率优先）
      replicate(K, sample.int(m, m), simplify = FALSE)
    }
    
    # 初始化距离矩阵与目标值
    dmat <- compute_dmat_from_H_parallel(H, m)
    cur_obj_info <- objective_scalar_from_dmat_fast(dmat, m)
    cur_obj <- cur_obj_info$val
    cur_metrics <- cur_obj_info$metrics
    
    # 局部最优初始化
    best_local <- list(obj = cur_obj, H = H, dmat = dmat, metrics = cur_metrics)
    
    # ---------------------- 2. SA 迭代 ----------------------
    T <- SA_PARAMS$T0  # 初始温度
    no_improve <- 0L   # 无改进迭代次数（早停用）
    
    for (it in seq_len(SA_PARAMS$max_iter)) {
      # （1）生成候选解（小跳跃/大跳跃）
      nb <- if (runif(1) < SA_PARAMS$large_move_prob) {
        neighbor_replace_row(H, m)  # 大跳跃
      } else {
        neighbor_swap_row(H)        # 小跳跃
      }
      
      # （2）增量更新距离矩阵与 H
      upd <- update_dmat_single_row_fast(dmat, H, nb$idx, nb$newperm, m)
      newdmat <- upd$dmat; newH <- upd$H
      
      # （3）计算候选解的目标值
      new_obj_info <- objective_scalar_from_dmat_fast(newdmat, m)
      delta <- new_obj_info$val - cur_obj
      
      # （4）Metropolis 接受准则
      accept <- FALSE
      if (delta >= 0) {
        accept <- TRUE  # 更优解：接受
      } else if (runif(1) < exp(delta / T)) {
        accept <- TRUE  # 较差解：按概率接受
      }
      
      # （5）更新状态
      if (accept) {
        H <- newH; dmat <- newdmat; cur_obj <- new_obj_info$val; cur_metrics <- new_obj_info$metrics
        no_improve <- 0L  # 重置无改进计数
        # 更新局部最优
        if (cur_obj > best_local$obj) {
          best_local <- list(obj = cur_obj, H = H, dmat = dmat, metrics = cur_metrics)
        }
      } else {
        no_improve <- no_improve + 1L  # 累积无改进计数
      }
      
      # （6）降温与早停
      T <- T * SA_PARAMS$alpha  # 几何降温
      if (no_improve > SA_PARAMS$no_improve_stop || T < 1e-8) {
        if (verbose) cat(sprintf("  迭代终止：无改进 %d 次 / 温度 < 1e-8\n", no_improve))
        break
      }
    }
    
    # ---------------------- 3. 更新全局最优 ----------------------
    if (verbose) {
      cat(sprintf("  重启 %d 最优：minK = %d, ave_sqK_all = %.4f\n",
                  r, best_local$metrics$minK, best_local$metrics$ave_sqK_all))
    }
    if (best_local$obj > best_global$obj) {
      best_global <- best_local
    }
  }
  
  best_global  # 返回全局最优解
}

message("✅ 模拟退火核心引擎加载完成")