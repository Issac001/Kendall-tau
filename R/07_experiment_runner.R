# 07_experiment_runner.R: 实验运行与结果保存
source("R/06_sa_engine.R")  # 依赖 SA 引擎

run_experiments <- function() {
  # 修复性检查：确保 results.csv 是文件，而不是目录
  results_csv_path <- OUTPUT_DIRS$results_csv
  if (dir.exists(results_csv_path)) {
    warning(sprintf("⚠️ 检测到 %s 是目录，已自动删除并重新创建为空文件。", results_csv_path))
    unlink(results_csv_path, recursive = TRUE)
  }
  if (!file.exists(results_csv_path)) {
    file.create(results_csv_path)
  }
  
  # 初始化结果列表
  results_list <- list()
  
  # 遍历所有 m（因子数量）
  for (m in m_range) {
    n <- total_n_from_m(m)  # 总设计行数（n = m，偶数）
    if (is.na(n)) {
      message(sprintf("⚠️  跳过 m=%d（需为偶数）", m))
      next
    }
    K <- n / 2  # 代表集大小
    
    cat("\n=======================================\n")
    cat(sprintf("开始实验：m=%d（因子数），n=%d（总设计行数），K=%d（代表集大小）\n", m, n, K))
    
    # ---------------------- 1. 运行 SA 优化 ----------------------
    best_design <- sa_optimize_H_fast(m = m, n = n, verbose = TRUE)
    met <- best_design$metrics  # 最优设计的指标
    
    # ---------------------- 2. 计算 S_full 与下界（仅 m ≤ 8） ----------------------
    # Sfull <- NA; ave_sqK_lower <- NA; minK_lower <- NA
    # if (m <= 8) {
    #   message("  计算 S_full 与理论下界...")
    #   sf <- compute_S_full(m)
    #   Sfull <- sf$S_full
    #   p <- sf$p
    #   mu <- met$ave_K_all
    #   
    #   # 平均平方距离下界（B）
    #   B <- (n^2 * Sfull - p^2 + 4 * p * mu) / 4
    #   ave_sqK_lower <- B
    #   
    #   # minK 下界（基于方差）
    #   mu2_lower <- B
    #   var_lower <- max(0, mu2_lower - mu^2)
    #   minK_lower <- max(0, mu - n * sqrt(var_lower))
    # }
    
    # ---------------------- 3. 保存结果 ----------------------
    # （1）结果行（用于 CSV）
    result_row <- tibble(
      m = m, n = n, K = K,
      minK = met$minK,
      ave_K_all = met$ave_K_all,
      ave_sqK_all = met$ave_sqK_all,
      var_K_all = met$var_K_all,
     # # S_full = Sfull,
     #  ave_sqK_lower = ave_sqK_lower,
     #  minK_lower = minK_lower,
       timestamp = Sys.time()
    )
    results_list[[length(results_list) + 1]] <- result_row
    
    # （2）汇总 CSV（增量写入）
    write_csv(
      x = result_row,
      file = OUTPUT_DIRS$results_csv,
      append = file.exists(OUTPUT_DIRS$results_csv)
    )
    
    # （3）最佳设计 RData（H_best + full_D）
    H_best <- best_design$H
    D_best <- do.call(rbind, H_best)  # 代表集矩阵
    D_fold <- t(apply(D_best, 1, foldover_perm, m = m))  # 折叠集矩阵
    full_D <- rbind(D_best, D_fold)  # 完整设计矩阵
    
    save(
      H_best, full_D, met,
      file = file.path(OUTPUT_DIRS$best_designs, sprintf("best_design_m%d_n%d.RData", m, n))
    )
    
    # （4）距离矩阵热图（可视化完整设计的距离）
    message("  生成距离矩阵热图...")
    nfull <- nrow(full_D)
    dist_mat_vis <- matrix(0L, nrow = nfull, ncol = nfull)
    
    # 计算完整设计的距离矩阵（可视化用，非并行）
    for (i in 1:(nfull-1)) {
      for (j in (i+1):nfull) {
        d <- kendall_tau_fast(full_D[i, ], full_D[j, ])
        dist_mat_vis[i, j] <- d; dist_mat_vis[j, i] <- d
      }
    }
    
    # 保存热图
    png(
      filename = file.path(OUTPUT_DIRS$heatmaps, sprintf("dist_heatmap_m%d_n%d.png", m, n)),
      width = 800, height = 800, res = 100
    )
    pheatmap(
      dist_mat_vis,
      cluster_rows = FALSE, cluster_cols = FALSE,
      main = sprintf("Kendall 距离矩阵（m=%d, n=%d）", m, n),
      xlab = "设计行索引", ylab = "设计行索引",
      color = viridis::viridis(50)  # 美观的颜色映射
    )
    dev.off()
    
    # ---------------------- 4. 输出实验总结 ----------------------
    cat(sprintf("\n✅ 实验完成（m=%d）：\n", m))
    cat(sprintf("  - 最小距离 minK = %d\n", met$minK))
    cat(sprintf("  - 平均距离 ave_K_all = %.3f\n", met$ave_K_all))
    cat(sprintf("  - 平方平均距离 ave_sqK_all = %.3f\n", met$ave_sqK_all))
   # if (!is.na(minK_lower)) {
   #   cat(sprintf("  - 理论下界：minK ≥ %.2f, ave_sqK ≥ %.4f\n", minK_lower, ave_sqK_lower))
   # }
  }
  
  # 汇总结果并返回
  res_df <- bind_rows(results_list)
  message(sprintf("\n🎉 所有实验完成！结果汇总：%s", OUTPUT_DIRS$results_csv))
  return(res_df)
}

message("✅ 实验运行器加载完成")