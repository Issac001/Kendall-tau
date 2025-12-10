# 01_config.R: 全局参数配置（可按需修改）

# ---------------------- 1. 实验范围配置 ----------------------
m_range <- 26:32  # 因子数量范围（仅处理偶数 m，因 n 需为偶数）
# 总设计行数 n 与 m 的关系：n = m（代表集大小 K = n/2 = m/2）
total_n_from_m <- function(m) {
  if (m %% 2 == 0) return(m) else return(NA_integer_)
}

# ---------------------- 2. 模拟退火（SA）参数 ----------------------
SA_PARAMS <- list(
  T0 = 1.0,                # 初始温度
  alpha = 0.997,           # 降温系数（每轮迭代 × alpha）
  max_iter = 6000,         # 每轮重启的最大迭代次数
  restarts = 8,            # 重启次数（避免局部最优）
  large_move_prob = 0.08,  # 大跳跃（替换行）概率
  eps_secondary = 1e-4,    # 次目标（ave_sqK）权重（破 ties）
  no_improve_stop = 800    # 早停阈值（无改进迭代次数）
)

# ---------------------- 3. 输出路径配置 ----------------------
# OUTPUT_DIRS <- list(
#   root = "outputs",                # 结果根目录
#   best_designs = "outputs/best_designs",  # 最佳设计 RData
#   heatmaps = "outputs/heatmaps",   # 热图 PNG
#   results_csv = "outputs/results.csv"     # 结果汇总 CSV
# )
PROJECT_ROOT <- "~/kendall_tau_design"
OUTPUT_DIRS <- list(
  base = file.path(PROJECT_ROOT, "outputs"),
  best_designs = file.path(PROJECT_ROOT, "outputs/best_designs"),
  heatmaps = file.path(PROJECT_ROOT, "outputs/heatmaps"),
  results_csv = file.path(PROJECT_ROOT, "outputs/results.csv")
)

# 创建输出目录（若不存在）
walk(OUTPUT_DIRS, ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))

message("✅ 全局参数配置完成，输出目录已创建")