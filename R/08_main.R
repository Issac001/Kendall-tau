# 08_main.R: 主脚本（一键启动所有实验）
cat("=======================================\n")
cat("          SA Foldover Design          \n")
cat("=======================================\n")

options(error = recover)
# 按依赖顺序加载所有模块
source("R/00_dependencies.R")
source("R/01_config.R")
source("R/02_utils.R")
source("R/03_distance_matrix.R")
source("R/04_metrics.R")
source("R/05_neighborhood.R")
source("R/06_sa_engine.R")
source("R/07_experiment_runner.R")

# 运行实验
res_df <- run_experiments()

# 打印最终结果汇总
cat("\n=======================================\n")
cat("            实验结果汇总               \n")
cat("=======================================\n")
print(res_df %>% select(-timestamp))  # 隐藏时间戳，简化输出

cat("\n✅ 所有结果已保存至 outputs/ 目录！\n")