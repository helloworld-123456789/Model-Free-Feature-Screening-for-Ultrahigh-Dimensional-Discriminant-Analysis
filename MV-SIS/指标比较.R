library(MASS)
library(stats)
library(pracma)
library(matrixStats)

# 导入函数文件
source("C:\\Users\\17366\\Desktop\\MV-SIS\\my_functions.R")


# 设置随机种子
set.seed(123)



# ==============================================================================
# 主程序
# ==============================================================================

# 参数设置
n_sim <- 100       # 模拟次数
p <- 2000          # 特征维数
R <- 2             # 类别数
n_values <- c(40, 80) # 样本量情况

# 真实参数 (对应 my_functions.R 中的生成逻辑)
# mu1: 第1位为3; mu2: 第2位为3
mu1_true <- numeric(p); mu1_true[1] <- 3
mu2_true <- numeric(p); mu2_true[2] <- 3
gamma0 <- mu1_true - mu2_true
active_true <- c(1, 2) # 真实变量索引
s_n <- length(active_true)

res_list <- list()

# 循环：样本量 -> 误差类型
for (n in n_values) {
  d_screen <- floor(n / log(n)) # 筛选阶段保留的变量数
  
  for (err_type in c("normal", "heavy")) {
    heavy_tail <- (err_type == "heavy")
    
    # 存储单次模拟的指标
    # 维度: n_sim x 8 (MS, CZ, IZ, CP, RSSE, CA, CA0, RCA)
    col_names <- c("MS", "CZ", "IZ", "CP", "RSSE", "CA", "CA0", "RCA")
    metrics_mv <- matrix(NA, n_sim, 8)
    metrics_psis <- matrix(NA, n_sim, 8)
    colnames(metrics_mv) <- col_names
    colnames(metrics_psis) <- col_names
    
    cat(sprintf("正在模拟: n=%d, Error=%s ...\n", n, err_type))
    
    for (sim in 1:n_sim) {
      # 1. 生成数据 (训练集 & 测试集)
      # Paper: "an independent testing dataset is generated with the same sample size"
      dat_train <- generate_data(n = n, p = p, R = R, balanced = TRUE, heavy_tail = heavy_tail)
      dat_test <- generate_data(n = n, p = p, R = R, balanced = TRUE, heavy_tail = heavy_tail)
      
      X_tr <- dat_train$X; Y_tr <- dat_train$Y
      X_te <- dat_test$X;  Y_te <- dat_test$Y
      
      # 计算基准准确率 CA0 (Oracle)
      ca0 <- compute_ca0(X_te, Y_te, mu1_true, mu2_true)
      
      # ==========================
      # 方法 1: MV-SIS
      # ==========================
      # (1) 筛选
      mv_scores <- compute_mv_index(X_tr, Y_tr)
      rank_mv <- order(mv_scores, decreasing = TRUE)
      top_screen_mv <- rank_mv[1:d_screen] # 初步筛选
      
      # (2) 模型选择 (BIC)
      sel_idx_mv <- select_model_bic(X_tr, Y_tr, top_screen_mv, max_d = d_screen)
      
      # (3) 估计与评估
      eval_mv <- evaluate_model(X_tr, Y_tr, X_te, Y_te, sel_idx_mv, p)
      
      # (4) 计算指标
      # MS: Model Size
      metrics_mv[sim, "MS"] <- length(sel_idx_mv)
      
      # CZ: Correct Zeros % = (p - |D U D*|) / (p - sn)
      # 正确排除的非活跃变量比例
      union_card <- length(union(active_true, sel_idx_mv))
      metrics_mv[sim, "CZ"] <- (p - union_card) / (p - s_n) * 100
      
      # IZ: Incorrect Zeros % = |D*c n D| / |D|
      # 错误排除的活跃变量比例 (False Negative Rate)
      missed_active <- sum(!(active_true %in% sel_idx_mv))
      metrics_mv[sim, "IZ"] <- missed_active / s_n * 100
      
      # CP: Coverage Probability (是否包含所有真实变量)
      metrics_mv[sim, "CP"] <- as.numeric(all(active_true %in% sel_idx_mv)) * 100
      
      # RSSE: ||gamma_hat - gamma0||
      metrics_mv[sim, "RSSE"] <- sqrt(sum((eval_mv$gamma_hat - gamma0)^2))
      
      # CA: Classification Accuracy
      metrics_mv[sim, "CA"] <- eval_mv$acc * 100
      
      # CA0: Oracle Accuracy
      metrics_mv[sim, "CA0"] <- ca0 * 100
      
      # RCA: Relative CA
      metrics_mv[sim, "RCA"] <- (eval_mv$acc / ca0) * 100
      
      
      # ==========================
      # 方法 2: PSIS
      # ==========================
      # (1) 筛选
      psis_scores <- compute_psis_index(X_tr, Y_tr, use_median = FALSE)
      rank_psis <- order(psis_scores, decreasing = TRUE)
      top_screen_psis <- rank_psis[1:d_screen]
      
      # (2) 模型选择 (BIC)
      sel_idx_psis <- select_model_bic(X_tr, Y_tr, top_screen_psis, max_d = d_screen)
      
      # (3) 估计与评估
      eval_psis <- evaluate_model(X_tr, Y_tr, X_te, Y_te, sel_idx_psis, p)
      
      # (4) 计算指标
      metrics_psis[sim, "MS"] <- length(sel_idx_psis)
      
      union_card_psis <- length(union(active_true, sel_idx_psis))
      metrics_psis[sim, "CZ"] <- (p - union_card_psis) / (p - s_n) * 100
      
      missed_active_psis <- sum(!(active_true %in% sel_idx_psis))
      metrics_psis[sim, "IZ"] <- missed_active_psis / s_n * 100
      
      metrics_psis[sim, "CP"] <- as.numeric(all(active_true %in% sel_idx_psis)) * 100
      
      metrics_psis[sim, "RSSE"] <- sqrt(sum((eval_psis$gamma_hat - gamma0)^2))
      
      metrics_psis[sim, "CA"] <- eval_psis$acc * 100
      metrics_psis[sim, "CA0"] <- ca0 * 100
      metrics_psis[sim, "RCA"] <- (eval_psis$acc / ca0) * 100
    }
    
    # 汇总结果
    # MS: Median (RSD)
    # Others: Mean
    
    calc_stats <- function(mat, method_name) {
      ms_med <- median(mat[, "MS"])
      ms_rsd <- IQR(mat[, "MS"]) / 1.349
      
      data.frame(
        Method = method_name,
        n = n,
        Error = err_type,
        `MS(RSD)` = sprintf("%.1f (%.1f)", ms_med, ms_rsd),
        `CZ(%)` = mean(mat[, "CZ"]),
        `IZ(%)` = mean(mat[, "IZ"]),
        `CP(%)` = mean(mat[, "CP"]),
        RSSE = mean(mat[, "RSSE"]),
        `CA(%)` = mean(mat[, "CA"]),
        `CA0(%)` = mean(mat[, "CA0"]),
        RCA = mean(mat[, "RCA"]),
        check.names = FALSE
      )
    }
    
    # 注意：Table 3 先列出 PSIS 再列出 MV-SIS
    res_list[[paste0("PSIS_", n, "_", err_type)]] <- calc_stats(metrics_psis, "PSIS")
    res_list[[paste0("MV-SIS_", n, "_", err_type)]] <- calc_stats(metrics_mv, "MV-SIS")
  }
}

# ==============================================================================
# 结果展示
# ==============================================================================
final_results <- do.call(rbind, res_list)
rownames(final_results) <- NULL

# 格式化输出
cat("\n================================================================================\n")
cat("Table 3 仿真结果复现 (n=40/80, p=2000, R=2, 500 reps)\n")
cat("================================================================================\n")
print(final_results, digits = 3)

cat("\n指标说明:\n")
cat("MS(RSD): 模型大小中位数 (鲁棒标准差)\n")
cat("CZ(%)  : 正确排除非重要变量比例 (Mean)\n")
cat("IZ(%)  : 错误排除重要变量比例 (Mean)\n")
cat("CP(%)  : 覆盖概率 (Mean)\n")
cat("RSSE   : 估计误差 (Mean)\n")
cat("CA(%)  : 分类准确率 (Mean)\n")
cat("CA0(%) : Oracle准确率 (Mean)\n")
cat("RCA    : 相对准确率 (Mean)\n")
