source("C:\\Users\\17366\\Desktop\\MV-SIS\\my_functions.R")

set.seed(123)

n_sim <- 100
n <- 200
p <- 2000
R <- 10
d <- floor(n / log(n))

res_list <- list()

for (err_type in c("normal", "heavy")) {
  heavy_tail <- (err_type == "heavy")
  
  mms_mv <- mms_psis <- numeric(n_sim)
  p_mv <- matrix(0, nrow = n_sim, ncol = R)
  p_psis <- matrix(0, nrow = n_sim, ncol = R)
  pall_mv <- pall_psis <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    dat <- generate_data(n = n, p = p, R = R, balanced = FALSE, heavy_tail = heavy_tail)
    X <- dat$X
    Y <- dat$Y
    active <- dat$active  # c(1,2,...,10)
    
    # MV-SIS
    mv_scores <- compute_mv_index(X, Y)
    rank_mv <- order(mv_scores, decreasing = TRUE)
    top_d_mv <- rank_mv[1:d]
    mms_mv[sim] <- max(sapply(active, function(a) which(rank_mv == a)))  # 包含所有活跃变量所需的最小模型大小
    
    for (r in 1:R) {
      p_mv[sim, r] <- as.numeric(r %in% top_d_mv)
    }
    pall_mv[sim] <- as.numeric(all(active %in% top_d_mv))
    
    # PSIS
    psis_scores <- compute_psis_star(X, Y, use_median = FALSE)
    rank_psis <- order(psis_scores, decreasing = TRUE)
    top_d_psis <- rank_psis[1:d]
    mms_psis[sim] <- max(sapply(active, function(a) which(rank_psis == a)))
    
    for (r in 1:R) {
      p_psis[sim, r] <- as.numeric(r %in% top_d_psis)
    }
    pall_psis[sim] <- as.numeric(all(active %in% top_d_psis))
  }
  
  # 汇总 MV-SIS 结果
  p_means_mv <- colMeans(p_mv)
  names(p_means_mv) <- paste0("p", 1:R)
  
  df_mv <- data.frame(
    method = "MV-SIS", error_type = err_type,
    mms = median(mms_mv)
  )
  df_mv <- cbind(df_mv, t(p_means_mv))
  df_mv$pall <- mean(pall_mv)
  res_list[[paste0("MV-SIS_", err_type)]] <- df_mv
  
  # 汇总 PSIS 结果
  p_means_psis <- colMeans(p_psis)
  names(p_means_psis) <- paste0("p", 1:R)
  
  df_psis <- data.frame(
    method = "PSIS", error_type = err_type,
    mms = median(mms_psis)
  )
  df_psis <- cbind(df_psis, t(p_means_psis))
  df_psis$pall <- mean(pall_psis)
  res_list[[paste0("PSIS_", err_type)]] <- df_psis
}

final_results <- do.call(rbind, res_list)
rownames(final_results) <- NULL

cat("\n========================================\n")
cat("实验四：十分类非平衡场景 (n=200, p=2000, R=10)\n")
cat("========================================\n")
cat("模拟次数:", n_sim, "\n")
cat("筛选变量数 d =", d, "\n\n")

print(final_results, row.names = FALSE)

cat("\n指标说明:\n")
cat("- mms: 包含所有活跃变量的最小模型大小(中位数)\n")
cat("- pall: 所有活跃变量都被选中的概率\n")
cat("- p1-p10: 各个活跃变量被选中的概率\n")

