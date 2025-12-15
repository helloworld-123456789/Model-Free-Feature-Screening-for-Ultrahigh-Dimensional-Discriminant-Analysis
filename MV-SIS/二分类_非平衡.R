library(matrixStats)  # 用于 rowMads 等

# 导入函数文件
source("C:\\Users\\17366\\Desktop\\MV-SIS\\my_functions.R")

set.seed(123)

n_sim <- 100
n <- 40
p <- 2000
R <- 2
d <- floor(n / log(n))

res_list <- list()

for (err_type in c("normal", "heavy")) {
  heavy_tail <- (err_type == "heavy")
  
  mms_mv <- mms_psis <- numeric(n_sim)
  p1_mv <- p2_mv <- p1_psis <- p2_psis <- numeric(n_sim)
  pall_mv <- pall_psis <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    dat <- generate_data(n = n, p = p, R = R, balanced = FALSE, heavy_tail = heavy_tail)
    X <- dat$X
    Y <- dat$Y
    active <- dat$active  # c(1,2)
    
    # MV-SIS
    mv_scores <- compute_mv_index(X, Y)
    rank_mv <- order(mv_scores, decreasing = TRUE)
    top_d_mv <- rank_mv[1:d]
    mms_mv[sim] <- max(sapply(active, function(a) which(rank_mv == a)))  # 包含所有活跃变量所需的最小模型大小
    p1_mv[sim] <- as.numeric(1 %in% top_d_mv)
    p2_mv[sim] <- as.numeric(2 %in% top_d_mv)
    pall_mv[sim] <- as.numeric(all(active %in% top_d_mv))
    
    # PSIS
    psis_scores <- compute_psis_index(X, Y, use_median = FALSE)
    rank_psis <- order(psis_scores, decreasing = TRUE)
    top_d_psis <- rank_psis[1:d]
    mms_psis[sim] <- max(sapply(active, function(a) which(rank_psis == a)))
    p1_psis[sim] <- as.numeric(1 %in% top_d_psis)
    p2_psis[sim] <- as.numeric(2 %in% top_d_psis)
    pall_psis[sim] <- as.numeric(all(active %in% top_d_psis))
  }
  
  res_list[[paste0("MV-SIS_", err_type)]] <- data.frame(
    method = "MV-SIS", error_type = err_type,
    mms = median(mms_mv), 
    p1 = mean(p1_mv), p2 = mean(p2_mv), pall = mean(pall_mv)
  )
  
  res_list[[paste0("PSIS_", err_type)]] <- data.frame(
    method = "PSIS", error_type = err_type,
    mms = median(mms_psis), 
    p1 = mean(p1_psis), p2 = mean(p2_psis), pall = mean(pall_psis)
  )
}

final_results <- do.call(rbind, res_list)
rownames(final_results) <- NULL

cat("\n========================================\n")
cat("实验一：二分类非平衡场景 (n=40, p=2000, R=2)\n")
cat("========================================\n")
cat("模拟次数:", n_sim, "\n")
cat("筛选变量数 d =", d, "\n\n")

print(final_results, row.names = FALSE)

cat("\n指标说明:\n")
cat("- mms: 包含所有活跃变量的最小模型大小(中位数)\n")
cat("- p1:  变量1被选中的概率\n")
cat("- p2:  变量2被选中的概率\n")
cat("- pall: 所有活跃变量都被选中的概率\n")

