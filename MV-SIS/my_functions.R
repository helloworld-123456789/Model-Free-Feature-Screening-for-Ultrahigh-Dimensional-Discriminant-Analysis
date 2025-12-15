# MV-SIS
#——————————————————————————————————————————————————————————————
# 计算 MV(X_k | Y) 的样本估计值
compute_mv_index <- function(X, Y) {
  # X: n x p 矩阵
  # Y: 长度为 n 的因子或整数向量（类别标签，从1开始）
  n <- length(Y)
  p <- ncol(X)
  R <- length(unique(Y))
  
  # 类别概率
  pr <- tabulate(Y) / n  # 假设 Y 是 1,2,...,R
  
  mv_scores <- numeric(p)
  
  for (k in 1:p) {
    x <- X[, k]
    Fx <- ecdf(x)(x)  # F(x_i) = P(X <= x_i)，经验分布
    
    # 对每个类别 r，计算条件分布 Fr(x_i)
    weighted_var <- 0
    for (r in 1:R) {
      idx_r <- which(Y == r)
      if (length(idx_r) == 0) next
      Fr_x <- ecdf(x[idx_r])(x)  # Fr(x_i) = P(X <= x_i | Y=r)
      diff_sq <- (Fr_x - Fx)^2
      weighted_var <- weighted_var + pr[r] * mean(diff_sq)
    }
    mv_scores[k] <- weighted_var
  }
  return(mv_scores)
}
#——————————————————————————————————————————————————————————————



# PSIS
#——————————————————————————————————————————————————————————————
compute_psis_index <- function(X, Y, use_median = FALSE) {
  # 对二分类，PSIS 就是 |mean(X[Y=1]) - mean(X[Y=2])|
  # use_median: 使用中位数代替均值（对重尾分布更鲁棒）
  n <- length(Y)
  p <- ncol(X)
  
  y1_idx <- which(Y == 1)
  y2_idx <- which(Y == 2)
  
  if (use_median) {
    mean1 <- colMedians(X[y1_idx, , drop = FALSE])
    mean2 <- colMedians(X[y2_idx, , drop = FALSE])
    # 中位数差 (鲁棒)
    psis_scores <- abs(mean1 - mean2)
  } else {
    # 使用 t-统计量 (标准 SIS 做二分类筛选的做法)
    # t = |mean1 - mean2| / sqrt(s1^2/n1 + s2^2/n2)
    n1 <- length(y1_idx)
    n2 <- length(y2_idx)
    
    mean1 <- colMeans(X[y1_idx, , drop = FALSE])
    mean2 <- colMeans(X[y2_idx, , drop = FALSE])
    
    # 需要 matrixStats::colVars
    var1 <- colVars(X[y1_idx, , drop = FALSE])
    var2 <- colVars(X[y2_idx, , drop = FALSE])
    
    # 计算标准误，加一个小常数防止除零
    se <- sqrt(var1/n1 + var2/n2 + 1e-10)
    
    psis_scores <- abs(mean1 - mean2) / se
  }
  return(psis_scores)
}
#——————————————————————————————————————————————————————————————


# 样本生成
#——————————————————————————————————————————————————————————————
generate_data <- function(n, p, R = 2, balanced = TRUE, heavy_tail = FALSE) {
  # 生成 Y
  if (balanced) {
    Y <- sample(1:R, size = n, replace = TRUE)
  } else {
    # 不平衡：等差数列概率，R=2 时为 1/3, 2/3
    if (R == 2) {
      probs <- c(1/3, 2/3)
    } else {
      # 一般化：p_r = 2[1 + (r-1)/(R-1)] / (3R)
      probs <- 2 * (1 + (0:(R-1)) / (R - 1)) / (3 * R)
    }
    Y <- sample(1:R, size = n, replace = TRUE, prob = probs)
  }
  
  # 生成 X
  X <- matrix(0, nrow = n, ncol = p)
  
  # 活跃变量：前 R 个（这里 R=2）
  active_idx <- 1:R
  
  # 均值向量：mu_r 只在第 r 位为 3
  mu_matrix <- matrix(0, nrow = R, ncol = p)
  for (r in 1:R) {
    mu_matrix[r, r] <- 3
  }
  
  # 误差项
  if (heavy_tail) {
    eps <- rt(n * p, df = 2)
  } else {
    eps <- rnorm(n * p)
  }
  eps <- matrix(eps, nrow = n, ncol = p)
  
  # 逐样本生成 X_i
  for (i in 1:n) {
    r <- Y[i]
    X[i, ] <- mu_matrix[r, ] + eps[i, ]
  }
  
  return(list(X = X, Y = Y, active = active_idx))
}
#——————————————————————————————————————————————————————————————

# 改进的PSIS (多分类)
#——————————————————————————————————————————————————————————————
compute_psis_star <- function(X, Y, use_median = FALSE) {
  # X: n x p matrix
  # Y: vector of class labels (e.g., 1, 2, ..., R)
  # use_median: 使用中位数代替均值（对重尾分布更鲁棒）
  
  n <- nrow(X)
  p <- ncol(X)
  classes <- sort(unique(Y))
  R <- length(classes)
  
  if (use_median) {
    # 中位数情况：使用最大成对中位数差
    class_medians <- sapply(classes, function(r) {
      colMedians(X[Y == r, , drop = FALSE])
    })  # p x R matrix
    
    psis_star_scores <- numeric(p)
    for (k in 1:p) {
      diffs <- abs(outer(class_medians[k, ], class_medians[k, ], "-"))
      psis_star_scores[k] <- max(diffs[upper.tri(diffs)])
    }
    
  } else {
    # 均值情况：使用最大成对 t-统计量 (Max Pairwise t-statistic)
    # 预计算每个类的统计量
    class_stats <- lapply(classes, function(r) {
      idx <- which(Y == r)
      list(
        n = length(idx),
        mean = colMeans(X[idx, , drop = FALSE]),
        var = colVars(X[idx, , drop = FALSE])
      )
    })
    
    # 初始化分数为0
    psis_star_scores <- numeric(p)
    
    # 遍历所有类别对 (i, j)
    for (i in 1:(R-1)) {
      for (j in (i+1):R) {
        stats_i <- class_stats[[i]]
        stats_j <- class_stats[[j]]
        
        # 计算 t 统计量向量 (长度为 p)
        se <- sqrt(stats_i$var / stats_i$n + stats_j$var / stats_j$n + 1e-10)
        diff <- abs(stats_i$mean - stats_j$mean)
        t_stats <- diff / se
        
        # 更新每个变量的最大 t 值
        psis_star_scores <- pmax(psis_star_scores, t_stats)
      }
    }
  }
  
  return(psis_star_scores)
}
#——————————————————————————————————————————————————————————————


