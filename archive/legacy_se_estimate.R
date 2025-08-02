library(survival)
library(dplyr)

# 함수: SE 비교를 위한 분석 함수
analyze_se <- function(data, scenario_name, cut_points) {
  cat("\n", rep("=", 50), "\n")
  cat("시나리오:", scenario_name, "\n")
  cat(rep("=", 50), "\n")
  
  # 원본 분석
  fit1 <- coxph(Surv(time, status) ~ age, data = data)
  
  # 구간분할
  d2 <- survSplit(Surv(time, status) ~ ., data = data, cut = cut_points, episode = "interval")
  fit2_no <- coxph(Surv(tstart, time, status) ~ age, data = d2)
  fit2_cluster <- coxph(Surv(tstart, time, status) ~ age + cluster(id), data = d2)
  
  # 클러스터 정보
  cluster_info <- d2 %>% 
    group_by(id) %>% 
    summarise(n_intervals = n(), .groups = "drop")
  
  # 결과 출력
  cat("원본 데이터 (n =", nrow(data), ") SE:", round(summary(fit1)$coefficients[1, "se(coef)"], 6), "\n")
  cat("구간분할 (n =", nrow(d2), ") 일반 SE:", round(summary(fit2_no)$coefficients[1, "se(coef)"], 6), "\n")
  cat("구간분할 (n =", nrow(d2), ") Robust SE:", round(summary(fit2_cluster)$coefficients[1, "robust se"], 6), "\n")
  cat("평균 구간수:", round(mean(cluster_info$n_intervals), 2), 
      "| 최대 구간수:", max(cluster_info$n_intervals),
      "| 구간 1개:", sum(cluster_info$n_intervals == 1), "명\n")
  
  # SE 비율 계산
  robust_ratio <- summary(fit2_cluster)$coefficients[1, "robust se"] / summary(fit1)$coefficients[1, "se(coef)"]
  cat("Robust SE / 원본 SE 비율:", round(robust_ratio, 3), 
      if(robust_ratio > 1) "(더 큼)" else "(더 작음)", "\n")
  
  return(list(
    original_se = summary(fit1)$coefficients[1, "se(coef)"],
    robust_se = summary(fit2_cluster)$coefficients[1, "robust se"],
    ratio = robust_ratio,
    avg_intervals = mean(cluster_info$n_intervals)
  ))
}

# 시나리오 1: 원래 케이스 (구간이 적음)
set.seed(123)
n <- 100
age1 <- rnorm(n, 50, 10)
ctime1 <- rexp(n, 0.1)
etime1 <- rexp(n, 0.12 * exp(-0.02 * age1))
otime1 <- pmin(etime1, ctime1)
status1 <- as.numeric(etime1 <= ctime1)
d1 <- data.frame(id = 1:n, age = age1, time = otime1, status = status1)

result1 <- analyze_se(d1, "구간 적음 (5단위 분할)", seq(5, max(otime1), by = 5))

# 시나리오 2: 더 세밀한 구간분할 (구간이 많음)
set.seed(123)
result2 <- analyze_se(d1, "구간 많음 (1단위 분할)", seq(1, max(otime1), by = 1))

# 시나리오 3: 생존시간이 긴 경우
set.seed(123)
age3 <- rnorm(n, 50, 10)
ctime3 <- rexp(n, 0.02)  # 검열시간을 늘림
etime3 <- rexp(n, 0.03 * exp(-0.02 * age3))  # 사건시간도 늘림
otime3 <- pmin(etime3, ctime3)
status3 <- as.numeric(etime3 <= ctime3)
d3 <- data.frame(id = 1:n, age = age3, time = otime3, status = status3)

result3 <- analyze_se(d3, "생존시간 긴 경우 (2단위 분할)", seq(2, max(otime3), by = 2))

# 시나리오 4: 강한 효과가 있는 경우
set.seed(123)
age4 <- rnorm(n, 50, 10)
ctime4 <- rexp(n, 0.1)
etime4 <- rexp(n, 0.12 * exp(-0.1 * age4))  # 연령 효과를 5배 크게
otime4 <- pmin(etime4, ctime4)
status4 <- as.numeric(etime4 <= ctime4)
d4 <- data.frame(id = 1:n, age = age4, time = otime4, status = status4)

result4 <- analyze_se(d4, "강한 연령효과 (2단위 분할)", seq(2, max(otime4), by = 2))

# 전체 비교표
cat("\n", rep("=", 60), "\n")
cat("전체 시나리오 비교\n")
cat(rep("=", 60), "\n")

comparison <- data.frame(
  시나리오 = c("구간 적음", "구간 많음", "생존시간 긴 경우", "강한 연령효과"),
  원본_SE = c(result1$original_se, result2$original_se, result3$original_se, result4$original_se),
  Robust_SE = c(result1$robust_se, result2$robust_se, result3$robust_se, result4$robust_se),
  비율 = c(result1$ratio, result2$ratio, result3$ratio, result4$ratio),
  평균구간수 = c(result1$avg_intervals, result2$avg_intervals, result3$avg_intervals, result4$avg_intervals)
)

comparison$원본_SE <- round(comparison$원본_SE, 4)
comparison$Robust_SE <- round(comparison$Robust_SE, 4)
comparison$비율 <- round(comparison$비율, 4)
comparison$평균구간수 <- round(comparison$평균구간수, 2)

print(comparison)
