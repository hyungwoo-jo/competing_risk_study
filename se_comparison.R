# SE 비교 및 시각화
library(survival)
library(dplyr)

set.seed(123)
n <- 100
age <- rnorm(n, 50, 10)
ctime <- rexp(n, 0.1)
etime <- rexp(n, 0.12 * exp(-0.02 * age))
otime <- pmin(etime, ctime)
status <- as.numeric(etime <= ctime)

# 원본 데이터
d1 <- data.frame(id = 1:n, age = age, time = otime, status = status)

# 구간분할 데이터
d2 <- survSplit(Surv(time, status) ~ ., data = d1, cut = seq(5, max(otime), by = 5), episode = "interval")

# 모델 비교
fit1 <- coxph(Surv(time, status) ~ age, data = d1)
fit2_no <- coxph(Surv(tstart, time, status) ~ age, data = d2)
fit2_cluster <- coxph(Surv(tstart, time, status) ~ age + cluster(id), data = d2)

# SE 비교표
se_comparison <- data.frame(
  Model = c("원본 (n=100)", "구간분할 (n=193)", "구간분할+클러스터 (n=193)"),
  SE = c(
    summary(fit1)$coefficients[1, "se(coef)"],
    summary(fit2_no)$coefficients[1, "se(coef)"],
    summary(fit2_cluster)$coefficients[1, "robust se"]
  ),
  Type = c("일반", "일반", "Robust")
)

print("=== 표준오차 비교 ===")
print(se_comparison)

# 클러스터 크기 확인
cluster_info <- d2 %>% 
  group_by(id) %>% 
  summarise(n_intervals = n(), .groups = "drop")

cat("\n=== 클러스터 정보 ===")
cat("\n평균 구간 수:", mean(cluster_info$n_intervals))
cat("\n최대 구간 수:", max(cluster_info$n_intervals))
cat("\n구간이 1개인 사람:", sum(cluster_info$n_intervals == 1), "명")