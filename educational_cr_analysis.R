# =============================================================================
# Educational Competing Risks Analysis: SE Differences & Robust SE Phenomena
# =============================================================================
# 이 코드는 Perplexity 대화에서 논의된 핵심 개념들을 실제로 체험할 수 있게 설계됨
# 주요 학습 목표:
# 1. Fine-Gray와 CRR의 차이점 이해
# 2. Robust SE가 더 작아지는 현상 관찰 
# 3. Tie 처리 방식의 영향 확인
# 4. Standard vs Robust SE의 실제적 차이 비교
# 5. True value 기반 교육적 시뮬레이션

library(survival)
library(cmprsk)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(knitr)
library(parallel)

# 시스템 코어 수 확인 및 설정
total_cores <- detectCores()
n_cores <- max(1, total_cores - 2)  # 전체 코어 - 2
cat(sprintf("🖥️  시스템 정보: 총 %d코어, 사용할 코어: %d개\n\n", total_cores, n_cores))

set.seed(2024)

# =============================================================================
# 1. 데이터 생성 함수 (교육적 목적으로 단순화)
# =============================================================================

generate_educational_data <- function(n = 300, 
                                    beta_age = 0.5, 
                                    beta_sex = -0.3, 
                                    beta_score = 0.2,
                                    tie_probability = 0.3,  # tie 발생 확률
                                    censoring_rate = 0.15) {
  # 공변량 생성 (표준화)
  age <- scale(rnorm(n, 65, 12))[,1]        # 나이 (표준화)
  sex <- rbinom(n, 1, 0.5)                  # 성별 (0=여성, 1=남성)
  score <- scale(rnorm(n, 75, 20))[,1]      # 건강점수 (표준화)
  
  # Primary event (사망) hazard
  lp_primary <- beta_age * age + beta_sex * sex + beta_score * score
  time_primary <- rexp(n, exp(lp_primary))
  
  # Competing event (전원) hazard  
  lp_competing <- 0.3 * age - 0.15 * sex + 0.1 * score
  time_competing <- rexp(n, exp(lp_competing))
  
  # 검열 시간
  cens_time <- rexp(n, censoring_rate)
  
  # 관찰 시간과 사건 상태 결정
  obs_time <- pmin(time_primary, time_competing, cens_time)
  status <- ifelse(time_primary <= pmin(time_competing, cens_time), 1,
                   ifelse(time_competing <= pmin(time_primary, cens_time), 2, 0))
  
  # Tie 인위적 생성 (교육 목적)
  tie_indices <- sample(1:n, round(n * tie_probability))
  if(length(tie_indices) > 1) {
    # 일부 관측치들을 같은 시간으로 설정
    tie_time <- median(obs_time[tie_indices])
    obs_time[tie_indices] <- tie_time
  }
  
  data.frame(
    id = 1:n,
    time = obs_time,
    status = status,
    age = age,
    sex = sex,
    score = score
  )
}

# =============================================================================
# True Value 계산 함수 (대용량 데이터셋 기반)
# =============================================================================

get_true_values <- function(beta_age = 0.5, 
                           beta_sex = -0.3, 
                           beta_score = 0.2,
                           large_n = 50000,
                           tie_method = "efron") {
  
  cat("🎯 True Values 계산 중...\n")
  cat(sprintf("   - 대용량 데이터셋 크기: %d\n", large_n))
  cat("   - Tie 없는 연속시간 사용\n")
  cat("   - CRR 방법으로 true value 추정\n\n")
  
  # 대용량 데이터 생성 (tie 없음)
  large_data <- generate_educational_data(
    n = large_n,
    beta_age = beta_age,
    beta_sex = beta_sex, 
    beta_score = beta_score,
    tie_probability = 0,  # tie 없음
    censoring_rate = 0.1  # 낮은 검열률
  )
  
  cat(sprintf("대용량 데이터 생성 완료:\n"))
  cat(sprintf("  - Primary events: %d (%.1f%%)\n", 
              sum(large_data$status == 1), mean(large_data$status == 1) * 100))
  cat(sprintf("  - Competing events: %d (%.1f%%)\n", 
              sum(large_data$status == 2), mean(large_data$status == 2) * 100))
  cat(sprintf("  - Censored: %d (%.1f%%)\n\n", 
              sum(large_data$status == 0), mean(large_data$status == 0) * 100))
  
  # CRR로 true value 추정
  cov_matrix <- as.matrix(large_data[, c("age", "sex", "score")])
  
  crr_fit <- crr(ftime = large_data$time, 
                 fstatus = large_data$status, 
                 cov1 = cov_matrix,
                 failcode = 1, 
                 cencode = 0)
  
  true_values <- crr_fit$coef
  names(true_values) <- c("age", "sex", "score")
  
  cat("📊 추정된 True Values:\n")
  for(i in 1:length(true_values)) {
    cat(sprintf("   %s: %.6f\n", names(true_values)[i], true_values[i]))
  }
  cat("\n")
  
  return(list(
    values = true_values,
    large_data = large_data,
    crr_fit = crr_fit
  ))
}

# =============================================================================
# 2. 모델 적합 함수 - 4가지 접근법 비교
# =============================================================================

fit_educational_models <- function(data, tie_method = "efron") {
  
  # 결과 저장용 리스트
  results <- list()
  
  cat("=== 4가지 접근법으로 Competing Risks 분석 ===\n\n")
  
  # -------------------------------
  # 1. cmprsk::crr() - 직접 구현
  # -------------------------------
  cat("1. cmprsk::crr() 접근법\n")
  cat("   - 직접적인 Fine-Gray 모델 구현\n")
  cat("   - Breslow tie 처리 방식 사용\n")
  
  cov_matrix <- as.matrix(data[, c("age", "sex", "score")])
  
  crr_fit <- crr(ftime = data$time, 
                 fstatus = data$status, 
                 cov1 = cov_matrix,
                 failcode = 1, 
                 cencode = 0)
  
  results$crr <- list(
    method = "cmprsk::crr",
    coef = crr_fit$coef,
    se = sqrt(diag(crr_fit$var)),
    tie_method = "breslow"
  )
  
  # -------------------------------
  # 2. finegray + coxph (standard SE)
  # -------------------------------
  cat("2. finegray + coxph (Standard SE)\n")
  cat("   - 데이터를 counting process 형태로 변환\n")
  cat("   - Model-based 표준오차 사용\n")
  
  status_factor <- factor(data$status, levels = c(0,1,2), 
                         labels = c("cens", "primary", "competing"))
  
  fg_data <- finegray(Surv(time, status_factor) ~ ., 
                      data = data, 
                      etype = "primary")
  
  # ID 추가 (clustering용)
  fg_data$subject_id <- fg_data$id
  
  fg_fit_standard <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + score,
                          data = fg_data,
                          weights = fgwt,
                          ties = tie_method,
                          id = subject_id)
  
  results$fg_standard <- list(
    method = "finegray + coxph (Standard SE)",
    coef = coef(fg_fit_standard),
    se = summary(fg_fit_standard)$coefficients[, "se(coef)"],
    tie_method = tie_method
  )
  
  # -------------------------------
  # 3. finegray + coxph (robust SE, no cluster)
  # -------------------------------
  cat("3. finegray + coxph (Robust SE, no cluster)\n")
  cat("   - Sandwich 추정량 사용\n")
  cat("   - Cluster 보정 없음\n")
  
  fg_fit_robust <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + score,
                        data = fg_data,
                        weights = fgwt,
                        ties = tie_method,
                        robust = TRUE,
                        id = subject_id)
  
  results$fg_robust <- list(
    method = "finegray + coxph (Robust SE, no cluster)",
    coef = coef(fg_fit_robust),
    se = summary(fg_fit_robust)$coefficients[, "robust se"],
    tie_method = tie_method
  )
  
  # -------------------------------
  # 4. finegray + coxph (robust SE + cluster)
  # -------------------------------
  cat("4. finegray + coxph (Robust SE + cluster)\n")
  cat("   - Sandwich 추정량 + Cluster 보정\n")
  cat("   - Subject-level 독립성 반영\n\n")
  
  fg_fit_cluster <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + score + cluster(subject_id),
                         data = fg_data,
                         weights = fgwt,
                         ties = tie_method)
  
  results$fg_cluster <- list(
    method = "finegray + coxph (Robust SE + cluster)",
    coef = coef(fg_fit_cluster),
    se = summary(fg_fit_cluster)$coefficients[, "robust se"],
    tie_method = tie_method
  )
  
  return(results)
}

# =============================================================================
# 3. 결과 비교 및 시각화 함수
# =============================================================================

compare_results <- function(results) {
  
  # 결과를 data.frame으로 정리
  comparison_df <- data.frame()
  
  for(method_name in names(results)) {
    method_result <- results[[method_name]]
    
    for(i in 1:length(method_result$coef)) {
      param_name <- names(method_result$coef)[i]
      
      comparison_df <- rbind(comparison_df, data.frame(
        Method = method_result$method,
        Parameter = param_name,
        Coefficient = method_result$coef[i],
        SE = method_result$se[i],
        Tie_Method = method_result$tie_method,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 결과 출력
  cat("=== 계수 추정치 비교 ===\n")
  print(comparison_df %>%
        select(Method, Parameter, Coefficient) %>%
        pivot_wider(names_from = Method, values_from = Coefficient) %>%
        kable(digits = 4))
  
  cat("\n=== 표준오차 비교 ===\n")
  se_comparison <- comparison_df %>%
    select(Method, Parameter, SE) %>%
    pivot_wider(names_from = Method, values_from = SE)
  
  print(kable(se_comparison, digits = 4))
  
  # SE 비율 계산 (CRR 대비)
  cat("\n=== 표준오차 비율 (CRR 대비) ===\n")
  se_wide <- se_comparison
  crr_se <- se_wide$`cmprsk::crr`
  
  se_ratios <- se_wide
  for(col in 2:ncol(se_wide)) {
    if(names(se_wide)[col] != "cmprsk::crr") {
      se_ratios[[col]] <- se_wide[[col]] / crr_se
    }
  }
  
  print(kable(se_ratios, digits = 3))
  
  # 시각화
  p1 <- ggplot(comparison_df, aes(x = Parameter, y = Coefficient, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    theme_minimal() +
    labs(title = "계수 추정치 비교", 
         subtitle = "4가지 접근법별 회귀계수") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(comparison_df, aes(x = Parameter, y = SE, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    theme_minimal() +
    labs(title = "표준오차 비교", 
         subtitle = "Robust SE가 더 작아지는 현상 관찰") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(
    data = comparison_df,
    plots = list(coef = p1, se = p2)
  ))
}

# =============================================================================
# 4. 교육적 해석 함수
# =============================================================================

educational_interpretation <- function(comparison_data) {
  
  cat("\n" %+% paste(rep("=", 70), collapse="") %+% "\n")
  cat("📚 교육적 해석 및 핵심 포인트\n")
  cat(paste(rep("=", 70), collapse="") %+% "\n\n")
  
  # SE 패턴 분석
  se_data <- comparison_data %>%
    select(Method, Parameter, SE) %>%
    pivot_wider(names_from = Method, values_from = SE)
  
  cat("🔍 관찰해야 할 핵심 현상들:\n\n")
  
  cat("1. 📉 ROBUST SE가 더 작아지는 현상\n")
  cat("   - 직관: Robust SE는 보수적이어야 함 (더 커야 함)\n")
  cat("   - 실제: finegray 데이터에서는 robust SE가 더 작을 수 있음\n")
  cat("   - 이유: Subject-level score aggregation으로 noise 상쇄\n\n")
  
  cat("2. 🎯 CRR vs Fine-Gray 차이\n")
  cat("   - CRR: 직접 구현, Breslow tie 처리\n")
  cat("   - Fine-Gray: Counting process 변환, Efron tie 처리\n")
  cat("   - 계수 차이는 tie 처리 방식과 구현 차이에서 발생\n\n")
  
  cat("3. 🔧 Cluster 보정의 중요성\n")
  cat("   - finegray는 한 subject가 여러 row로 분할됨\n")
  cat("   - cluster(id) 없으면 분산이 과소추정될 수 있음\n")
  cat("   - cluster(id) 사용이 통계적으로 올바른 접근\n\n")
  
  cat("4. ⚖️ Tie 처리 방식의 영향\n")
  cat("   - Breslow (CRR): 단순, tie 많을 때 부정확\n")
  cat("   - Efron (Fine-Gray): 복잡, tie 있어도 정확\n\n")
  
  # 실제 수치로 현상 확인
  crr_se <- se_data$`cmprsk::crr`[1]  # age 계수의 SE
  robust_se <- se_data$`finegray + coxph (Robust SE, no cluster)`[1]
  cluster_se <- se_data$`finegray + coxph (Robust SE + cluster)`[1]
  
  cat("📊 실제 관찰된 현상 (age 계수 기준):\n")
  cat(sprintf("   - CRR SE: %.4f\n", crr_se))
  cat(sprintf("   - Robust SE (no cluster): %.4f (%.1f%%)\n", 
              robust_se, (robust_se/crr_se - 1) * 100))
  cat(sprintf("   - Robust SE (cluster): %.4f (%.1f%%)\n", 
              cluster_se, (cluster_se/crr_se - 1) * 100))
  
  if(robust_se < crr_se) {
    cat("   ✅ Robust SE가 CRR SE보다 작음 - 정상적 현상!\n")
  }
  
  cat("\n💡 실무 권장사항:\n")
  cat("   1. 한 가지 방법을 선택해서 일관되게 사용\n")
  cat("   2. finegray 사용 시 cluster(id) 포함 권장\n")
  cat("   3. Robust SE가 작아져도 우려할 필요 없음\n")
  cat("   4. 중요한 것은 통계적 가정의 충족\n\n")
}

# =============================================================================
# 5. 시뮬레이션 실행 함수 (병렬 처리 + True Value 기반)
# =============================================================================

sample_from_large_data <- function(large_data, sample_size = 300) {
  # 대용량 데이터에서 샘플링
  sample_indices <- sample(1:nrow(large_data), sample_size, replace = FALSE)
  sampled_data <- large_data[sample_indices, ]
  
  # ID 재배정
  sampled_data$id <- 1:nrow(sampled_data)
  
  return(sampled_data)
}

run_educational_simulation <- function(true_values_list, 
                                     n_sim = 50, 
                                     sample_size = 300,
                                     tie_method = "efron") {
  
  cat("🔬 교육용 시뮬레이션 시작\n")
  cat(sprintf("   - 시뮬레이션 횟수: %d\n", n_sim))
  cat(sprintf("   - 샘플 크기: %d\n", sample_size))
  cat(sprintf("   - 사용할 코어: %d개\n", n_cores))
  cat("   - True value 기반 샘플링 방식\n\n")
  
  # 병렬 처리를 위한 함수
  simulate_single <- function(i) {
    # 대용량 데이터에서 샘플링
    sampled_data <- sample_from_large_data(true_values_list$large_data, sample_size)
    
    # 모델 적합
    results <- fit_educational_models(sampled_data, tie_method = tie_method)
    
    # 결과 정리
    sim_result <- list()
    for(method in names(results)) {
      method_result <- results[[method]]
      sim_result[[method]] <- data.frame(
        sim_id = i,
        parameter = names(method_result$coef),
        true_value = true_values_list$values[names(method_result$coef)],
        coef = method_result$coef,
        se = method_result$se,
        bias = method_result$coef - true_values_list$values[names(method_result$coef)],
        stringsAsFactors = FALSE
      )
    }
    
    return(sim_result)
  }
  
  # 병렬 실행
  cat("병렬 시뮬레이션 실행 중...\n")
  sim_list <- mclapply(1:n_sim, simulate_single, mc.cores = n_cores)
  
  # 결과 통합
  sim_results <- list()
  for(method in names(sim_list[[1]])) {
    sim_results[[method]] <- do.call(rbind, lapply(sim_list, function(x) x[[method]]))
  }
  
  cat("✅ 시뮬레이션 완료!\n\n")
  
  return(sim_results)
}

# =============================================================================
# 6. 시뮬레이션 결과 분석 함수
# =============================================================================

analyze_simulation_results <- function(sim_results, true_values) {
  
  cat("📈 시뮬레이션 결과 분석\n")
  cat("=" %+% paste(rep("=", 40), collapse="") %+% "\n\n")
  
  summary_stats <- list()
  
  for(method in names(sim_results)) {
    method_data <- sim_results[[method]]
    
    method_summary <- method_data %>%
      group_by(parameter) %>%
      summarise(
        true_val = first(true_value),
        n_sims = n(),
        mean_est = mean(coef, na.rm = TRUE),
        bias = mean(bias, na.rm = TRUE),
        rmse = sqrt(mean(bias^2, na.rm = TRUE)),
        mean_se = mean(se, na.rm = TRUE),
        emp_se = sd(coef, na.rm = TRUE),
        coverage = mean(abs(bias) <= 1.96 * se, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(method = method)
    
    summary_stats[[method]] <- method_summary
  }
  
  # 전체 요약 통계
  all_stats <- do.call(rbind, summary_stats)
  
  # 결과 출력
  cat("🎯 Bias 비교:\n")
  bias_table <- all_stats %>%
    select(method, parameter, bias) %>%
    pivot_wider(names_from = method, values_from = bias)
  print(kable(bias_table, digits = 6))
  
  cat("\n📊 RMSE 비교:\n")
  rmse_table <- all_stats %>%
    select(method, parameter, rmse) %>%
    pivot_wider(names_from = method, values_from = rmse)
  print(kable(rmse_table, digits = 6))
  
  cat("\n📏 Standard Error 비교:\n")
  se_table <- all_stats %>%
    select(method, parameter, mean_se) %>%
    pivot_wider(names_from = method, values_from = mean_se)
  print(kable(se_table, digits = 6))
  
  cat("\n🎯 Coverage Probability (목표: 0.95):\n")
  coverage_table <- all_stats %>%
    select(method, parameter, coverage) %>%
    pivot_wider(names_from = method, values_from = coverage)
  print(kable(coverage_table, digits = 3))
  
  # SE 비율 분석 (CRR 대비)
  cat("\n⚖️  Standard Error 비율 (CRR 대비):\n")
  crr_se <- all_stats %>% 
    filter(method == "cmprsk::crr") %>% 
    select(parameter, mean_se)
  
  se_ratios <- all_stats %>%
    left_join(crr_se, by = "parameter", suffix = c("", "_crr")) %>%
    mutate(se_ratio = mean_se / mean_se_crr) %>%
    select(method, parameter, se_ratio) %>%
    pivot_wider(names_from = method, values_from = se_ratio)
  
  print(kable(se_ratios, digits = 3))
  
  return(all_stats)
}

# =============================================================================
# 7. 메인 실행 부분 (개선된 교육용 워크플로우)
# =============================================================================

main_educational_analysis <- function(beta_age = 0.5, 
                                    beta_sex = -0.3, 
                                    beta_score = 0.2,
                                    large_n = 50000,
                                    sample_size = 300,
                                    n_sim = 100,
                                    tie_method = "efron") {
  
  cat("🎓 개선된 Competing Risks 교육용 분석\n")
  cat("=" %+% paste(rep("=", 60), collapse="") %+% "\n\n")
  
  # STEP 1: True Value 계산
  cat("📊 STEP 1: True Values 계산\n")
  cat("-" %+% paste(rep("-", 40), collapse="") %+% "\n")
  
  true_values_list <- get_true_values(
    beta_age = beta_age,
    beta_sex = beta_sex, 
    beta_score = beta_score,
    large_n = large_n,
    tie_method = tie_method
  )
  
  # STEP 2: 샘플 데이터로 단일 분석
  cat("🔍 STEP 2: 샘플 데이터 단일 분석\n")
  cat("-" %+% paste(rep("-", 40), collapse="") %+% "\n")
  
  # 대용량 데이터에서 샘플링
  sample_data <- sample_from_large_data(true_values_list$large_data, sample_size)
  
  cat(sprintf("샘플 데이터 요약:\n"))
  cat(sprintf("  - 샘플 크기: %d\n", nrow(sample_data)))
  cat(sprintf("  - Primary event: %d (%.1f%%)\n", 
              sum(sample_data$status == 1), mean(sample_data$status == 1) * 100))
  cat(sprintf("  - Competing event: %d (%.1f%%)\n", 
              sum(sample_data$status == 2), mean(sample_data$status == 2) * 100))
  cat(sprintf("  - 검열: %d (%.1f%%)\n", 
              sum(sample_data$status == 0), mean(sample_data$status == 0) * 100))
  
  # Tie 정보
  tie_times <- table(sample_data$time)
  n_ties <- sum(tie_times > 1)
  cat(sprintf("  - Tie가 있는 시점: %d개\n\n", n_ties))
  
  # 모델 적합
  results <- fit_educational_models(sample_data, tie_method = tie_method)
  
  # 결과 비교 (True values와 함께)
  comparison <- compare_results_with_truth(results, true_values_list$values)
  
  # 교육적 해석
  educational_interpretation_with_truth(comparison$data, true_values_list$values)
  
  # STEP 3: 시뮬레이션 분석
  cat("🔬 STEP 3: 시뮬레이션 분석\n")
  cat("-" %+% paste(rep("-", 40), collapse="") %+% "\n")
  
  sim_results <- run_educational_simulation(
    true_values_list = true_values_list,
    n_sim = n_sim,
    sample_size = sample_size,
    tie_method = tie_method
  )
  
  # 시뮬레이션 결과 분석
  sim_analysis <- analyze_simulation_results(sim_results, true_values_list$values)
  
  # STEP 4: 최종 요약
  cat("🎯 STEP 4: 최종 교육적 요약\n")
  cat("-" %+% paste(rep("-", 40), collapse="") %+% "\n")
  
  final_educational_summary(sim_analysis, true_values_list$values)
  
  return(list(
    true_values = true_values_list,
    sample_data = sample_data,
    single_results = results,
    single_comparison = comparison,
    sim_results = sim_results,
    sim_analysis = sim_analysis
  ))
}

# =============================================================================
# 8. 보조 함수들 (True value 기반)
# =============================================================================

compare_results_with_truth <- function(results, true_values) {
  
  # 결과를 data.frame으로 정리
  comparison_df <- data.frame()
  
  for(method_name in names(results)) {
    method_result <- results[[method_name]]
    
    for(i in 1:length(method_result$coef)) {
      param_name <- names(method_result$coef)[i]
      
      comparison_df <- rbind(comparison_df, data.frame(
        Method = method_result$method,
        Parameter = param_name,
        True_Value = true_values[param_name],
        Coefficient = method_result$coef[i],
        SE = method_result$se[i],
        Bias = method_result$coef[i] - true_values[param_name],
        Tie_Method = method_result$tie_method,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 결과 출력
  cat("=== True Values vs 추정치 비교 ===\n")
  truth_table <- comparison_df %>%
    select(Method, Parameter, True_Value, Coefficient, Bias) %>%
    pivot_wider(names_from = Method, values_from = c(Coefficient, Bias))
  
  print(kable(truth_table, digits = 6))
  
  cat("\n=== 표준오차 비교 ===\n")
  se_comparison <- comparison_df %>%
    select(Method, Parameter, SE) %>%
    pivot_wider(names_from = Method, values_from = SE)
  
  print(kable(se_comparison, digits = 6))
  
  return(list(data = comparison_df))
}

educational_interpretation_with_truth <- function(comparison_data, true_values) {
  
  cat("\n" %+% paste(rep("=", 70), collapse="") %+% "\n")
  cat("📚 True Value 기반 교육적 해석\n")
  cat(paste(rep("=", 70), collapse="") %+% "\n\n")
  
  cat("🎯 True Values:\n")
  for(i in 1:length(true_values)) {
    cat(sprintf("   %s: %.6f\n", names(true_values)[i], true_values[i]))
  }
  cat("\n")
  
  # 각 방법별 bias 계산
  age_data <- comparison_data %>% filter(Parameter == "age")
  
  cat("🔍 AGE 계수에 대한 분석:\n")
  for(i in 1:nrow(age_data)) {
    method <- age_data$Method[i]
    bias <- age_data$Bias[i]
    se <- age_data$SE[i]
    
    cat(sprintf("   %s:\n", method))
    cat(sprintf("     - Bias: %.6f\n", bias))
    cat(sprintf("     - SE: %.6f\n", se))
    cat(sprintf("     - |Bias|/SE: %.3f\n", abs(bias)/se))
  }
  
  cat("\n💡 핵심 학습 포인트:\n")
  cat("   1. True value와의 차이로 실제 성능 평가 가능\n")
  cat("   2. Bias와 SE의 trade-off 관계 관찰\n")
  cat("   3. 방법별 강점과 약점 명확히 파악\n")
  cat("   4. 단일 분석의 한계와 시뮬레이션의 필요성\n\n")
}

final_educational_summary <- function(sim_analysis, true_values) {
  
  cat("🏆 최종 교육적 요약\n")
  cat("=" %+% paste(rep("=", 50), collapse="") %+% "\n\n")
  
  cat("📋 검증된 주요 현상들:\n\n")
  
  # SE 비율 분석
  age_stats <- sim_analysis %>% filter(parameter == "age")
  crr_se <- age_stats$mean_se[age_stats$method == "cmprsk::crr"]
  
  robust_methods <- age_stats %>% 
    filter(grepl("Robust", method)) %>%
    arrange(mean_se)
  
  cat("1. 📉 Robust SE 현상 (AGE 계수 기준):\n")
  for(i in 1:nrow(robust_methods)) {
    method <- robust_methods$method[i]
    se <- robust_methods$mean_se[i]
    ratio <- se / crr_se
    
    cat(sprintf("   %s: %.6f (CRR 대비 %.1f%%)\n", 
                method, se, (ratio - 1) * 100))
  }
  
  cat("\n2. 📊 Coverage 성능:\n")
  coverage_stats <- sim_analysis %>% 
    filter(parameter == "age") %>%
    select(method, coverage) %>%
    arrange(desc(coverage))
  
  for(i in 1:nrow(coverage_stats)) {
    method <- coverage_stats$method[i]
    cov <- coverage_stats$coverage[i]
    status <- ifelse(abs(cov - 0.95) < 0.02, "✅", "⚠️")
    
    cat(sprintf("   %s %s: %.3f\n", status, method, cov))
  }
  
  cat("\n3. 🎯 RMSE 성능:\n")
  rmse_stats <- sim_analysis %>% 
    filter(parameter == "age") %>%
    select(method, rmse) %>%
    arrange(rmse)
  
  for(i in 1:nrow(rmse_stats)) {
    method <- rmse_stats$method[i]
    rmse <- rmse_stats$rmse[i]
    
    cat(sprintf("   %d. %s: %.6f\n", i, method, rmse))
  }
  
  cat("\n🎓 교육적 결론:\n")
  cat("   ✅ Robust SE가 작아지는 현상은 정상적임\n")
  cat("   ✅ 방법별로 고유한 특성과 장단점 존재\n")
  cat("   ✅ 시뮬레이션을 통한 체계적 평가의 중요성\n")
  cat("   ✅ True value 기반 학습의 교육적 효과\n\n")
}

# =============================================================================
# 실행
# =============================================================================

# 메인 분석 실행 (개선된 버전)
cat("🚀 교육용 Competing Risks 분석 시작\n\n")

educational_results <- main_educational_analysis(
  beta_age = 0.5,
  beta_sex = -0.3, 
  beta_score = 0.2,
  large_n = 30000,    # 대용량 데이터 크기 (원래대로)
  sample_size = 400,  # 샘플 크기
  n_sim = 200,        # 시뮬레이션 횟수 (원래대로)
  tie_method = "efron"
)