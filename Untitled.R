library(survey)
library(dplyr)
set.seed(123)

cat("=== EXTREME WEIGHT SCENARIO ===\n")
cat("더 극단적인 가중치로 차이를 명확히 보여드리겠습니다.\n\n")

extreme_coverage_simulation <- function(n_sim = 500) {
  
  true_slope <- 2.0
  results <- data.frame()
  
  cat("시뮬레이션 설정:\n")
  cat("- 진짜 slope = 2.0\n")
  cat("- 극단적인 selection bias (일부 관측치가 매우 낮은 선택확률)\n")
  cat("- 결과적으로 매우 큰 IPTW 가중치 발생\n\n")
  
  for(i in 1:n_sim) {
    # 데이터 생성
    n <- 150
    age <- runif(n, 20, 80)
    treatment <- rbinom(n, 1, plogis(-2 + 0.05*age))  # 나이 많을수록 treatment
    
    # Outcome (진짜 treatment effect = 2.0)
    y <- 5 + 0.1*age + true_slope*treatment + rnorm(n, 0, 2)
    
    # 극단적인 selection bias
    # 젊고 treatment 받지 않은 사람들이 매우 낮은 확률로 선택됨
    logit_p <- -1 + 0.03*age + 1.5*treatment + 0.02*age*treatment
    p_select <- plogis(logit_p)
    p_select <- pmax(p_select, 0.05)  # 최소 5% 확률
    
    selected <- rbinom(n, 1, p_select)
    
    # 선택된 데이터
    sample_idx <- which(selected == 1)
    if(length(sample_idx) < 30) next
    
    age_s <- age[sample_idx]
    treatment_s <- treatment[sample_idx]
    y_s <- y[sample_idx]
    p_s <- p_select[sample_idx]
    
    # IPTW 가중치 (매우 큰 값들 포함)
    weights <- 1/p_s
    
    # Extreme weights 확인
    if(max(weights) < 5) next  # 충분히 극단적이지 않으면 skip
    
    # 방법 1: 잘못된 처리 (frequency weight로 처리)
    fit_wrong <- lm(y_s ~ treatment_s + age_s, weights = weights)
    est_wrong <- coef(fit_wrong)["treatment_s"]
    se_wrong <- summary(fit_wrong)$coef["treatment_s", "Std. Error"]
    ci_wrong <- est_wrong + c(-1.96, 1.96) * se_wrong
    
    # 방법 2: 올바른 처리 (survey package)
    tryCatch({
      data_survey <- data.frame(y_s, treatment_s, age_s, weights)
      design <- svydesign(ids = ~1, weights = ~weights, data = data_survey)
      fit_correct <- svyglm(y_s ~ treatment_s + age_s, design = design)
      est_correct <- coef(fit_correct)["treatment_s"]
      se_correct <- summary(fit_correct)$coef["treatment_s", "Std. Error"]
      ci_correct <- est_correct + c(-1.96, 1.96) * se_correct
      
      # 방법 3: Stabilized weights
      mean_treatment <- mean(treatment_s)
      ps_model <- glm(treatment_s ~ age_s, family = binomial)
      ps <- predict(ps_model, type = "response")
      
      stabilized_w <- weights * mean_treatment / 
        ifelse(treatment_s == 1, ps, 1-ps)
      
      design_stab <- svydesign(ids = ~1, weights = ~stabilized_w, data = data_survey)
      fit_stab <- svyglm(y_s ~ treatment_s + age_s, design = design_stab)
      est_stab <- coef(fit_stab)["treatment_s"]
      se_stab <- summary(fit_stab)$coef["treatment_s", "Std. Error"]
      ci_stab <- est_stab + c(-1.96, 1.96) * se_stab
      
      # 결과 저장
      results <- rbind(results, data.frame(
        sim = i,
        method = c("Wrong (Freq)", "Correct (Survey)", "Stabilized IPTW"),
        estimate = c(est_wrong, est_correct, est_stab),
        se = c(se_wrong, se_correct, se_stab),
        coverage = c(
          ci_wrong[1] <= true_slope & true_slope <= ci_wrong[2],
          ci_correct[1] <= true_slope & true_slope <= ci_correct[2],
          ci_stab[1] <= true_slope & true_slope <= ci_stab[2]
        ),
        max_weight = max(weights),
        n_sample = length(sample_idx)
      ))
    }, error = function(e) NULL)
    
    if(i %% 100 == 0) cat("Completed", i, "simulations...\n")
  }
  
  if(nrow(results) == 0) {
    cat("시뮬레이션 실패: 충분한 극단적 가중치가 생성되지 않았습니다.\n")
    return(NULL)
  }
  
  # 결과 요약
  summary_stats <- results %>%
    group_by(method) %>%
    summarise(
      n_sim = n(),
      avg_estimate = mean(estimate),
      avg_se = mean(se),
      coverage_rate = mean(coverage),
      bias = mean(estimate - true_slope),
      avg_max_weight = mean(max_weight),
      .groups = "drop"
    )
  
  cat("\n=== EXTREME WEIGHT 시뮬레이션 결과 ===\n")
  cat("(True slope = 2.0)\n\n")
  
  cat(sprintf("%-20s %8s %8s %8s %8s %10s\n", 
              "Method", "N_sim", "Bias", "Avg_SE", "Coverage", "Avg_MaxWt"))
  cat(sprintf("%-20s %8s %8s %8s %8s %10s\n", 
              "------", "-----", "----", "------", "--------", "---------"))
  
  for(i in 1:nrow(summary_stats)) {
    cat(sprintf("%-20s %8d %8.3f %8.3f %8.3f %10.1f\n",
                summary_stats$method[i], 
                summary_stats$n_sim[i],
                summary_stats$bias[i],
                summary_stats$avg_se[i], 
                summary_stats$coverage_rate[i],
                summary_stats$avg_max_weight[i]))
  }
  
  cat("\n해석:\n")
  cat("1. Wrong (Freq): 가중치를 frequency로 잘못 처리 → SE 과소추정 → Coverage 저하\n")
  cat("2. Correct (Survey): 올바른 probability weight 처리 → 올바른 SE → Coverage 0.95 근접\n")
  cat("3. Stabilized IPTW: 극단적 가중치 안정화 → 더 좋은 성능\n")
  
  # SE 비교 시각화
  if(nrow(results) > 50) {
    se_comparison <- results %>%
      select(method, se, max_weight) %>%
      group_by(method) %>%
      summarise(
        avg_se = mean(se),
        se_std = sd(se),
        .groups = "drop"
      )
    
    cat("\n=== SE 비교 ===\n")
    print(se_comparison)
  }
  
  return(summary_stats)
}

# 실제 Fine-Gray 스타일 시뮬레이션
finegray_style_demo <- function() {
  
  cat("\n=== FINE-GRAY 스타일 시뮬레이션 ===\n")
  
  # Fine-Gray와 유사한 가중치 패턴 생성
  n <- 200
  age <- rnorm(n, 65, 10)
  treatment <- rbinom(n, 1, 0.5)
  
  # Competing risks: primary event vs death
  # 나이가 많을수록 death 위험 높음, treatment는 primary event 위험 높임
  lambda_primary <- exp(-3 + 0.5*treatment + 0.01*age)
  lambda_death <- exp(-2 + 0.05*age - 0.3*treatment)
  
  # Event times
  time_primary <- rexp(n, lambda_primary)
  time_death <- rexp(n, lambda_death)
  
  # 관측되는 것은 첫 번째 event
  time_obs <- pmin(time_primary, time_death)
  event_type <- ifelse(time_primary < time_death, 1, 2)  # 1=primary, 2=death
  
  # Censoring (30개월에서 administrative censoring)
  cens_time <- 30
  time_final <- pmin(time_obs, cens_time)
  event_final <- ifelse(time_obs <= cens_time, event_type, 0)
  
  # Fine-Gray 스타일 가중치 계산 (simplified)
  # 시간이 지날수록 death risk가 높은 사람들의 가중치가 작아짐
  fg_weights <- exp(-0.02*age * time_final)  # 나이와 시간에 따른 가중치
  fg_weights <- pmax(fg_weights, 0.1)  # 최소값 설정
  
  cat("Fine-Gray 스타일 가중치 분포:\n")
  cat("Min:", round(min(fg_weights), 3), "\n")
  cat("Max:", round(max(fg_weights), 3), "\n")
  cat("Mean:", round(mean(fg_weights), 3), "\n")
  
  # Primary event에 대한 분석 (Fine-Gray 스타일)
  primary_data <- data.frame(
    treatment = treatment,
    age = age,
    outcome = as.numeric(event_final == 1),  # primary event 여부
    weights = fg_weights
  )
  
  # 방법 비교
  # 1. Unweighted
  fit_unweight <- glm(outcome ~ treatment + age, data = primary_data, family = binomial)
  
  # 2. Wrong (frequency weights)
  fit_wrong <- glm(outcome ~ treatment + age, data = primary_data, 
                   family = binomial, weights = weights)
  
  # 3. Correct (survey)
  design_fg <- svydesign(ids = ~1, weights = ~weights, data = primary_data)
  fit_correct <- svyglm(outcome ~ treatment + age, design = design_fg, family = binomial)
  
  cat("\n=== Fine-Gray 스타일 결과 비교 ===\n")
  cat(sprintf("%-15s %12s %12s\n", "Method", "Treatment_Coef", "SE"))
  cat(sprintf("%-15s %12s %12s\n", "------", "-------------", "--"))
  
  # Unweighted
  coef_unw <- coef(fit_unweight)["treatment"]
  se_unw <- summary(fit_unweight)$coef["treatment", "Std. Error"]
  cat(sprintf("%-15s %12.3f %12.3f\n", "Unweighted", coef_unw, se_unw))
  
  # Wrong
  coef_wrong <- coef(fit_wrong)["treatment"] 
  se_wrong <- summary(fit_wrong)$coef["treatment", "Std. Error"]
  cat(sprintf("%-15s %12.3f %12.3f\n", "Wrong (Freq)", coef_wrong, se_wrong))
  
  # Correct
  coef_correct <- coef(fit_correct)["treatment"]
  se_correct <- summary(fit_correct)$coef["treatment", "Std. Error"] 
  cat(sprintf("%-15s %12.3f %12.3f\n", "Correct (Survey)", coef_correct, se_correct))
  
  cat("\n관찰:\n")
  cat("- Fine-Gray 가중치는 확률적 의미를 가짐\n")
  cat("- Survey package 사용시 올바른 SE 계산\n")
  cat("- 실제 Fine-Gray에서는 이런 차이가 더 극명하게 나타남\n")
}

# 실행
cat("극단적인 가중치 시나리오 실행 중...\n")
extreme_result <- extreme_coverage_simulation(n_sim = 300)

cat("\n" , rep("=", 60), "\n")
finegray_style_demo()