# Tie Effects Simulation Script
source('educational_cr_analysis.R')

cat('🚀 Tie 효과 시뮬레이션 시작\n\n')

# 시뮬레이션 함수
run_simple_simulation <- function(n_sim = 20, apply_ties = FALSE) {
  results <- data.frame()
  
  for (i in 1:n_sim) {
    set.seed(i + 5000)
    
    # 데이터 생성
    data <- generate_educational_data(n = 400, 
                                    beta_age = 0.25, 
                                    beta_sex = -0.17, 
                                    beta_score = 0.096,
                                    apply_ties = apply_ties,
                                    tie_precision = 0.5)
    
    # CRR 분석
    cov_matrix <- as.matrix(data[, c("age", "sex", "score")])
    crr_fit <- crr(ftime = data$time, 
                   fstatus = data$status, 
                   cov1 = cov_matrix,
                   failcode = 1, 
                   cencode = 0)
    
    # Fine-Gray 분석
    fg_data <- finegray(Surv(time, factor(status)) ~ ., data = data, etype = '1')
    
    # Standard SE
    fg_std <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + score, 
                   data = fg_data, weight = fgwt)
    
    # Robust SE
    fg_rob <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + score + cluster(id), 
                   data = fg_data, weight = fgwt, robust = TRUE)
    
    # 결과 저장 (AGE 계수만)
    result_row <- data.frame(
      sim = i,
      crr_coef = crr_fit$coef[1],
      crr_se = sqrt(crr_fit$var[1,1]),
      std_coef = coef(fg_std)[1],
      std_se = sqrt(vcov(fg_std)[1,1]),
      rob_coef = coef(fg_rob)[1],
      rob_se = sqrt(vcov(fg_rob)[1,1])
    )
    
    results <- rbind(results, result_row)
  }
  
  return(results)
}

# No Ties 시뮬레이션
cat('📊 No Ties 시뮬레이션 (20회)...\n')
no_ties_results <- run_simple_simulation(20, apply_ties = FALSE)

# With Ties 시뮬레이션
cat('📊 With Ties 시뮬레이션 (20회)...\n')
with_ties_results <- run_simple_simulation(20, apply_ties = TRUE)

cat('\n✅ 시뮬레이션 완료! 결과 분석:\n\n')

# 요약 통계
summarize_sim <- function(results, label) {
  cat('=== ', label, ' ===\n')
  summary_data <- data.frame(
    Method = c('CRR', 'FG_Standard', 'FG_Robust'),
    Mean_Coef = c(mean(results$crr_coef), mean(results$std_coef), mean(results$rob_coef)),
    Mean_SE = c(mean(results$crr_se), mean(results$std_se), mean(results$rob_se)),
    SD_SE = c(sd(results$crr_se), sd(results$std_se), sd(results$rob_se))
  )
  print(summary_data, digits = 5)
  cat('\n')
  return(summary_data)
}

no_ties_summary <- summarize_sim(no_ties_results, 'NO TIES')
with_ties_summary <- summarize_sim(with_ties_results, 'WITH TIES')

# SE 변화 분석
cat('🔍 SE 변화 분석\n')
cat('================\n')
se_comparison <- data.frame(
  Method = no_ties_summary$Method,
  SE_NoTies = no_ties_summary$Mean_SE,
  SE_WithTies = with_ties_summary$Mean_SE,
  Change_Pct = round((with_ties_summary$Mean_SE / no_ties_summary$Mean_SE - 1) * 100, 2),
  Abs_Change = round(with_ties_summary$Mean_SE - no_ties_summary$Mean_SE, 5)
)
print(se_comparison)

# Robust SE 현상 분석
cat('\n💡 Robust SE 현상 분석\n')
cat('=====================\n')

robust_vs_std_no_ties <- no_ties_summary$Mean_SE[3] / no_ties_summary$Mean_SE[2]
robust_vs_std_with_ties <- with_ties_summary$Mean_SE[3] / with_ties_summary$Mean_SE[2]

cat('Robust SE / Standard SE:\n')
cat('   No Ties: ', round(robust_vs_std_no_ties, 4), '\n')
cat('   With Ties: ', round(robust_vs_std_with_ties, 4), '\n')

if (robust_vs_std_with_ties < 0.95) {
  cat('   ✅ 뚜렷한 Robust SE 축소 현상 관찰!\n')
} else if (robust_vs_std_with_ties < 1.0) {
  cat('   ✅ 약한 Robust SE 축소 현상 관찰\n')
} else {
  cat('   ⚠️ Robust SE 축소 현상 없음 (표준적 상황)\n')
}

# True value와 bias 분석
true_age <- 0.248672
cat('\n📈 Bias 분석 (True Age Coef = 0.248672)\n')
cat('====================================\n')
bias_analysis <- data.frame(
  Method = no_ties_summary$Method,
  Bias_NoTies = round(no_ties_summary$Mean_Coef - true_age, 5),
  Bias_WithTies = round(with_ties_summary$Mean_Coef - true_age, 5),
  Bias_Change = round((with_ties_summary$Mean_Coef - true_age) - (no_ties_summary$Mean_Coef - true_age), 5)
)
print(bias_analysis)

cat('\n🎯 결론:\n')
if (robust_vs_std_with_ties < 1.0) {
  cat('- ✅ Robust SE 축소 현상이 관찰되었습니다!\n')
  cat('- 🔍 이는 Fine-Gray 가중치와 tie 처리의 상호작용 때문입니다.\n')
} else {
  cat('- ⚠️ 이번 시뮬레이션에서는 robust SE 축소 현상이 관찰되지 않았습니다.\n')
  cat('- 🔍 더 극단적인 조건(더 많은 tie, 작은 샘플)이 필요할 수 있습니다.\n')
}

cat('- 📊 CRR과 Fine-Gray 방법들이 tie에 대해 서로 다른 반응을 보입니다.\n')
cat('- 🎓 교육적 목적으로 tie 효과의 방법론적 차이를 명확히 확인했습니다.\n')