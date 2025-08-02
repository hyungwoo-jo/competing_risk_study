library(survival)
library(cmprsk)
library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)
# library(doParallel) # No longer needed
library(data.table)
set.seed(12345)

# Flexible data generator
generate_data <- function(n = 500,
                          beta1 = 0.5, beta2 = -0.3, beta3 = 0.2, beta4 = 0,
                          apply_ties = FALSE, tie_precision = 0.25) {
  age <- scale(rnorm(n, 60, 10))[,1]
  sex <- rbinom(n, 1, 0.5)
  score <- scale(rnorm(n, 50, 15))[,1]
  
  lp1 <- beta1*age + beta2*sex + beta3*score + beta4*age*sex
  lp2 <- 0.3*age - 0.2*sex + 0.15*score
  
  time1 <- rexp(n, exp(lp1))
  time2 <- rexp(n, exp(lp2))
  cens_time <- rexp(n, 0.08)
  
  obs_time <- pmin(time1, time2, cens_time)
  status <- ifelse(time1 <= pmin(time2, cens_time), 1,
                   ifelse(time2 <= pmin(time1, cens_time), 2, 0))
  
  if (apply_ties) {
    obs_time <- ceiling(obs_time / tie_precision) * tie_precision
  }
  
  data.frame(id = 1:n, time = obs_time, status = status, age = age, sex = sex, score = score)
}

# Model fitting
fit_models <- function(data, include_interaction = FALSE, use_robust_se = TRUE, tie_method = "efron") {
  tryCatch({
    status_factor <- factor(data$status, levels=c(0,1,2), labels=c("cens","primary","competing"))
    
    fg_formula <- Surv(time, status_factor) ~ .
    
    
    fg_data <- finegray(fg_formula, data = data, etype = "primary")
    if (include_interaction) {
      cox_formula_str <- "Surv(fgstart, fgstop, fgstatus) ~ age + sex + score + age:sex"
      cov_matrix <- model.matrix(~ age + sex + score + age:sex, data = data)[, -1]
    } else {
      cox_formula_str <- "Surv(fgstart, fgstop, fgstatus) ~ age + sex + score"
      cov_matrix <- model.matrix(~ age + sex + score, data = data)[, -1]
    }
    
    
    
    fg_fit <- coxph(as.formula(cox_formula_str),
                    data = fg_data,
                    weights = fgwt,
                    ties = tie_method)
    
    if (use_robust_se) {
      fg_summary <- summary(fg_fit)
      fg_se <- fg_summary$coefficients[, "robust se"]
    } else {
      fg_summary <- summary(fg_fit)
      fg_se <- fg_summary$coefficients[, "se(coef)"]
    }
    crr_fit <- crr(ftime = data$time, fstatus = data$status, cov1 = cov_matrix,
                   failcode = 1, cencode = 0)
    
    n_coef <- ifelse(include_interaction, 4, 3)
    list(
      fg = list(coef = coef(fg_fit)[1:n_coef], se = fg_se[1:n_coef]),
      crr = list(coef = crr_fit$coef[1:n_coef], se = sqrt(diag(crr_fit$var))[1:n_coef])
    )
  }, error = function(e) {
    n_coef <- ifelse(include_interaction, 4, 3)
    return(list(
      fg = list(coef = rep(NA, n_coef), se = rep(NA, n_coef)),
      crr = list(coef = rep(NA, n_coef), se = rep(NA, n_coef))
    ))
  })
}

# Get true values (always from no-tie data)
get_true_values <- function(true_beta1 = 0.5, true_beta2 = -0.3, true_beta3 = 0.2, true_beta4 = 0,
                            model_interaction = FALSE, use_robust_se = TRUE, tie_method = "efron") {
  large_data <- generate_data(n = 50000,
                              beta1 = true_beta1, beta2 = true_beta2,
                              beta3 = true_beta3, beta4 = true_beta4,
                              apply_ties = FALSE)
  
  fits <- fit_models(large_data, include_interaction = model_interaction, use_robust_se = use_robust_se, tie_method = tie_method)
  
  param_names <- if(model_interaction) c("age", "sex", "score", "age:sex") else c("age", "sex", "score")
  
  cat("True Values (No Ties, Large N):\n")
  for(i in 1:length(fits$crr$coef)) {
    cat(sprintf("%s - FG: %.4f, CRR: %.4f\n",
                param_names[i], fits$fg$coef[i], fits$crr$coef[i]))
  }
  cat("\n")
  
  return(fits$crr$coef)
}

# Flexible simulation runner
run_simulation <- function(scenario_name,
                           true_values, # true_values를 인자로 받음
                           n_sim = 500,
                           n = 500,
                           true_beta1 = 0.5, true_beta2 = -0.3, true_beta3 = 0.2, true_beta4 = 0,
                           model_interaction = FALSE,
                           data_has_ties = FALSE,
                           tie_precision = 0.25,
                           use_robust_se_coxph = TRUE,
                           tie_method = "efron") {
  
  n_cores <- 16
  
  cat(sprintf("\nRunning: %s\n", scenario_name))
  cat(sprintf("True interaction: %.2f, Model includes interaction: %s, Data has ties: %s, Coxph uses robust SE: %s, Tie method: %s\n",
              true_beta4, model_interaction, data_has_ties, use_robust_se_coxph, tie_method))
  
  # Using mclapply for parallel processing
  results_list <- mclapply(1:n_sim, function(i) {
    set.seed(12345 + i)
    
    data <- generate_data(n = n,
                          beta1 = true_beta1, beta2 = true_beta2,
                          beta3 = true_beta3, beta4 = true_beta4,
                          apply_ties = data_has_ties, tie_precision = tie_precision)
    fits <- fit_models(data, include_interaction = model_interaction, use_robust_se = use_robust_se_coxph, tie_method = tie_method)
    
    if (!any(is.na(fits$fg$coef)) && !any(is.na(fits$crr$coef))) {
      param_names <- if(model_interaction) c("age", "sex", "score", "age:sex") else c("age", "sex", "score")
      sim_results <- data.frame()
      for (j in seq_along(param_names)) {
        sim_results <- rbind(sim_results, data.frame(
          sim = i,
          parameter = param_names[j],
          true_value = true_values[j],
          fg_est = fits$fg$coef[j],
          fg_se = fits$fg$se[j],
          crr_est = fits$crr$coef[j],
          crr_se = fits$crr$se[j]
        ))
      }
      return(sim_results)
    }
    return(NULL)
  }, mc.cores = n_cores)
  
  results <- rbindlist(results_list[!sapply(results_list, is.null)])
  
  results$scenario <- scenario_name
  
  # Summary
  summary_stats <- results %>%
    group_by(parameter) %>%
    summarise(
      true_val = first(true_value),
      fg_bias = mean(fg_est - true_value),
      fg_rmse = sqrt(mean((fg_est - true_value)^2)),
      fg_coverage = mean(abs(fg_est - true_value) <= 1.96 * fg_se),
      crr_bias = mean(crr_est - true_value),
      crr_rmse = sqrt(mean((crr_est - true_value)^2)),
      crr_coverage = mean(abs(crr_est - true_value) <= 1.96 * crr_se),
      .groups = "drop"
    )
  
  summary_stats$scenario <- scenario_name
  
  print(summary_stats %>%
          select(parameter, fg_bias, crr_bias, fg_rmse, crr_rmse),
        digits = 3)
  
  return(list(results = results, summary = summary_stats))
}

# true_values를 외부에서 한 번만 계산하거나 직접 지정
true_values <- get_true_values(
  true_beta1 = 0.5,
  true_beta2 = -0.3,
  true_beta3 = 0.2,
  true_beta4 = 0,        # interaction 없음
  model_interaction = F   # 모델도 interaction 없이
)
print(true_values)
# true_values <- c(0.2216, -0.1006, 0.0601)
# fwrite(as.data.frame(as.list(true_values)), 'tv_interaction.csv')
# true_values<-fread('tv_interaction.csv')
# print('go')
# true_values
true_values<-c(0.2235835, -0.09834393, 0.05931059)
true_values <- c(0.2216, -0.1006, 0.0601)
# 시뮬레이션 실행 시 true_values를 인자로 전달
result <- run_simulation(
  scenario_name = "Hospital_Daily_Check",
  true_values = true_values, # 계산된 true_values 전달
  n_sim = 300,
  n = 500,
  true_beta1 = 0.5,
  true_beta2 = -0.3,
  true_beta3 = 0.2,
  true_beta4 = 0,           # 실제 interaction 없음
  model_interaction = FALSE,  # 모델도 interaction 없이
  data_has_ties = TRUE,       # tie 있음
  tie_precision = 1,          # 일 단위 측정
  use_robust_se_coxph = F, # FALSE: non-robust SE
  tie_method = "breslow"      # Added tie_method argument, default to "efron"
)


fwrite(result$results, 'rr_breslow.csv')
fwrite(result$summary, 'rs_breslow.csv')
print(result$summary)
