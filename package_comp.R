library(survival)
library(cmprsk)
library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)
library(doParallel)
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
  
  data.frame(time = obs_time, status = status, age = age, sex = sex, score = score)
}

# Model fitting
fit_models <- function(data, include_interaction = FALSE) {
  tryCatch({
    status_factor <- factor(data$status, levels=c(0,1,2), labels=c("cens","primary","competing"))
    
    if (include_interaction) {
      fg_formula <- Surv(time, status_factor) ~ age + sex + score + age:sex
    } else {
      fg_formula <- Surv(time, status_factor) ~ age + sex + score
    }
    
    fg_data <- finegray(fg_formula, data = data, etype = "primary")
    fg_data$id <- 1:nrow(fg_data)
    
    if (include_interaction) {
      fg_fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + score + age:sex + cluster(id),
                      data = fg_data, weights = fgwt)
      cov_matrix <- model.matrix(~ age + sex + score + age:sex, data = data)[, -1]
    } else {
      fg_fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + score + cluster(id),
                      data = fg_data, weights = fgwt)
      cov_matrix <- model.matrix(~ age + sex + score, data = data)[, -1]
    }
    
    crr_fit <- crr(ftime = data$time, fstatus = data$status, cov1 = cov_matrix,
                   failcode = 1, cencode = 0)
    
    n_coef <- ifelse(include_interaction, 4, 3)
    list(
      fg = list(coef = coef(fg_fit)[1:n_coef], se = sqrt(diag(fg_fit$var))[1:n_coef]),
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
                            model_interaction = FALSE) {
  # Generate large data WITHOUT ties
  large_data <- generate_data(n = 50000, 
                              beta1 = true_beta1, beta2 = true_beta2, 
                              beta3 = true_beta3, beta4 = true_beta4,
                              apply_ties = FALSE)
  
  # Fit model according to specification
  fits <- fit_models(large_data, include_interaction = model_interaction)
  
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
                           n_sim = 500,
                           n = 500,
                           true_beta1 = 0.5, true_beta2 = -0.3, true_beta3 = 0.2, true_beta4 = 0,
                           model_interaction = FALSE,
                           data_has_ties = FALSE,
                           tie_precision = 0.25) {
  
  # Get true values (always without ties)
  true_values <- get_true_values(true_beta1, true_beta2, true_beta3, true_beta4, model_interaction)
  
  n_cores <- 2
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    library(survival)
    library(cmprsk)
  })
  clusterExport(cl, c("generate_data", "fit_models"), envir = environment())
  
  cat(sprintf("\nRunning: %s\n", scenario_name))
  cat(sprintf("True interaction: %.2f, Model includes interaction: %s, Data has ties: %s\n", 
              true_beta4, model_interaction, data_has_ties))
  
  results <- parLapply(cl, 1:n_sim, function(i) {
    set.seed(12345 + i)
    
    data <- generate_data(n = n, 
                          beta1 = true_beta1, beta2 = true_beta2,
                          beta3 = true_beta3, beta4 = true_beta4,
                          apply_ties = data_has_ties, tie_precision = tie_precision)
    
    fits <- fit_models(data, include_interaction = model_interaction)
    
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
  })
  
  stopCluster(cl)
  
  results <- do.call(rbind, results[!sapply(results, is.null)])
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







true_values <- get_true_values(
  true_beta1 = 0.5, 
  true_beta2 = -0.3, 
  true_beta3 = 0.2, 
  true_beta4 = 0,      # interaction 없음
  model_interaction = FALSE  # 모델도 interaction 없이
)

# 시뮬레이션 실행 (매일 측정으로 tie 발생)
result <- run_simulation(
  scenario_name = "Hospital_Daily_Check",
  n_sim = 500,
  n = 500,
  true_beta1 = 0.5,
  true_beta2 = -0.3, 
  true_beta3 = 0.2,
  true_beta4 = 0,           # 실제 interaction 없음
  model_interaction = FALSE,  # 모델도 interaction 없이
  data_has_ties = TRUE,      # tie 있음
  tie_precision = 1          # 일 단위 측정
)