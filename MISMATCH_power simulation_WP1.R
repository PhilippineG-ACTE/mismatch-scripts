---
title: "Power Simulation for MISMATCH Project-Sample size estimation for WP1"
author: "Philippine Geelhand"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
```{r c}
## Power simulation to determine the target sample size of WP1
## Dyad matching procedure: 2 groups of 6 participants (2 NT, 2 AU, 2 SZ):
# 3 same-neurotype interactions: NT1-NT2, SZ1-SZ2, AU1-AU2
# 4 mixed-neurotype interactions: NT1-AU1, NT1-SZ1, NT2-AU2, NT2-SZ2
# 2 mixed-neurodivergence: AU1-SZ2, AU2-SZ1
# Every participant gets exactly 3 partners -> 9 dyads per block

## Within each dyad, both members serve as Director and Matcher (5 turns each)
## Sample size scales by replicating this block (n_blocks x 6 participants)
## Primary effect of interest: status_match (Match vs Mismatch)
## Random effects (to control for data interdependence): (1|director_id) + (1|dyad_id)

## Required packages
library(lme4)
library(MASS)

#### CONFIGURATION SECTION ####

n_blocks <- 14             # 1 block = 6 participants, 9 dyads, start testing with 14 blocks (84 participants)
turns_per_dyad <- 10       # 5 as Director + 5 as Matcher per dyad
n_simulations <- 600       # Simulations per sample size

## Effect sizes
# Diagnosis main effects (relative to NT reference) - control variables
effect_director_AU <- 0.2
effect_director_SZ <- 0.2
effect_matcher_AU  <- 0.2
effect_matcher_SZ  <- 0.2

effect_turn <- 0.1         # Linear trend across turns within a dyad

## Count data parameters (starting point from original study: Geelhand et al.(2025)
baseline_rate <- 5.25
overdispersion <- 3.68     # Fixed theta derived from Geelhand et al. (2025)


# Geelhand et al. (2025): large effect size (chi-sq = 20.7) but in a simpler between-group design (Director*Matcher)
# Three scenarios testing small, medium and large effects of neurotype (mis)match

scenarios <- list(
  conservative = list(effect_match_bonus = 0.2, sd_director = 0.35, sd_dyad = 0.35),
  moderate     = list(effect_match_bonus = 0.5, sd_director = 0.25, sd_dyad = 0.25),
  optimistic   = list(effect_match_bonus = 0.8, sd_director = 0.17, sd_dyad = 0.17)
)


build_block_dyads <- function(offset) {
  # Local indices: 1=NT1, 2=NT2, 3=AU1, 4=AU2, 5=SZ1, 6=SZ2
  local_diag <- c("NT", "NT", "AU", "AU", "SZ", "SZ")
  dyad_pairs <- list(c(1,2), c(1,3), c(1,5), c(2,4), c(2,6),
                     c(3,4), c(3,6), c(4,5), c(5,6))
  
  list(
    participant_id = offset + 1:6,
    diag = local_diag,
    dyad_pairs = lapply(dyad_pairs, function(p) offset + p)
  )
}

#### SIMULATION FUNCTION ####

simulate_multiPartner_study <- function(n_blocks, turns_per_dyad,
                                        effect_match_bonus,
                                        random_effect_sd_director,
                                        random_effect_sd_dyad) {
  
  half_turns <- turns_per_dyad / 2
  all_rows <- list()
  row_counter <- 0
  dyad_counter <- 0
  participant_diag_lookup <- character(n_blocks * 6)
  
  for (b in 1:n_blocks) {
    offset <- (b - 1) * 6
    block <- build_block_dyads(offset)
    participant_diag_lookup[block$participant_id] <- block$diag
    
    for (pair in block$dyad_pairs) {
      dyad_counter <- dyad_counter + 1
      member_a <- pair[1]
      member_b <- pair[2]
      
      for (t in 1:half_turns) {
        row_counter <- row_counter + 1
        all_rows[[row_counter]] <- data.frame(
          dyad_id = dyad_counter, turn = t,
          director_id = member_a, matcher_id = member_b
        )
      }
      for (t in 1:half_turns) {
        row_counter <- row_counter + 1
        all_rows[[row_counter]] <- data.frame(
          dyad_id = dyad_counter, turn = half_turns + t,
          director_id = member_b, matcher_id = member_a
        )
      }
    }
  }
  
  data <- do.call(rbind, all_rows)
  data$director_diag <- participant_diag_lookup[data$director_id]
  data$matcher_diag <- participant_diag_lookup[data$matcher_id]
  data$status_match <- ifelse(data$director_diag == data$matcher_diag, 1, 0)
  
  n_participants <- n_blocks * 6
  n_dyads <- dyad_counter
  
  # Random effects
  director_re <- rnorm(n_participants, 0, random_effect_sd_director)
  dyad_re <- rnorm(n_dyads, 0, random_effect_sd_dyad)
  
  data$director_re <- director_re[data$director_id]
  data$dyad_re <- dyad_re[data$dyad_id]
  
  # Fixed effects design vectors
  X_director_AU <- as.numeric(data$director_diag == "AU")
  X_director_SZ <- as.numeric(data$director_diag == "SZ")
  X_matcher_AU  <- as.numeric(data$matcher_diag == "AU")
  X_matcher_SZ  <- as.numeric(data$matcher_diag == "SZ")
  
  turn_scaled <- (data$turn - 1) / (turns_per_dyad - 1)
  
  linear_pred <- log(baseline_rate) +
    effect_director_AU * X_director_AU + effect_director_SZ * X_director_SZ +
    effect_matcher_AU * X_matcher_AU + effect_matcher_SZ * X_matcher_SZ +
    effect_match_bonus * data$status_match +
    effect_turn * turn_scaled +
    data$director_re + data$dyad_re
  
  lambda <- exp(linear_pred)
  data$outcome <- rnbinom(nrow(data), mu = lambda, size = overdispersion)
  
  # Format factors for model fitting
  data$director_id <- factor(data$director_id)
  data$dyad_id <- factor(data$dyad_id)
  data$director_diag <- factor(data$director_diag, levels = c("NT", "AU", "SZ"))
  data$matcher_diag <- factor(data$matcher_diag, levels = c("NT", "AU", "SZ"))
  data$status_match <- factor(data$status_match, levels = c(0, 1),
                              labels = c("Mismatch", "Match"))
  
  return(data)
}

#### MODEL FITTING AND POWER CALCULATION ####

run_single_simulation <- function(sim_num, n_blocks, turns_per_dyad,
                                  effect_match_bonus,
                                  random_effect_sd_director,
                                  random_effect_sd_dyad) {
  
  tryCatch({
    data <- simulate_multiPartner_study(n_blocks, turns_per_dyad,
                                        effect_match_bonus,
                                        random_effect_sd_director,
                                        random_effect_sd_dyad)
    
    # PRIMARY MODEL: tests status_match (1 df) while controlling for each
    # partner's own diagnosis and turn position
    model <- suppressWarnings(
      glmer(outcome ~ status_match + director_diag + matcher_diag + turn +
              (1 | director_id) + (1 | dyad_id),
            data = data,
            family = negative.binomial(theta = overdispersion),
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 30000)))
    )
    
    if (length(model@optinfo$conv$lme4) > 0) {
      return(list(significant = NA, p_value = NA, converged = FALSE, chi_square = NA))
    }
    
    coef_summary <- summary(model)$coefficients
    match_row <- grep("status_matchMatch", rownames(coef_summary))
    
    if (length(match_row) == 0) {
      return(list(significant = NA, p_value = NA, converged = FALSE, chi_square = NA))
    }
    
    p_value <- coef_summary[match_row[1], "Pr(>|z|)"]
    z_value <- coef_summary[match_row[1], "z value"]
    chi_square <- z_value^2
    
    return(list(significant = p_value < 0.05, p_value = p_value,
                converged = TRUE, chi_square = chi_square))
    
  }, error = function(e) {
    if (sim_num == 1) cat("DEBUG - Error on sim 1:", conditionMessage(e), "\n")
    return(list(significant = NA, p_value = NA, converged = FALSE,
                chi_square = NA, error = as.character(e)))
  })
}

### ALTERNATIVE (with lower power): full 3x3 categorical interaction (director_diag x matcher_diag)

# run_single_simulation_full_interaction <- function(sim_num, n_blocks, turns_per_dyad) {
#   tryCatch({
#     data <- simulate_multiPartner_study(n_blocks, turns_per_dyad)
#     full_model <- suppressWarnings(glmer(
#       outcome ~ director_diag * matcher_diag + turn +
#         (1 | director_id) + (1 | dyad_id),
#       data = data, family = negative.binomial(theta = overdispersion),
#       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 30000))))
#     reduced_model <- suppressWarnings(glmer(
#       outcome ~ director_diag + matcher_diag + turn +
#         (1 | director_id) + (1 | dyad_id),
#       data = data, family = negative.binomial(theta = overdispersion),
#       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 30000))))
#     lrt <- anova(reduced_model, full_model)
#     p_value <- lrt$`Pr(>Chisq)`[2]
#     chi_square <- lrt$Chisq[2]
#     return(list(significant = p_value < 0.05, p_value = p_value,
#                 converged = TRUE, chi_square = chi_square))
#   }, error = function(e) {
#     return(list(significant = NA, p_value = NA, converged = FALSE, chi_square = NA))
#   })
# }

#### POWER ANALYSIS ACROSS SCENARIOS (total sample size) ####

n_blocks <- 10

cat("Testing", n_blocks, "blocks (", n_blocks * 6, "participants ) across",
    length(scenarios), "scenarios...\n\n")

set.seed(123)

scenario_single_results <- lapply(names(scenarios), function(scenario_name) {
  s <- scenarios[[scenario_name]]
  cat("--- Scenario:", scenario_name, "---\n")
  cat("  effect_match_bonus =", s$effect_match_bonus,
      "| sd_director =", s$sd_director, "| sd_dyad =", s$sd_dyad, "\n")
  
  results <- lapply(1:n_simulations, function(i) {
    run_single_simulation(i, n_blocks, turns_per_dyad,
                          s$effect_match_bonus, s$sd_director, s$sd_dyad)
  })
  
  significant <- sapply(results, function(x) x$significant)
  converged <- sapply(results, function(x) x$converged)
  chi_sq <- sapply(results, function(x) x$chi_square)
  
  power <- mean(significant[converged], na.rm = TRUE)
  conv_rate <- mean(converged, na.rm = TRUE)
  mean_chisq <- mean(chi_sq[converged], na.rm = TRUE)
  
  cat("  Power:", round(power * 100, 1), "% | Convergence:", round(conv_rate * 100, 1),
      "% | Mean chi-sq:", round(mean_chisq, 2), "\n\n")
  
  data.frame(scenario = scenario_name, n_blocks = n_blocks,
             n_participants = n_blocks * 6, power = power,
             convergence_rate = conv_rate, mean_chi_square = mean_chisq)
})

scenario_single_table <- do.call(rbind, scenario_single_results)
print(scenario_single_table)

#### POWER CURVE ACROSS SCENARIOS (varies n_blocks) ####

run_scenario_power_curve <- function(block_sizes = seq(2, 20, by = 2)) {
  
  all_results <- list()
  
  for (scenario_name in names(scenarios)) {
    s <- scenarios[[scenario_name]]
    cat("\n=== Scenario:", scenario_name, "===\n")
    
    power_results <- numeric(length(block_sizes))
    convergence_results <- numeric(length(block_sizes))
    
    for (i in seq_along(block_sizes)) {
      nb <- block_sizes[i]
      cat("Testing", nb, "blocks (", nb * 6, "participants )...\n")
      
      results <- lapply(1:n_simulations, function(j) {
        run_single_simulation(j, nb, turns_per_dyad,
                              s$effect_match_bonus, s$sd_director, s$sd_dyad)
      })
      
      significant <- sapply(results, function(x) x$significant)
      converged <- sapply(results, function(x) x$converged)
      
      power_results[i] <- mean(significant[converged], na.rm = TRUE)
      convergence_results[i] <- mean(converged, na.rm = TRUE)
    }
    
    all_results[[scenario_name]] <- data.frame(
      scenario = scenario_name,
      n_blocks = block_sizes,
      n_participants = block_sizes * 6,
      power = power_results,
      convergence_rate = convergence_results
    )
  }
  
  combined <- do.call(rbind, all_results)
  
  # Plot: overlay all three scenarios
  scenario_colors <- c(conservative = "firebrick", moderate = "goldenrod",
                       optimistic = "darkgreen")
  
  plot(NULL, xlim = range(combined$n_participants), ylim = c(0, 100),
       xlab = "Number of Participants", ylab = "Power (%)", las = 1,
       main = "Power curve testing effect of (mis)mtch across effect-size scenarios")
  abline(h = 80, col = "gray40", lty = 2, lwd = 2)
  grid()
  
  for (scenario_name in names(scenarios)) {
    sub <- combined[combined$scenario == scenario_name, ]
    lines(sub$n_participants, sub$power * 100, type = "b", pch = 19,
          col = scenario_colors[scenario_name], lwd = 2)
  }
  
  legend("bottomright",
         legend = paste0(names(scenarios), " (effect=",
                         sapply(scenarios, function(s) s$effect_match_bonus),
                         ", sd=", sapply(scenarios, function(s) s$sd_director), ")"),
         col = scenario_colors[names(scenarios)], lwd = 2, pch = 19,
         cex = 0.8, bg = "white")
  
  # Recommendations per scenario
  cat("\n--- SAMPLE SIZE RECOMMENDATIONS BY SCENARIO ---\n")
  for (scenario_name in names(scenarios)) {
    sub <- combined[combined$scenario == scenario_name, ]
    target_80 <- sub$n_participants[which(sub$power >= 0.80)[1]]
    cat(scenario_name, ": ", ifelse(is.na(target_80), "not reached in tested range",
                                    paste(target_80, "participants")), "\n")
  }
  
  return(combined)
}

scenario_curve_data <- run_scenario_power_curve()
print(scenario_curve_data)
```

