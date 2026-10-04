---
  title: "Power Simulation for MISMATCH Project - Sample size estimation for WP2"
author: "Philippine Geelhand"
date: "`r Sys.Date()`"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
```{r car}
## Within-participant design: AU1-AU2, NT1-NT2, AU1-NT1, AU2-NT2
## Four dyads per block (2 same-diagnosis, 2 mixed)
## Within each dyad, both members alternate as Director and Matcher.
## Child sample aged 8-13, recruited in two balanced groups (8-10 and 11-13)

# Sample size scales by replicating this block (n_blocks x 4 participants).
# Random effects: (1|director_id) + (1|dyad_id)

# Required libraries
library(lme4)
library(MASS)

#### CONFIGURATION SECTION #####

test_target <- "main"      # "main" tests only interaction of partner neurotype; "age_interaction" tests interaction of age with neurotype
n_blocks <- 35             # Test sample size (1 block = 4 participants, 4 dyads)
turns_per_dyad <- 10       # 5 as director + 5 as matcher per dyad
n_simulations <- 400       

# Age (continuous control)
effect_director_age <- 0.07
effect_matcher_age  <- 0.03
age_center <- 11           # Age is centered here (midpoint of the 8-14 span)

effect_director_AU <- 0.2
effect_matcher_AU  <- 0.2
effect_turn <- 0.1         # Linear trend across turns within a dyad

# Count data parameters, Geelhand et al.(2025) used as baseline in this WP2 as well
baseline_rate <- 5.25
overdispersion <- 3.68     # Fixed theta

#### Scenarios testing different effect sizes

scenarios <- list(
  conservative = list(effect_match_bonus = 0.2, sd_director = 0.35, sd_dyad = 0.35,
                      age_x_match = 0.03),
  moderate     = list(effect_match_bonus = 0.5, sd_director = 0.25, sd_dyad = 0.25,
                      age_x_match = 0.06),
  optimistic   = list(effect_match_bonus = 0.8, sd_director = 0.17, sd_dyad = 0.17,
                      age_x_match = 0.10)
)

build_block_dyads <- function(offset) {
  # Local indices: 1=AU1, 2=AU2, 3=NT1, 4=NT2
  local_diag <- c("AU", "AU", "NT", "NT")
  dyad_pairs <- list(c(1, 2),   # AU1-AU2 (same)
                     c(3, 4),   # NT1-NT2 (same)
                     c(1, 3),   # AU1-NT1 (mixed)
                     c(2, 4))   # AU2-NT2 (mixed)
  
  list(
    participant_id = offset + 1:4,
    diag = local_diag,
    dyad_pairs = lapply(dyad_pairs, function(p) offset + p)
  )
}

####  SIMULATION FUNCTION ####

simulate_two_partner_study <- function(n_blocks, turns_per_dyad,
                                       effect_match_bonus,
                                       sd_director, sd_dyad,
                                       effect_age_x_match) {
  
  half_turns <- turns_per_dyad / 2
  n_participants <- n_blocks * 4
  diag_lookup <- character(n_participants)
  
  all_rows <- list()
  dyad_counter <- 0
  
  for (b in 1:n_blocks) {
    block <- build_block_dyads((b - 1) * 4)
    diag_lookup[block$participant_id] <- block$diag
    
    for (pair in block$dyad_pairs) {
      dyad_counter <- dyad_counter + 1
      member_a <- pair[1]
      member_b <- pair[2]
      
      # First half: A directs / B matches. Second half: roles swap.
      all_rows[[dyad_counter]] <- data.frame(
        dyad_id = dyad_counter,
        turn = 1:turns_per_dyad,
        director_id = c(rep(member_a, half_turns), rep(member_b, half_turns)),
        matcher_id  = c(rep(member_b, half_turns), rep(member_a, half_turns))
      )
    }
  }
  
  data <- do.call(rbind, all_rows)
  n_dyads <- dyad_counter
  
  # Diagnosis and match status
  data$director_diag <- diag_lookup[data$director_id]
  data$matcher_diag  <- diag_lookup[data$matcher_id]
  data$status_match  <- ifelse(data$director_diag == data$matcher_diag, 1, 0)
  
  # Age: balanced recruitment bands, exact age uniform within each band.
  # Within each diagnosis, half the children come from 8.0-10.99 and half from 11.0-13.99
  age <- numeric(n_participants)
  for (d in c("AU", "NT")) {
    ids <- which(diag_lookup == d)
    band <- sample(rep(c("Young", "Old"), length.out = length(ids)))
    age[ids] <- ifelse(band == "Young",
                       runif(length(ids), 8, 11),
                       runif(length(ids), 11, 14))
  }
  age_c <- age - age_center
  data$director_age_c <- age_c[data$director_id]
  data$matcher_age_c  <- age_c[data$matcher_id]
  
  # Turn position scaled 0-1
  data$turn_c <- (data$turn - 1) / (turns_per_dyad - 1)
  
  # Random effects
  director_re <- rnorm(n_participants, 0, sd_director)
  dyad_re <- rnorm(n_dyads, 0, sd_dyad)
  
  # Linear predictor
  linear_pred <- log(baseline_rate) +
    effect_director_AU * as.numeric(data$director_diag == "AU") +
    effect_matcher_AU  * as.numeric(data$matcher_diag == "AU") +
    effect_match_bonus * data$status_match +
    effect_director_age * data$director_age_c +
    effect_matcher_age  * data$matcher_age_c +
    effect_age_x_match * data$status_match * data$director_age_c +
    effect_turn * data$turn_c +
    director_re[data$director_id] +
    dyad_re[data$dyad_id]
  
  data$outcome <- rnbinom(nrow(data), mu = exp(linear_pred), size = overdispersion)
  
  # Format factors for model fitting
  data$director_id <- factor(data$director_id)
  data$dyad_id <- factor(data$dyad_id)
  data$director_diag <- factor(data$director_diag, levels = c("NT", "AU"))
  data$matcher_diag  <- factor(data$matcher_diag,  levels = c("NT", "AU"))
  data$status_match <- factor(data$status_match, levels = c(0, 1),
                              labels = c("Mismatch", "Match"))
  
  return(data)
}

#### MODEL FITTING AND POWER CALCULATION ####

run_single_simulation <- function(sim_num, n_blocks, turns_per_dyad,
                                  effect_match_bonus, sd_director, sd_dyad,
                                  age_x_match) {
  
  tryCatch({
    # The age interaction is only built into the data when it is the target
    true_age_x_match <- if (test_target == "age_interaction") age_x_match else 0
    
    data <- simulate_two_partner_study(n_blocks, turns_per_dyad,
                                       effect_match_bonus, sd_director, sd_dyad,
                                       true_age_x_match)
    
    if (test_target == "age_interaction") {
      fixed_part <- paste("status_match * director_age_c + director_diag +",
                          "matcher_diag + matcher_age_c + turn_c")
      target_term <- "status_matchMatch:director_age_c"
    } else {
      fixed_part <- paste("status_match + director_diag + matcher_diag +",
                          "director_age_c + matcher_age_c + turn_c")
      target_term <- "status_matchMatch"
    }
    model_formula <- as.formula(paste("outcome ~", fixed_part,
                                      "+ (1 | director_id) + (1 | dyad_id)"))
    
    model <- suppressWarnings(
      glmer(model_formula,
            data = data,
            family = negative.binomial(theta = overdispersion),
            control = glmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 30000)))
    )
    
    if (length(model@optinfo$conv$lme4) > 0) {
      return(list(significant = NA, p_value = NA, converged = FALSE, chi_square = NA))
    }
    
    coef_summary <- summary(model)$coefficients
    target_row <- which(rownames(coef_summary) == target_term)
    
    if (length(target_row) == 0) {
      return(list(significant = NA, p_value = NA, converged = FALSE, chi_square = NA))
    }
    
    p_value <- coef_summary[target_row[1], "Pr(>|z|)"]
    z_value <- coef_summary[target_row[1], "z value"]
    
    return(list(significant = p_value < 0.05, p_value = p_value,
                converged = TRUE, chi_square = z_value^2))
    
  }, error = function(e) {
    if (sim_num == 1) cat("DEBUG - Error on sim 1:", conditionMessage(e), "\n")
    return(list(significant = NA, p_value = NA, converged = FALSE,
                chi_square = NA, error = as.character(e)))
  })
}

cat("Test target:", test_target, "\n")
cat("Testing", n_blocks, "blocks (", n_blocks * 4, "participants,",
    n_blocks * 4, "dyads ) across", length(scenarios), "scenarios...\n\n")

set.seed(123)

scenario_single_results <- lapply(names(scenarios), function(scenario_name) {
  s <- scenarios[[scenario_name]]
  cat("--- Scenario:", scenario_name, "---\n")
  cat("  effect_match_bonus =", s$effect_match_bonus,
      "| age_x_match =", s$age_x_match,
      "| sd_director =", s$sd_director, "| sd_dyad =", s$sd_dyad, "\n")
  
  results <- lapply(1:n_simulations, function(i) {
    run_single_simulation(i, n_blocks, turns_per_dyad,
                          s$effect_match_bonus, s$sd_director, s$sd_dyad,
                          s$age_x_match)
  })
  
  significant <- sapply(results, function(x) x$significant)
  converged <- sapply(results, function(x) x$converged)
  chi_sq <- sapply(results, function(x) x$chi_square)
  
  power <- mean(significant[converged], na.rm = TRUE)
  conv_rate <- mean(converged, na.rm = TRUE)
  
  cat("  Power:", round(power * 100, 1), "% | Convergence:",
      round(conv_rate * 100, 1), "%\n\n")
  
  data.frame(scenario = scenario_name, n_blocks = n_blocks,
             n_participants = n_blocks * 4, power = power,
             convergence_rate = conv_rate,
             mean_chi_square = mean(chi_sq[converged], na.rm = TRUE))
})

print(do.call(rbind, scenario_single_results))

#### POWER CURVE ACROSS SCENARIOS ####

run_scenario_power_curve <- function(block_sizes = seq(5, 50, by = 5)) {
  
  all_results <- list()
  
  for (scenario_name in names(scenarios)) {
    s <- scenarios[[scenario_name]]
    cat("\n=== Scenario:", scenario_name, "===\n")
    
    power_results <- numeric(length(block_sizes))
    convergence_results <- numeric(length(block_sizes))
    
    for (i in seq_along(block_sizes)) {
      nb <- block_sizes[i]
      cat("Testing", nb, "blocks (", nb * 4, "participants )...\n")
      
      results <- lapply(1:n_simulations, function(j) {
        run_single_simulation(j, nb, turns_per_dyad,
                              s$effect_match_bonus, s$sd_director, s$sd_dyad,
                              s$age_x_match)
      })
      
      significant <- sapply(results, function(x) x$significant)
      converged <- sapply(results, function(x) x$converged)
      
      power_results[i] <- mean(significant[converged], na.rm = TRUE)
      convergence_results[i] <- mean(converged, na.rm = TRUE)
    }
    
    all_results[[scenario_name]] <- data.frame(
      scenario = scenario_name,
      n_blocks = block_sizes,
      n_participants = block_sizes * 4,
      power = power_results,
      convergence_rate = convergence_results
    )
  }
  
  combined <- do.call(rbind, all_results)
  
  # Plot: overlay all scenarios
  scenario_colors <- c(conservative = "firebrick", moderate = "goldenrod",
                       optimistic = "darkgreen")
  plot_title <- ifelse(test_target == "age_interaction",
                       "Power Curve Across Scenarios: Match x Age Interaction",
                       "Power Curve Across Scenarios: Match vs Mismatch Effect")
  
  plot(NULL, xlim = range(combined$n_participants), ylim = c(0, 100),
       xlab = "Number of Participants", ylab = "Power (%)", las = 1,
       main = plot_title)
  abline(h = 80, col = "gray40", lty = 2, lwd = 2)
  grid()
  
  for (scenario_name in names(scenarios)) {
    sub <- combined[combined$scenario == scenario_name, ]
    lines(sub$n_participants, sub$power * 100, type = "b", pch = 19,
          col = scenario_colors[scenario_name], lwd = 2)
  }
  
  effect_labels <- if (test_target == "age_interaction") {
    sapply(scenarios, function(s) s$age_x_match)
  } else {
    sapply(scenarios, function(s) s$effect_match_bonus)
  }
  legend("bottomright",
         legend = paste0(names(scenarios), " (effect=", effect_labels,
                         ", sd=", sapply(scenarios, function(s) s$sd_director), ")"),
         col = scenario_colors[names(scenarios)], lwd = 2, pch = 19,
         cex = 0.8, bg = "white")
  
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
