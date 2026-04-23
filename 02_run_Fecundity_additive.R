library(readr)
library(dplyr)
library(tidyr)
library(cmdstanr)
library(posterior)

setwd(dirname(normalizePath(sys.frame(1)$ofile)))

key_mat_levels <- c("AA", "AN", "NA", "NN")

d_all <- read_csv("fec_data.csv", show_col_types = FALSE)

d_all <- d_all %>%
  mutate(
    EE_seq = case_when(
      is.na(EE_seq) ~ "none",
      EE_seq == "NA" ~ "none",
      TRUE ~ EE_seq
    )
  )

d_fs <- d_all %>%
  filter(!is.na(fecundity)) %>%
  mutate(
    block = as.factor(block),
    pop = recode(
      pop,
      UNREL22 = "UNREL_22",
      UNREL120 = "UNREL_120"
    )
  ) %>%
  separate(pop, c("EE_type", "pop_id_raw"), sep = "_", remove = FALSE) %>%
  separate(block, c("B", "block_id_raw"), sep = 1, remove = FALSE) %>%
  mutate(
    block_id = as.integer(block_id_raw)
  ) %>%
  select(-B, -block_id_raw) %>%
  mutate(
    mat = case_when(
      grandmother == "A" & mother == "A" ~ "AA",
      grandmother == "A" & mother == "N" ~ "AN",
      grandmother == "N" & mother == "A" ~ "NA",
      grandmother == "N" & mother == "N" ~ "NN",
      TRUE ~ NA_character_
    ),
    mat = factor(mat, levels = key_mat_levels),
    mat_id = as.integer(mat),
    Ci = if_else(cue == "light", 1L, -1L),
    state_combo = interaction(mat, Ci, drop = TRUE),
    state_id = as.integer(factor(state_combo)),
    EE_code = case_when(
      EE == "ANC" ~ 1L,
      EE == "REL" ~ 2L,
      EE == "UNREL" ~ 3L,
      TRUE ~ 1L
    ),
    is_evo = if_else(EE_code == 1L, 0L, 1L),
    evo_type = case_when(
      EE_code == 1L ~ 0L,
      EE_code == 2L ~ 1L,
      EE_code == 3L ~ 2L
    ),
    pop_label = case_when(
      EE_code == 1L ~ "ANC",
      EE_code == 2L ~ paste0("REL_", pop_id_raw),
      EE_code == 3L ~ paste0("UNREL_", pop_id_raw)
    ),
    pop_id = as.integer(factor(pop_label)),
    block_id = as.integer(factor(block_id))
  ) %>%
  select(
    fecundity,
    EE_seq, EE, EE_code, is_evo, evo_type,
    grandmother, mother, cue,
    mat, mat_id, Ci, state_id,
    block, block_id,
    pop, pop_id, pop_label
  )

d_Evolved <- d_fs %>%
  droplevels()

d_Evolved_stan <- d_Evolved %>%
  mutate(
    evo_type = as.integer(evo_type),
    mat_id = as.integer(mat_id),
    state_id = as.integer(state_id),
    Ci = as.integer(Ci),
    block_id = as.integer(block_id)
  )

rel_levels <- d_Evolved_stan %>%
  filter(evo_type == 1L) %>%
  pull(pop_label) %>%
  unique()

unrel_levels <- d_Evolved_stan %>%
  filter(evo_type == 2L) %>%
  pull(pop_label) %>%
  unique()

d_Evolved_stan <- d_Evolved_stan %>%
  mutate(
    pop_rel_id = if_else(
      evo_type == 1L,
      as.integer(factor(pop_label, levels = rel_levels)),
      0L
    ),
    pop_unrel_id = if_else(
      evo_type == 2L,
      as.integer(factor(pop_label, levels = unrel_levels)),
      0L
    )
  )

J_rel <- length(rel_levels)
J_unrel <- length(unrel_levels)
B <- length(unique(d_Evolved_stan$block_id))
M <- length(unique(d_Evolved_stan$mat_id))
S <- length(unique(d_Evolved_stan$state_id))

stan_data <- list(
  N = nrow(d_Evolved_stan),
  fecundity = d_Evolved_stan$fecundity,
  evo_type = d_Evolved_stan$evo_type,
  M = M,
  mat_id = d_Evolved_stan$mat_id,
  S = S,
  state_id = d_Evolved_stan$state_id,
  Ci = d_Evolved_stan$Ci,
  B = B,
  block_id = d_Evolved_stan$block_id,
  J_rel = J_rel,
  pop_rel_id = d_Evolved_stan$pop_rel_id,
  J_unrel = J_unrel,
  pop_unrel_id = d_Evolved_stan$pop_unrel_id,
  use_state = c(0L, 0L, 0L)
)

stan_code <- '
functions {
  vector center_sum_to_zero(vector x_raw) {
    return x_raw - mean(x_raw);
  }
}

data {
  int<lower=1> N;
  vector<lower=0>[N] fecundity;

  array[N] int<lower=0, upper=2> evo_type;

  int<lower=1> M;
  array[N] int<lower=1, upper=M> mat_id;
  array[N] int<lower=-1, upper=1> Ci;

  int<lower=1> S;
  array[N] int<lower=1, upper=S> state_id;

  int<lower=1> B;
  array[N] int<lower=1, upper=B> block_id;

  int<lower=0> J_rel;
  int<lower=0> J_unrel;
  array[N] int<lower=0, upper=J_rel> pop_rel_id;
  array[N] int<lower=0, upper=J_unrel> pop_unrel_id;

  array[3] int<lower=0, upper=1> use_state;
}

parameters {
  vector[3] alpha;
  matrix[3, M] mat_raw;
  vector[3] beta_cue;
  matrix[3, S] state_raw;

  vector[B] z_block;
  real<lower=0> sigma_block;

  vector[J_rel] z_pop_rel;
  vector[J_unrel] z_pop_unrel;
  real<lower=0> sigma_pop;

  real<lower=1e-6> phi;
}

transformed parameters {
  vector[B] block_eff = center_sum_to_zero(z_block) * sigma_block;
  vector[J_rel] pop_rel_eff = rep_vector(0, J_rel);
  vector[J_unrel] pop_unrel_eff = rep_vector(0, J_unrel);
  matrix[3, M] mat_eff;
  matrix[3, S] state_eff;

  for (t in 1:3) {
    mat_eff[t] = to_row_vector(center_sum_to_zero(to_vector(mat_raw[t])));
    state_eff[t] = to_row_vector(center_sum_to_zero(to_vector(state_raw[t])));
  }

  if (J_rel > 0)
    pop_rel_eff = center_sum_to_zero(z_pop_rel) * sigma_pop;

  if (J_unrel > 0)
    pop_unrel_eff = center_sum_to_zero(z_pop_unrel) * sigma_pop;
}

model {
  alpha ~ normal(0, 1);
  to_vector(mat_raw) ~ normal(0, 0.3);
  to_vector(state_raw) ~ normal(0, 0.3);
  beta_cue ~ normal(0, 0.3);

  z_block ~ normal(0, 1);
  sigma_block ~ exponential(1);

  if (J_rel > 0)
    z_pop_rel ~ normal(0, 1);

  if (J_unrel > 0)
    z_pop_unrel ~ normal(0, 1);

  sigma_pop ~ exponential(1);
  phi ~ lognormal(1, 0.3);

  for (n in 1:N) {
    int t = evo_type[n] + 1;
    real eta = alpha[t];
    real rate;

    eta += block_eff[block_id[n]];

    if (evo_type[n] == 1) {
      if (pop_rel_id[n] > 0)
        eta += pop_rel_eff[pop_rel_id[n]];
    }

    if (evo_type[n] == 2) {
      if (pop_unrel_id[n] > 0)
        eta += pop_unrel_eff[pop_unrel_id[n]];
    }

    if (use_state[t] == 1) {
      eta += state_eff[t, state_id[n]];
    } else {
      eta += mat_eff[t, mat_id[n]] + beta_cue[t] * Ci[n];
    }

    rate = phi * exp(-eta);
    fecundity[n] ~ gamma(phi, rate);
  }
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    int t = evo_type[n] + 1;
    real eta = alpha[t];
    real rate;

    eta += block_eff[block_id[n]];

    if (evo_type[n] == 1) {
      if (pop_rel_id[n] > 0)
        eta += pop_rel_eff[pop_rel_id[n]];
    }

    if (evo_type[n] == 2) {
      if (pop_unrel_id[n] > 0)
        eta += pop_unrel_eff[pop_unrel_id[n]];
    }

    if (use_state[t] == 1) {
      eta += state_eff[t, state_id[n]];
    } else {
      eta += mat_eff[t, mat_id[n]] + beta_cue[t] * Ci[n];
    }

    rate = phi * exp(-eta);
    log_lik[n] = gamma_lpdf(fecundity[n] | phi, rate);
  }
}
'

stan_file <- write_stan_file(stan_code)
mod <- cmdstan_model(stan_file, force_recompile = TRUE)

fit <- mod$sample(
  data = stan_data,
  chains = 6,
  parallel_chains = 6,
  iter_warmup = 2000,
  iter_sampling = 4000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  seed = 42,
  refresh = 0
)

saveRDS(as_draws_df(fit$draws()), "stan_output/fit_fecundity_additive_draws_df.rds")

fit$cmdstan_diagnose()
