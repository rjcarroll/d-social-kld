// dynamic_irt.stan
//
// Dynamic two-parameter IRT model for D-SOCIAL-KLD scores.
//
// This implements the model in Carroll, Primo & Richter (2016, SMJ), following
// the dynamic ideal-point framework of Martin & Quinn (2002). The model treats
// KLD indicators as test items and firms as respondents, estimating an
// underlying latent dimension of corporate social responsibility (CSR) for each
// firm in each year.
//
// ── The probability model ───────────────────────────────────────────────────
//
//   For each observed (firm i, item j, year t) triple:
//
//     y[i,j,t] ~ Bernoulli(Phi(beta[j,t] * theta[i,t] - alpha[j,t]))
//
//   where Phi() is the standard-normal CDF (i.e., this is a probit link).
//
//   theta[i,t]  — latent CSR score for firm i in year t  ("responsibility")
//   alpha[j,t]  — difficulty of indicator j in year t: higher values mean
//                 fewer firms adopt the policy, all else equal
//   beta[j,t]   — discrimination of indicator j in year t: how well the
//                 indicator separates high- from low-CSR firms
//
//   Strength indicators tend to have beta > 0 (more-responsible firms are more
//   likely to have the strength coded 1). Concern indicators tend to have
//   beta < 0 (more-responsible firms are less likely to have the concern coded
//   1). The model estimates this from the data rather than imposing it.
//
// ── Dynamic prior on theta ──────────────────────────────────────────────────
//
//   First year a firm appears in the data:
//     theta[i, first_t] ~ Normal(0, sigma_init)
//
//   Each subsequent year:
//     theta[i,t] | theta[i,t-1] ~ Normal(theta[i,t-1], sigma_firm[i])
//
//   sigma_firm[i] is a firm-specific innovation standard deviation, estimated
//   with partial pooling via a hierarchical prior:
//     sigma_firm[i] ~ HalfNormal(0, sigma_global)
//     sigma_global   ~ HalfNormal(0, 1)
//     sigma_init     ~ HalfNormal(0, 1)
//
//   This structure lets the model learn how quickly each firm's CSR position
//   moves over time, while sharing information across firms about the typical
//   rate of change.
//
// ── Item priors ─────────────────────────────────────────────────────────────
//
//   alpha[j,t] ~ Normal(0, 5)
//   beta[j,t]  ~ Normal(0, 2.5)
//
//   These are mildly informative, regularizing priors. They are wider than
//   default Stan priors to allow the item parameters room to vary.
//
// ── Identification ──────────────────────────────────────────────────────────
//
//   IRT models are identified only up to a reflection and a scale. We fix
//   both by imposing sign constraints on two anchor firms:
//
//     Halliburton       (HAL): theta[HAL, t] < 0  for all t  →  low-CSR anchor
//     Abbott Labs       (ABT): theta[ABT, t] > 0  for all t  →  high-CSR anchor
//
//   HAL matches the negative anchor from Carroll, Primo & Richter (2016).
//   ABT replaces MSFT/JNJ as the positive anchor because MSFT and JNJ have no
//   KLD indicator data before 2001; ABT has full coverage for 1991–2018.
//   The sign constraints are encoded directly in the parameter declarations
//   below; Stan applies the appropriate Jacobian adjustments automatically.
//
//   Note: the scale of the theta distribution is determined by the item priors
//   and the data. Scores are interval-scaled (differences are meaningful;
//   ratios are not).
//
// ── Computational notes ──────────────────────────────────────────────────────
//
//   This model has O(N_firm_years) + O(N_items) parameters and an O(N)
//   likelihood, where N can be in the millions for the full KLD dataset. A
//   few things to keep in mind:
//
//   1. Non-centered parameterization: the dynamic prior is written in centered
//      form here for clarity. For very deep chains or slow mixing, switching
//      the transition terms to a non-centered form (z ~ Normal(0,1),
//      theta[curr] = theta[prev] + sigma_firm * z) can improve geometry.
//
//   2. Within-chain parallelism: the likelihood loop is the bottleneck. The
//      reduce_sum() function (Stan 2.23+) can parallelize it across CPU cores
//      within a single chain, giving a near-linear speedup. See the commented
//      skeleton at the bottom of this file.
//
//   3. For exploratory runs or testing, use a temporal or firm-count subset
//      of the data (see R/02-estimate.R for the `demo_mode` flag).

// ============================================================================
data {

  // ── dimensions ────────────────────────────────────────────────────────────

  int<lower=1> N;           // total non-missing observations
  int<lower=1> N_fy;        // total unique firm-year pairs
  int<lower=1> N_fy_free;   // firm-year pairs NOT belonging to either anchor
  int<lower=1> N_fy_hal;    // HAL (Halliburton) firm-year pairs  — low-CSR anchor
  int<lower=1> N_fy_abt;    // ABT (Abbott Labs) firm-year pairs  — high-CSR anchor
  int<lower=1> N_items;     // unique items (indicator × year combinations)
  int<lower=1> N_firms;     // unique firms

  // ── observations ──────────────────────────────────────────────────────────

  array[N] int<lower=0, upper=1> y;            // binary response (0/1)
  array[N] int<lower=1, upper=N_fy>    fy_obs; // firm-year index for obs n
  array[N] int<lower=1, upper=N_items> it_obs; // item index for obs n

  // ── firm-year → firm mapping ───────────────────────────────────────────────
  // Needed so the dynamic prior can look up sigma_firm[firm] for a given FY.

  array[N_fy] int<lower=1, upper=N_firms> fy_firm;

  // ── index maps: anchor/free slots → full FY index ─────────────────────────
  // These let us assemble theta from its three sub-vectors in transformed
  // parameters without any branching in the likelihood.

  array[N_fy_free] int<lower=1, upper=N_fy> free_to_fy;
  array[N_fy_hal]  int<lower=1, upper=N_fy> hal_to_fy;
  array[N_fy_abt]  int<lower=1, upper=N_fy> abt_to_fy;

  // ── dynamic model structure ───────────────────────────────────────────────

  // Initial firm-years: first time a firm appears in the data.
  int<lower=1>                                   N_init;
  array[N_init] int<lower=1, upper=N_fy>         init_fy;

  // Consecutive transitions: (prev, curr) pairs within the same firm.
  int<lower=0>                                   N_trans;
  array[N_trans] int<lower=1, upper=N_fy>        trans_prev;
  array[N_trans] int<lower=1, upper=N_fy>        trans_curr;
  array[N_trans] int<lower=1, upper=N_firms>     trans_firm;

}

// ============================================================================
parameters {

  // ── latent CSR scores ─────────────────────────────────────────────────────
  // Split into three sub-vectors so we can declare sign constraints on the
  // anchor firms. Stan applies Jacobian adjustments for the truncation
  // automatically.

  vector[N_fy_free]           theta_free; // unconstrained firm-years
  vector<upper=0>[N_fy_hal]   theta_hal;  // HAL (Halliburton): constrained negative
  vector<lower=0>[N_fy_abt]   theta_abt;  // ABT (Abbott Labs): constrained positive

  // ── item parameters ───────────────────────────────────────────────────────

  vector[N_items] alpha; // difficulty
  vector[N_items] beta;  // discrimination

  // ── dynamic variance parameters ───────────────────────────────────────────

  vector<lower=1e-8>[N_firms] sigma_firm;   // per-firm innovation SD
  real<lower=1e-8>            sigma_global; // hierarchical scale for sigma_firm
  real<lower=1e-8>            sigma_init;   // SD for firms' first-year draw

}

// ============================================================================
transformed parameters {

  // Assemble the full theta vector in firm-year order. This is the quantity
  // that enters the likelihood and the dynamic prior.
  vector[N_fy] theta;
  theta[free_to_fy] = theta_free;
  theta[hal_to_fy]  = theta_hal;
  theta[abt_to_fy]  = theta_abt;

}

// ============================================================================
model {

  // ── hyperpriors ───────────────────────────────────────────────────────────

  sigma_global ~ normal(0, 1);   // half-normal (positive constraint above)
  sigma_init   ~ normal(0, 1);   // half-normal

  // ── item priors ───────────────────────────────────────────────────────────

  alpha ~ normal(0, 5);
  beta  ~ normal(0, 2.5);

  // ── hierarchical prior on firm innovation SDs ─────────────────────────────

  sigma_firm ~ normal(0, sigma_global); // half-normal (positive constraint)

  // ── dynamic prior on theta ────────────────────────────────────────────────

  // First appearance: draw from a zero-centered normal.
  // (For anchor firms, the constraint on theta makes this a truncated normal;
  //  Stan accounts for the Jacobian automatically.)
  theta[init_fy] ~ normal(0, sigma_init);

  // Year-to-year transitions: random walk with firm-specific innovation SD.
  for (k in 1:N_trans) {
    theta[trans_curr[k]] ~ normal(theta[trans_prev[k]],
                                  sigma_firm[trans_firm[k]]);
  }

  // ── likelihood ────────────────────────────────────────────────────────────

  {
    // Compute the linear predictor for every observation, then evaluate the
    // Bernoulli-probit log-likelihood. Pulling this into a local block avoids
    // allocating `mu` in the global scope.
    // Vectorized Bernoulli-probit likelihood. Phi() can underflow to exactly
    // 0 for |mu| > ~37 during extreme HMC proposals; this simply causes the
    // proposal to be rejected (harmless). The initialization issue is handled
    // by starting all chains at init = 0 on the unconstrained scale
    // (see 02-estimate.R) so that mu ≈ 0 at startup and Phi(0) = 0.5.
    vector[N] mu;
    for (n in 1:N) {
      mu[n] = beta[it_obs[n]] * theta[fy_obs[n]] - alpha[it_obs[n]];
    }
    y ~ bernoulli(Phi(mu));
  }

}

// ============================================================================
generated quantities {

  // ── log-likelihood (for LOO-CV) ───────────────────────────────────────────
  // Storing pointwise log-likelihood enables posterior predictive checks and
  // model comparison via the `loo` package without re-running the model.
  //
  // Memory note: this is an N-length vector per posterior draw, which can be
  // large. For the full KLD dataset consider subsetting or dropping this block
  // and recomputing offline if memory is a constraint.

  vector[N] log_lik;
  for (n in 1:N) {
    real mu_n = beta[it_obs[n]] * theta[fy_obs[n]] - alpha[it_obs[n]];
    log_lik[n] = bernoulli_lpmf(y[n] | Phi(mu_n));
  }

}

// ============================================================================
// Parallelization sketch (reduce_sum, Stan 2.23+)
// ────────────────────────────────────────────────
// Uncomment and adapt this section if you want within-chain parallelism.
// You will also need to add `using stan::math::reduce_sum;` and compile with
// the STAN_THREADS flag. See: https://mc-stan.org/docs/stan-users-guide/
//
// functions {
//   real partial_log_lik(array[] int y_slice,
//                        int start, int end,
//                        vector theta, vector alpha, vector beta,
//                        array[] int fy_obs_slice,
//                        array[] int it_obs_slice) {
//     int n_slice = end - start + 1;
//     vector[n_slice] mu;
//     for (n in 1:n_slice)
//       mu[n] = beta[it_obs_slice[n]] * theta[fy_obs_slice[start + n - 1]]
//               - alpha[it_obs_slice[n]];
//     real ll = 0;
//     for (n in 1:n_slice)
//       ll += y_slice[n] == 1 ? std_normal_lcdf(mu[n] | ) : std_normal_lccdf(mu[n] | );
//     return ll;
//   }
// }
//
// In model block, replace the likelihood loop with:
//   target += reduce_sum(partial_log_lik, y,
//                        grainsize,
//                        theta, alpha, beta,
//                        fy_obs, it_obs);
