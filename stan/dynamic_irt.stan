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
//   both by pinning two anchor firms' first-year scores to fixed values:
//
//     Halliburton (HAL): theta[HAL, first_t] = -1  →  low-CSR anchor
//     Abbott Labs (ABT): theta[ABT, first_t] = +1  →  high-CSR anchor
//
//   HAL matches the negative anchor from Carroll, Primo & Richter (2016).
//   ABT replaces MSFT/JNJ as the positive anchor because MSFT and JNJ have no
//   KLD indicator data before 2001; ABT has full coverage for 1991–2018.
//
//   The fixed values are placed only on the first-year anchor scores;
//   subsequent years evolve freely via the NCP transition mechanism (see
//   below). Fixing at ±1 rather than using sign constraints eliminates the
//   scale degeneracy that arises when anchors can hover near zero.
//
//   Scores are interval-scaled (differences are meaningful; ratios are not).
//
// ── Parameterization ────────────────────────────────────────────────────────
//
//   This model uses a *non-centered parameterization* (NCP) for the dynamic
//   prior on theta and for the firm-level innovation SDs. The centered form
//   of the dynamic prior creates a funnel geometry (sigma_global → sigma_firm
//   → theta transitions) that causes HMC to hit max treedepth and produces
//   near-zero E-BFMI. NCP decouples these dependencies:
//
//     sigma_firm[i]  = sigma_global * sigma_firm_raw[i]
//                      sigma_firm_raw[i] ~ HalfNormal(0, 1)
//
//     theta[first_t] = sigma_init * z_init           (non-anchor firms)
//                      z_init ~ Normal(0, 1)
//
//     theta[t]       = theta[t-1] + sigma_firm[i] * z_trans[k]
//                      z_trans[k] ~ Normal(0, 1)
//
//   theta is then a deterministic function in transformed parameters, built
//   sequentially from the z innovations. The anchor first-year scores are
//   fixed constants (not sampled parameters).
//
// ── Computational notes ──────────────────────────────────────────────────────
//
//   1. Within-chain parallelism: the likelihood loop is the bottleneck. The
//      reduce_sum() function (Stan 2.23+) can parallelize it across CPU cores
//      within a single chain, giving a near-linear speedup. See the commented
//      skeleton at the bottom of this file.
//
//   2. log_lik (LOO-CV): pointwise log-likelihoods are not stored here because
//      for the full KLD dataset the N-length vector per draw causes memory
//      exhaustion. To compute LOO, re-run with a generated quantities block
//      added, or use moment matching on a subsample.
//
//   3. For exploratory runs or testing, use a temporal or firm-count subset
//      of the data (see R/02-estimate.R for the `demo_mode` flag).

// ============================================================================
data {

  // ── dimensions ────────────────────────────────────────────────────────────

  int<lower=1> N;        // total non-missing observations
  int<lower=1> N_fy;     // total unique firm-year pairs
  int<lower=1> N_items;  // unique items (indicator × year combinations)
  int<lower=1> N_firms;  // unique firms

  // ── observations ──────────────────────────────────────────────────────────

  array[N] int<lower=0, upper=1>    y;      // binary response (0/1)
  array[N] int<lower=1, upper=N_fy>    fy_obs; // firm-year index for obs n
  array[N] int<lower=1, upper=N_items> it_obs; // item index for obs n

  // ── firm-year → firm mapping ───────────────────────────────────────────────
  // Needed so the dynamic prior can look up sigma_firm[firm] for a given FY.

  array[N_fy] int<lower=1, upper=N_firms> fy_firm;

  // ── anchor first-year firm-years (for identification) ─────────────────────
  // These are the fy_idx values of HAL's and ABT's first appearances.
  // Only the first year is constrained; subsequent years follow from NCP
  // transitions.

  int<lower=1, upper=N_fy> hal_init_fy;  // HAL's first-year fy_idx
  int<lower=1, upper=N_fy> abt_init_fy;  // ABT's first-year fy_idx

  // ── non-anchor initial firm-years ─────────────────────────────────────────
  // fy_idx values of all other firms' first appearances.

  int<lower=0> N_init_other;
  array[N_init_other] int<lower=1, upper=N_fy> other_init_fy;

  // ── dynamic model structure ───────────────────────────────────────────────
  // Consecutive transitions: (prev, curr) pairs within the same firm,
  // ordered chronologically within each firm so that theta[trans_prev[k]] is
  // always set before theta[trans_curr[k]] in the transformed parameters loop.

  int<lower=0>                               N_trans;
  array[N_trans] int<lower=1, upper=N_fy>    trans_prev;
  array[N_trans] int<lower=1, upper=N_fy>    trans_curr;
  array[N_trans] int<lower=1, upper=N_firms> trans_firm;

}

// ============================================================================
parameters {

  // ── NCP innovations ───────────────────────────────────────────────────────
  // These are standard-normal raw parameters. The actual theta values are
  // built deterministically in transformed parameters.

  vector[N_init_other] z_init;   // first-year innovations for non-anchor firms
  vector[N_trans]      z_trans;  // year-to-year transition innovations

  // ── item parameters ───────────────────────────────────────────────────────

  vector[N_items] alpha;  // difficulty
  vector[N_items] beta;   // discrimination

  // ── NCP variance parameters ───────────────────────────────────────────────

  vector<lower=0>[N_firms] sigma_firm_raw;  // NCP: sigma_firm = sigma_global * this
  real<lower=0>            sigma_global;    // hierarchical scale for sigma_firm
  real<lower=0>            sigma_init;      // SD for firms' first-year draw

}

// ============================================================================
transformed parameters {

  // Actual firm innovation SDs, recovered from the NCP decomposition.
  vector<lower=0>[N_firms] sigma_firm = sigma_global * sigma_firm_raw;

  // ── Build theta[N_fy] from NCP innovations ────────────────────────────────
  //
  // Order of construction:
  //   1. Anchor first-year scores (constrained scalars set directly).
  //   2. All other firms' first-year scores (sigma_init * z_init).
  //   3. Consecutive transitions for all firms (theta[prev] + sigma_firm * z).
  //
  // The transition array is sorted (ticker, year_idx), so for each firm the
  // earlier year is always populated before the later year — the loop is
  // well-defined without any graph-topology checks.

  vector[N_fy] theta;

  // Anchor first-year scores fixed for identification (reflection + scale).
  theta[hal_init_fy] = -1.0;
  theta[abt_init_fy] =  1.0;

  if (N_init_other > 0)
    theta[other_init_fy] = sigma_init * z_init;

  for (k in 1:N_trans)
    theta[trans_curr[k]] = theta[trans_prev[k]]
                           + sigma_firm[trans_firm[k]] * z_trans[k];

}

// ============================================================================
model {

  // ── hyperpriors ───────────────────────────────────────────────────────────

  sigma_global ~ normal(0, 1);   // half-normal (positive constraint above)
  sigma_init   ~ normal(0, 1);   // half-normal

  // ── NCP priors ────────────────────────────────────────────────────────────

  sigma_firm_raw ~ normal(0, 1); // half-normal (positive constraint above)
  z_init         ~ normal(0, 1);
  z_trans        ~ normal(0, 1);

  // ── item priors ───────────────────────────────────────────────────────────

  alpha ~ normal(0, 5);
  beta  ~ normal(0, 2.5);

  // ── likelihood ────────────────────────────────────────────────────────────

  {
    vector[N] mu;
    for (n in 1:N)
      mu[n] = beta[it_obs[n]] * theta[fy_obs[n]] - alpha[it_obs[n]];
    y ~ bernoulli(Phi(mu));
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
