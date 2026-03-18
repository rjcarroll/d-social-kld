# D-SOCIAL-KLD

Replication and extension code for **D-SOCIAL-KLD scores** — a Bayesian
item-response-theory measure of corporate social responsibility (CSR) at the
firm-year level, covering 1991–2018.

The measure is introduced in:

> Carroll, R. J., Primo, D. M., & Richter, B. K. (2016). Using item response
> theory to improve measurement in strategic management research: An application
> to corporate social responsibility. *Strategic Management Journal*, 37(1),
> 66–85. https://doi.org/10.1002/smj.2463

Pre-computed scores are in [`data/output/`](data/output/). If you just want the
numbers, start there. If you want to understand or reproduce the estimation, read
on.

---

## What is D-SOCIAL-KLD?

Standard CSR indices (including the widely-used KLD Index) add up binary
indicators — strengths and concerns — as if each indicator were equally
informative. Item response theory (IRT) relaxes that assumption. It asks: given
all the binary choices a firm makes across dozens of CSR indicators, what is the
most likely underlying level of CSR responsibility? Indicators that better
discriminate between high- and low-CSR firms are weighted more heavily; easier
or noisier indicators are weighted less.

The **dynamic** extension allows a firm's latent CSR score to evolve from year
to year rather than treating each year as an independent cross-section. This
makes within-firm comparisons across time meaningful.

The result is a continuous, interval-scaled score for each firm in each year it
appears in the KLD data. Scores are centered roughly around zero; positive
scores indicate above-median CSR, negative scores below-median.

---

## Data outputs

| File | Description |
|------|-------------|
| `data/output/d-social-kld_scores.csv` | Firm-year D-SOCIAL-KLD scores: posterior mean, SD, and percentiles (5th through 95th) for the latent CSR parameter θ |
| `data/output/d-social-kld_item-params.csv` | Item-level parameters: posterior summaries for difficulty (α) and discrimination (β) for each KLD indicator in each year |

### Score file columns

| Column | Description |
|--------|-------------|
| `name` | Company name (from KLD) |
| `ticker` | Ticker symbol (see notes on disambiguation below) |
| `year` | Calendar year |
| `mean` | Posterior mean of θ — the primary D-SOCIAL-KLD score |
| `sd` | Posterior standard deviation — a direct measure of estimation uncertainty |
| `pct05`–`pct95` | Posterior percentiles |
| `min`, `max` | Extremes of the posterior sample |

The `sd` column is important. Unlike additive indices, the IRT approach
produces an explicit uncertainty estimate for every firm-year score. Researchers
should account for this uncertainty in downstream analyses (e.g., by using the
full posterior distribution rather than just the mean).

### Item parameters

The α (difficulty) and β (discrimination) columns characterize each KLD
indicator:

- **α (difficulty)**: Higher values mean fewer firms adopt this indicator at a
  given CSR level. High-α indicators are "harder" to earn (strengths) or harder
  to avoid (concerns).
- **β (discrimination)**: How steeply the probability of adoption changes with
  the firm's CSR level. Strength indicators typically have β > 0; concern
  indicators typically have β < 0. Near-zero β means the indicator does not
  effectively separate high- from low-CSR firms.

---

## Reproducing the scores

### Prerequisites

- R ≥ 4.2
- The following R packages (managed via `renv`):
  - `tidyverse`, `foreach`
  - `cmdstanr` (≥ 0.6) and a local [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) installation
  - `posterior`, `loo`
- The raw KLD STATS data (see below)

Install R packages:
```r
install.packages("renv")
renv::restore()
```

Install CmdStan (one-time):
```r
cmdstanr::install_cmdstan()
```

### Raw KLD data

The raw KLD STATS data (`KLD_stats_all_through_2018.csv`) is proprietary and
not distributed with this repository. It is available through
[Wharton Research Data Services (WRDS)](https://wrds-www.wharton.upenn.edu/)
under the MSCI ESG (formerly KLD) dataset.

Place the file at:
```
data/kld/KLD_stats_all_through_2018.csv
```

### Running the pipeline

Each script is self-contained and should be run in order from the project root:

```bash
Rscript R/01-clean.R      # ~5 min
Rscript R/02-estimate.R   # many hours — see note below
Rscript R/03-extract.R    # ~10 min
```

**Runtime note**: The estimation step (`02-estimate.R`) is computationally
intensive. On the full 1991–2018 dataset with 4 chains, expect runtimes of
several hours to a few days depending on hardware. A `demo_mode` flag in
`02-estimate.R` enables a fast subset run (1991–2012 only, fewer iterations)
for testing the pipeline end-to-end before committing to a full run. The default
demo subset is 1991–2000.

---

## The model

The core model is a dynamic two-parameter probit IRT model:

```
y[i,j,t] ~ Bernoulli(Φ(β[j,t] · θ[i,t] − α[j,t]))
```

where:
- **θ[i,t]** — latent CSR score for firm *i* in year *t* (the D-SOCIAL-KLD score)
- **α[j,t]** — difficulty of KLD indicator *j* in year *t*
- **β[j,t]** — discrimination of indicator *j* in year *t*
- **Φ(·)**   — the standard normal CDF (probit link)

The dynamic prior links a firm's score across years:

```
θ[i, first year] ~ Normal(0, σ_init)
θ[i,t] | θ[i,t−1] ~ Normal(θ[i,t−1], σ_firm[i])
```

where σ_firm[i] is a firm-specific innovation standard deviation estimated with
partial pooling across firms. This allows the model to learn how quickly each
firm's CSR position changes, while sharing information about the typical rate of
change across the full sample.

**Identification**: the latent scale is identified by constraining Halliburton
(HAL) to have a negative score in all years (low-CSR anchor) and Abbott
Laboratories (ABT) to have a positive score in all years (high-CSR anchor).
HAL is the same negative anchor used in the original paper. ABT replaces MSFT
as the positive anchor because MSFT and JNJ lack KLD indicator data before
2001; ABT has complete coverage for 1991–2018.

The full model is in [`stan/dynamic_irt.stan`](stan/dynamic_irt.stan), which
contains extensive inline documentation.

---

## Differences from the original paper

| Feature | Carroll et al. (2016) | This repo |
|---------|-----------------------|-----------|
| Coverage | 1991–2012 | 1991–2018 |
| Estimation engine | `MCMCpack::MCMCdynamicIRT1d` (Gibbs) | Stan (HMC/NUTS) |
| Diagnostics | Limited | R-hat, ESS, LOO-CV |
| Item parameters | Not reported | Included in output |
| Chains | 1 | 4 (default) |

Stan's Hamiltonian Monte Carlo sampler is generally more efficient than Gibbs
sampling — it produces less autocorrelated draws and provides richer convergence
diagnostics (R-hat, bulk and tail effective sample size, divergence checks).
The model specification is otherwise the same.

---

## Ticker disambiguation

Several firms share a ticker symbol across different KLD coverage universes. The
convention used here: the firm that first appeared in the dataset keeps its
ticker; later entrants receive a deterministic numeric suffix (e.g., `IBM.2`,
`IBM.3`). For years where `issuerid` is available (2017–2018), that field is
used to distinguish truly different firms sharing a ticker. The full
disambiguation logic is in [`R/01-clean.R`](R/01-clean.R) and
[`R/utils.R`](R/utils.R).

---

## Citation

If you use these scores, please cite the original paper:

```bibtex
@article{carroll2016irt,
  author  = {Carroll, Robert J. and Primo, David M. and Richter, Brian K.},
  title   = {Using item response theory to improve measurement in strategic
             management research: {An} application to corporate social
             responsibility},
  journal = {Strategic Management Journal},
  year    = {2016},
  volume  = {37},
  number  = {1},
  pages   = {66--85},
  doi     = {10.1002/smj.2463}
}
```

---

## Questions and contributions

Issues and pull requests are welcome. For questions about the underlying data or
methodology, the original paper and its online supplement are the primary
references. For questions about this codebase, open an issue.
