# Additional numerical evidence

- [The complete R market model](#the-complete-r-market-model)
- [1. Experimental scenarios and optimizer
  diagnostics](#1-experimental-scenarios-and-optimizer-diagnostics)
- [2. Distribution of high-risk households and cumulative
  floors](#2-distribution-of-high-risk-households-and-cumulative-floors)
- [3. Quote versus in-force accounting and presence
  constraints](#3-quote-versus-in-force-accounting-and-presence-constraints)
- [4. Guarded feasibility and corrective joint
  projection](#4-guarded-feasibility-and-corrective-joint-projection)
- [5. Renewal intensity, damping and
  convergence](#5-renewal-intensity-damping-and-convergence)
- [6. Actuarial benchmark deviations, not insurer-to-insurer
  transfers](#6-actuarial-benchmark-deviations-not-insurer-to-insurer-transfers)
- [7. Recalculated finite unilateral-deviation
  diagnostic](#7-recalculated-finite-unilateral-deviation-diagnostic)
- [8. Independent static best-response check](#8-Independent-static-best-response-check) 
- [10. Environment and replication
  audit](#10-environment-and-replication-audit)


# The complete R market model

The following executable R chunk is the *full* market model. 
It constructs municipalities and three risk
classes from a seeded Dirichlet distribution; maintains the in-force
premium recursively; solves each insurer’s locally constrained pricing
problem; and records optimization, feasibility, exposure, and
quote-versus-book diagnostics. 

# 1. Experimental scenarios and optimizer diagnostics

The main six screens cross two monitoring units (commune and risk class)
with no cumulative floor, a 50% retained-initial-share floor, and a 75%
retained-initial-share floor. All have the same underlying generated
municipalities, $\lambda=0.35$, $\nu=0.7$, periodwise relative cap
$\varepsilon=0.01$, and seeded five-start local optimization in `full`
mode. The `PRICE_UB=Inf` source still translates to a **finite numerical
box of $10^6$** inside the SQP solver; we explicitly test finite price
caps below. A local solver success code is not a global optimum
certificate.

<div id="tbl-runs">

Table 1

**Table S0. Actual model runs performed in this render; failed optional
runs are reported, not replaced by archived data.**

| Run                     | Status    |   t | Warnings | Seconds |
|:------------------------|:----------|----:|---------:|--------:|
| strict_commune          | COMPLETED | 200 |        0 |  4446.9 |
| strict_risk_class       | COMPLETED | 200 |        0 |  8125.7 |
| cum50_commune           | COMPLETED | 149 |        0 |  2687.7 |
| cum50_risk_class        | COMPLETED | 156 |        0 |  3933.3 |
| cum25_commune           | COMPLETED | 200 |        0 |  3595.2 |
| cum25_risk_class        | COMPLETED | 101 |        0 |  2559.1 |
| lambda_unconstrained_20 | COMPLETED |  66 |        0 |   891.6 |
| lambda_commune_20       | COMPLETED |  95 |        0 |  1330.9 |
| lambda_risk_class_20    | COMPLETED |  90 |        0 |  2170.8 |
| lambda_unconstrained_35 | COMPLETED |  35 |        0 |   514.5 |
| lambda_commune_35       | COMPLETED | 250 |        0 |  5375.5 |
| lambda_risk_class_35    | COMPLETED | 101 |        0 |  2558.8 |
| lambda_unconstrained_70 | COMPLETED |  28 |        0 |   359.8 |
| lambda_commune_70       | COMPLETED | 250 |        0 |  4721.3 |
| lambda_risk_class_70    | COMPLETED | 192 |        0 |  5013.2 |
| guard_risk_0.5          | COMPLETED | 120 |        0 |  2997.2 |
| guard_risk_0.7          | COMPLETED |  90 |        0 |  2275.7 |
| guard_risk_1.0          | COMPLETED | 120 |        0 |  3067.7 |
| guard_commune           | COMPLETED |  45 |        0 |   721.6 |
| cap_unconstrained_50    | COMPLETED |  35 |        0 |   369.1 |
| cap_unconstrained_200   | COMPLETED |  35 |        0 |   399.6 |

</div>

**Table S1. Price bounds and optimizer diagnostics recalculated from
fresh simulations.**

| Run               | Updates | Starts | Max_quote | Fallbacks | Near_upper |
|:------------------|--------:|-------:|----------:|----------:|-----------:|
| cum25_commune     |     600 |      5 |     5.122 |         0 |          0 |
| cum25_risk_class  |     303 |      5 |     4.933 |         0 |          0 |
| cum50_commune     |     447 |      5 |     5.950 |         0 |          0 |
| cum50_risk_class  |     468 |      5 |     5.806 |         0 |          0 |
| strict_commune    |     600 |      5 |     6.382 |         0 |          0 |
| strict_risk_class |     600 |      5 |     6.384 |         0 |          0 |

![Figure S1. Maximum selected quote by iteration in the six newly
computed synchronous experiments.](figures/S1_optimizer_price_paths.png)

<div id="tbl-price-ceiling-stress">

Table 2

**Table S2. Independent unconstrained reruns with explicit finite upper
price boxes; results are local numerical diagnostics.**

| Cap     | Status    | Max_quote | C_profit | C_high_N |   t |
|:--------|:----------|----------:|---------:|---------:|----:|
| cap_50  | COMPLETED |     6.384 |   7505.4 |    344.7 |  35 |
| cap_200 | COMPLETED |     6.384 |   7505.4 |    344.7 |  35 |

</div>

# 2. Distribution of high-risk households and cumulative floors

The floor is 75% of each
insurer’s initial share in the monitored unit; a uniform initial share
of one third therefore implies a 25% class floor. The following figures
show why a commune-level floor need not prohibit **within-commune
selection**. The experiment uses synthetic households and does not
estimate French market-wide behavior.

<div id="tbl-cumulative-results">

Table 3

**Table S3. High-risk outcomes of insurer C; premiums quoted versus
premiums on the in-force book are distinguished.**

| Regime            |   t |     N | Quote_LR | Book_LR | Book_margin |
|:------------------|----:|------:|---------:|--------:|------------:|
| strict_commune    | 200 | 344.2 |    0.761 |   0.761 |       399.7 |
| strict_risk_class | 200 | 344.7 |    0.761 |   0.761 |       400.4 |
| cum50_commune     | 149 | 214.3 |    0.768 |   0.768 |       232.2 |
| cum50_risk_class  | 156 | 474.8 |    0.859 |   0.859 |       290.4 |
| cum25_commune     | 200 | 187.2 |    0.894 |   0.894 |        85.4 |
| cum25_risk_class  | 101 | 712.2 |    1.192 |   1.192 |      -426.5 |

</div>

![Figure S2. High-risk insurer shares under six simulated constraint
regimes.](figures/S2_high_risk_allocation.png)

![Figure S3. In-force underwriting contribution of granular insurer C
from the high-risk class.](figures/S3_high_risk_underwriting.png)

<div id="tbl-selected-market-shares">

Table 4

**Table S4. High-risk policies and insurer shares in the three main
comparison regimes.**

| Run               | Insurer |      N | Share_pct | Book_profit |
|:------------------|:--------|-------:|----------:|------------:|
| strict_risk_class | A       | 1372.5 |      48.2 |        -740 |
| strict_risk_class | B       | 1131.3 |      39.7 |         -63 |
| strict_risk_class | C       |  344.7 |      12.1 |         400 |
| cum25_commune     | A       | 1716.5 |      60.3 |       -3179 |
| cum25_commune     | B       |  944.9 |      33.2 |       -1081 |
| cum25_commune     | C       |  187.2 |       6.6 |          85 |
| cum25_risk_class  | A       | 1424.3 |      50.0 |       -2299 |
| cum25_risk_class  | B       |  712.2 |      25.0 |        -433 |
| cum25_risk_class  | C       |  712.2 |      25.0 |        -427 |

</div>

# 3. Quote versus in-force accounting and presence constraints

Expected claims are evaluated
on insured risk exposures; portfolio accounting uses **in-force average
premiums** (`*_pbar`) rather than current quotes (`*_prix`). A segment
loss ratio exceeding one means expected claims exceed premium revenue in
that segment. It is not equivalent to total insurer insolvency.

<div id="tbl-book-classes">

Table 5

**Table S5. Book-based price, pure expected cost, loss ratio and margin
for C by risk class.**

| Regime            | Risk     |      N | Book_price |  Cost | Book_LR | Margin |
|:------------------|:---------|-------:|-----------:|------:|--------:|-------:|
| strict_risk_class | Low      | 2060.5 |      2.398 | 0.087 |   0.036 | 4762.9 |
| strict_risk_class | Moderate | 1394.1 |      3.304 | 1.624 |   0.492 | 2342.1 |
| strict_risk_class | High     |  344.7 |      4.862 | 3.700 |   0.761 |  400.4 |
| cum25_commune     | Low      | 1689.8 |      1.790 | 0.087 |   0.049 | 2877.3 |
| cum25_commune     | Moderate |  997.8 |      2.605 | 1.636 |   0.628 |  967.7 |
| cum25_commune     | High     |  187.2 |      4.314 | 3.858 |   0.894 |   85.4 |
| cum25_risk_class  | Low      | 1759.3 |      2.012 | 0.086 |   0.043 | 3388.5 |
| cum25_risk_class  | Moderate | 1004.7 |      3.011 | 1.598 |   0.531 | 1419.5 |
| cum25_risk_class  | High     |  712.2 |      3.126 | 3.725 |   1.192 | -426.5 |

</div>

<div id="tbl-synchronous-constraint-slacks">

Table 6

**Table S6. Diagnosed shares in the synchronous price-update screens,
which may include negative post-update slacks.**

| Run              | Rule       | Firm | Binding_pct | Min_slack |
|:-----------------|:-----------|:-----|------------:|:----------|
| cum25_commune    | cumulative | A    |         4.1 | -1.10e-04 |
| cum25_commune    | period     | A    |         0.1 | -1.49e-03 |
| cum25_commune    | cumulative | B    |        45.5 | -3.91e-04 |
| cum25_commune    | period     | B    |         0.0 | 2.32e-04  |
| cum25_commune    | cumulative | C    |        34.8 | -3.56e-04 |
| cum25_commune    | period     | C    |         0.2 | -1.08e-03 |
| cum25_risk_class | cumulative | A    |         8.6 | -1.11e-05 |
| cum25_risk_class | period     | A    |         1.3 | -6.84e-03 |
| cum25_risk_class | cumulative | B    |         8.9 | -3.06e-05 |
| cum25_risk_class | period     | B    |         0.7 | -3.94e-04 |
| cum25_risk_class | cumulative | C    |         8.9 | -2.91e-05 |
| cum25_risk_class | period     | C    |         0.7 | -4.90e-04 |

</div>

# 4. Guarded feasibility and corrective joint projection

For the guarded regime, we record
the original jointly proposed candidate and the actually accepted joint
vector, including backtracking and coordinated repairs. **An accepted
feasible update does not establish individual optimality.** A future
iteration after a guarded path’s recorded horizon is not certified
merely because the previous iterations were feasible.

<div id="tbl-guard-diagnostics">

Table 7

**Table S7. Feasibility accounting for independent guarded runs, with
their actual horizons.**

| Nu | Status | t | Proposed_infeasible | Accepted_infeasible | Backtracks | Joint_repairs | Min_slack |
|---:|:---|---:|---:|---:|---:|---:|:---|
| 0.5 | COMPLETED | 120 | 12 | 0 | 12 | 0 | 1.23e-06 |
| 0.7 | COMPLETED | 90 | 20 | 0 | 16 | 4 | -9.25e-09 |
| 1.0 | COMPLETED | 120 | 68 | 0 | 44 | 24 | -1.00e-08 |

**Table S7b. Additional guarded commune pilot; feasibility only, not
long-run convergence.**

| Protocol        |   t | Min_accepted | Projections | Maximum_price_residual |
|:----------------|----:|-------------:|------------:|-----------------------:|
| Guarded commune |  45 |            0 |           3 |                      0 |

</div>

![Figure S4. Proposed versus accepted vectors for
nu=0.7.](figures/S4_proposed_vs_accepted_slack.png)

# 5. Renewal intensity, damping and convergence

Renewal fraction $\lambda$ changes
the economic reallocation of customers and the mix of premiums collected
on the retained book. Damping $\nu$ changes the numerical step toward
quoted best responses. These are **not the same parameter**. Experiments
ending at the iteration cap are explicitly identified as not converged
to the prescribed numerical threshold.

<div id="tbl-lambda">

Table 8

**Table S8. Independent renewal-sensitivity simulations; stopping is a
local numerical criterion only.**

| Regime        | Lambda |   t | Strict_stop | N_high | LR_high | Profit |
|:--------------|-------:|----:|:------------|-------:|--------:|-------:|
| unconstrained |   0.20 |  66 | TRUE        |  344.7 |   0.761 |   7505 |
| commune       |   0.20 |  95 | TRUE        |  187.2 |   0.894 |   3930 |
| risk_class    |   0.20 |  90 | TRUE        |  712.2 |   1.192 |   4381 |
| unconstrained |   0.35 |  35 | TRUE        |  344.7 |   0.761 |   7505 |
| commune       |   0.35 | 250 | FALSE       |  187.2 |   0.894 |   3930 |
| risk_class    |   0.35 | 101 | TRUE        |  712.2 |   1.192 |   4381 |
| unconstrained |   0.70 |  28 | TRUE        |  344.7 |   0.761 |   7505 |
| commune       |   0.70 | 250 | FALSE       |  187.2 |   0.894 |   3930 |
| risk_class    |   0.70 | 192 | TRUE        |  712.2 |   1.192 |   4381 |

</div>

<div id="tbl-nu">

Table 9

**Table S9. Guarded damping sensitivity; unequal run horizons are
reported.**

|  Nu |   t |  N_high | LR_high | Price_residual |
|----:|----:|--------:|--------:|:---------------|
| 0.5 | 120 | 712.160 |   1.192 | 1.74e-05       |
| 0.7 |  90 | 712.156 |   1.191 | 4.30e-07       |
| 1.0 | 120 | 712.156 |   1.191 | 2.39e-07       |

</div>

![Figure S5. Stopping iteration by renewal fraction; open symbols
represent trajectories truncated at the simulation
cap.](figures/S5_lambda_convergence.png)

![Figure S6. Full price steps, scalar backtracking, and coordinated
price projections by numerical
damping.](figures/S6_nu_correction_frequencies.png)

# 6. Actuarial benchmark deviations, not insurer-to-insurer transfers

A premium above the loaded actuarial
benchmark is not a legislated tax paid to another insurer. A premium
below that benchmark is not a cash subsidy. The originally named
`tax_paid`/`subsidy_received` components are better described as
**above-reference**/**below-reference** amounts; their aggregate signed
difference need not be zero.

<div id="tbl-reference-deviations">

Table 10

**Table S10. Premium differences relative to the loaded actuarial
benchmark, NOT balanced transfers.**

| Regime            |   Above |  Below |     Net |
|:------------------|--------:|-------:|--------:|
| strict_commune    | 16453.0 | 1190.9 | 15262.1 |
| strict_risk_class | 16461.2 | 1200.3 | 15260.9 |
| cum50_commune     | 13396.3 | 2954.9 | 10441.4 |
| cum50_risk_class  | 15590.6 | 1529.8 | 14060.8 |
| cum25_commune     | 10168.2 | 4544.3 |  5623.9 |
| cum25_risk_class  | 12236.4 | 3402.9 |  8833.6 |

</div>

![Figure S7. Quoted premiums above and below the loaded actuarial
pricing benchmark.](figures/S7_actuarial_reference_deviations.png)

# 7. Recalculated finite unilateral-deviation diagnostic

This is a **new, finite
perturbation screen**, calculated at the final guarded risk-class state
rather than loaded from the earlier archived deviation CSV. Each insurer
is allowed to shift either a single pricing component or all of its
pricing components by -10%, -5%, +5%, or +10%. The other insurers’
quotes are fixed, while next-period shares respond through the same
logit and inertia equations. Constraints are evaluated jointly over all
firms and risk classes. Even if no feasible profitable perturbation is
found, this is not a proof of GNE or global unilateral optimality.

<div id="tbl-deviation-grid">

Table 11

**Table S11. New finite one-coordinate and simultaneous-within-insurer
deviation scan; relaxed tolerance 5e-8.**

| Firm |   N | Strict_feasible | Relaxed_feasible | Max_gain | Baseline_next_slack |
|:-----|----:|----------------:|-----------------:|---------:|:--------------------|
| A    |   8 |               0 |                0 |       NA | -1.21e-08           |
| B    |  84 |               0 |                0 |       NA | -1.21e-08           |
| C    | 244 |              40 |              120 |   -0.001 | -1.21e-08           |

</div>

# 8. Independent static best-response check 

Run:

```tet
source("check_best_responses.R")
```

It initializes 20 synthetic communes and 60 risk cells using seed 123; it solves 81 one-dimensional best-response problems (one national, 20 commune, 60 cell prices), independently of the SQP solver. 

```text
R1 pooled quote 3.2298046700; 20 communes; 60 risk cells
R1 iter 1: price residual 2.18; C high-risk 0.13427
R1 iter 5: price residual 0.0834; C high-risk 0.12103
R1 iter 10: price residual 0.00147; C high-risk 0.12104
R1 iter 15: price residual 9.51e-05; C high-risk 0.12102
R1 iter 20: price residual 5.33e-06; C high-risk 0.12102
R1 iter 25: price residual 3.07e-07; C high-risk 0.12102
R1 iter 30: price residual 1.77e-08; C high-risk 0.12102
R1 iter 35: price residual 1.03e-09; C high-risk 0.12102
R1 iter 36: price residual 5.8e-10; C high-risk 0.12102
No SQP RDS provided; independent static profile was validated.
```

The validation files report (approximate display precision):

| Insurer | Profit at profile | Numerically best-response profit | Unilateral gain | Maximum quote/BR discrepancy |
|:--|--:|--:|--:|--:|
| A | 3,131.860 | 3,131.860 | 0 | 3.04e-11 |
| B | 4,969.640 | 4,969.640 | 0 | 4.41e-10 |
| C | 7,505.362 | 7,505.362 | 0 | 4.69e-10 |

</div>

# 10. Environment and replication audit

The CSV exports in `tables/` are derived from the new run. They can be
used for independent figure checks, but **are not read by any
computation above**. No observed French insurer data or historic Figure
1–4 inputs are claimed to have been regenerated here. The model is
seeded identically in each run, but local SQP outputs may vary with
solver/library versions. Missing/failed guarded computations remain
marked as such; no archived values are inserted in their place.

<div id="tbl-manifest">

Table 12

**Table S12. Files produced by this render.**

| Items           | Count |
|:----------------|------:|
| Simulation runs |    21 |
| CSV tables      |    18 |
| PNG figures     |     7 |
| PDF figures     |     7 |

</div>

