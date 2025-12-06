# CXCL13 Deconvolution Simulations for TIRTL-seq

This repository contains R code used to simulate and evaluate CXCL13 expression deconvolution from TIRTL-seq style plate data. The code implements a weighted non-negative least squares (NNLS) framework to infer clonotype-level CXCL13 expression from well-level mixture measurements, as described in the accompanying internship report.

The focus is on understanding identifiability limits under realistic TIRTL loading conditions (number of wells, clonotypes per well, and read depth heterogeneity).

---

## Mathematical model

When running TIRTL at high density (e.g. full 384-well plates), each well contains a mixture of clonotypes with known relative abundances (fractions or normalized read counts). If we additionally measure CXCL13 expression per well, we can model the well-level signal as a linear combination of clonotype-level expression.

Let

- \( y_w \): CXCL13 read count (or normalized abundance) in well \( w \)
- \( A_{w,c} \): fraction or estimated cell count of clonotype \( c \) in well \( w \)
- \( x_c \): CXCL13 expression per clonotype (unknown)
- \( s_w \): optional per-well scale factor (e.g. based on total reads)
- \( b \): background term
- \( \varepsilon_w \): noise

### Linear mixture model

For each well \( w \):

\[
y_w \approx s_w \sum_{c=1}^{C} A_{w,c} x_c + b + \varepsilon_w.
\]

In the simplest case we set \( s_w = 1 \) and solve

\[
y_w \approx \sum_{c=1}^{C} A_{w,c} x_c + b.
\]

This is a standard non-negative regression problem where each well contributes one equation.


### Weighted NNLS

To account for differences in sequencing depth across wells, we use weighted NNLS:

\[
\min_{x \ge 0} \sum_{w=1}^{W} w_w \left( y_w - (A x)_w \right)^2,
\]

where the weights are derived from well-level TIRTL read counts \( r_w \):

\[
w_w = \left( \frac{r_w}{\operatorname{median}(r)} \right)^{1/2},
\]

with clipping to avoid extreme dominance by a few wells:

\[
10^{-3} \le w_w \le 10^{3}.
\]

The solution vector \( x_c \) provides per-clonotype CXCL13 expression, while the fitted background \( b \) captures plate-level baseline signal.



### Simulation framework

Ground-truth clonotype-level CXCL13 values are simulated on the log\(_2(\text{CXCL13} + 1)\) scale using two distributions:

- **Non-reactive clonotypes:** low-intensity distribution (e.g. mean = 0.3, SD = 0.4).
- **Reactive clonotypes:** higher-intensity distribution (e.g. mean = 3.5, SD = 0.7).

These parameters are motivated by bimodal CXCL13 expression patterns observed in human T cell scRNA-seq datasets.

Simulations then:

1. Use an observed TIRTL plate layout (clonotypes per well) as the mixing matrix \( A \).
2. Assign ground-truth \( x_c \) from the two log-normal distributions.
3. Optionally add:
   - row-level multiplicative noise (biological variability),
   - well-level Gaussian measurement noise at the aggregated \( y_w \) level.
4. Fit the NNLS model and compare estimated vs true \( x_c \) via correlation and recall.

Additional experiments vary:

- the expression gap between reactive and non-reactive states,
- the number of clonotypes per well,
- the number of wells used in the fit,

to explore identifiability limits of the deconvolution problem.

---



## Repository structure

- `main_simulation_logic/`
  - `Ground_Truth_Generator.R`  
    Core script to construct the ground-truth clonotype-level CXCL13 expression and the corresponding well-level mixtures based on a TIRTL-style plate layout.
  - `Simulating_and_Fitting.R`  
    Main simulation and analysis pipeline. Runs repeated simulations, applies weighted NNLS, and computes evaluation metrics (e.g. correlations, recall vs expression gap, dependence on clonotype density and well count).

- `testing_scripts/`
  - R scripts used for exploratory tests, parameter sweeps, and sanity checks referenced in the report (e.g. testing different ranges of clonotypes per well or numbers of wells). These are not required for the main pipeline but are included for transparency and reproducibility.

---

## Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/echo-logist/cxcl13-tirtl-simulations.git
   cd cxcl13-tirtl-simulations
