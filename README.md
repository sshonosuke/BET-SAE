# Bayesian Benchmarking Small Area Estimation via Entropic Tilting

R code for the simulation studies in:

Shonosuke Sugasawa, Genya Kobayashi, and Yuki Kawakubo.  
**Bayesian Benchmarking Small Area Estimation via Entropic Tilting.**  
https://arxiv.org/abs/2407.17848

The code implements Bayesian benchmarking via entropic tilting and compares it with alternative benchmarking methods under Fay–Herriot and Poisson–gamma models.

## Files

- `FH-function.R`: Estimation and benchmarking functions for the Fay–Herriot model with normal and horseshoe random effects.
- `PG-function.R`: Estimation and benchmarking functions for the Poisson–gamma model.
- `Sim-FH-A.R` and `Sim-FH-B.R`: Fay–Herriot simulations under settings A and B.
- `Sim-PG-A.R` and `Sim-PG-B.R`: Poisson–gamma simulations under settings A and B.

## Requirements

Install the required R packages:

```r
install.packages(c("MASS", "MCMCpack", "statmod", "progress"))
```

## Usage

Place all files in the same directory and set it as the R working directory. Run a simulation script, for example:

```r
source("Sim-FH-A.R")
```

Each script allows the data-generating scenario (`scenario=1, 2, 3`), number of Monte Carlo replications (`R`), and number of areas (`m`) to be specified near the beginning. The defaults are `scenario=3`, `R=2000`, and `m=30`. For a quick trial, reduce `R`.

Results, including mean squared errors, coverage probabilities, and average interval lengths, are saved as `.RData` files in the working directory.
