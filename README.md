
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BINDER

<!-- badges: start -->

<!-- badges: end -->

## Introduction

BINDER (BayesIan gene regulatory Networks inferreD via gene coExpression
and compaRative genomics) is a hybrid approach, integrating coexpression
and comparative genomics to infer prokaryotic regulons.

BINDER infers gene regulatory networks by constructing two strata: a
primary stratum and an auxiliary stratum. The primary stratum is
composed of coexpression data from a primary organism of interest and is
supplemented by auxiliary data in the form of motif predictions and
known orthologous regulator-target interactions from a proxy organism.
BINDER implements a Bayesian hierarchical model that appositely models
the type and structure of both this primary and auxiliary data to infer
the probability of a regulatory interaction between a regulator-target
pair. The auxiliary data inform the prior distributions and the
posterior distributions are updated by accounting for the primary
coexpression data in a novel, apposite bivariate likelihood function.

## BINDER R package

### Install BINDER

Perhaps the easiest way to install the `BINDER` R package from GitHub is
via the `install_github` function from the `devtools` (see `?devtools`
and `?devtools::install_github`) R
    package:

    if(!require("devtools")){ install.packages("devtools") } # Install devtools R package.
    devtools::install_github("ptrcksn/BINDER") # Install BINDER R package from GitHub.

### Simulate some data

Let’s define our data such that we have `N` candidate features of
interest and `M` samples/experimental conditions:

``` r
N <- 500 # Number of putative target features of interest.
M <- 100 # Number of samples/experimental conditions.
```

The main function of the `BINDER` package is the `BINDER::binder`
function; this function runs the BINDER model implementation with the
data provided.

The `BINDER::binder` function comprises several function arguments (see
`?BINDER::binder`); two of these arguments have no default values and
require user-defined values: `proxy_regulon` and `expression`.

`proxy_regulon` requires a data structure of class `data.frame` with the
following columns:

\[`regulator` `target_candidate` `ortholog_module_status` `ME` `PE`\]

Each column of `proxy_regulon` should be named and ordered as above.

Let’s simulate a `proxy_regulon` structure where Feature1 and Feature2
represent putative regulators and Feature1, Feature2, …, Feature500
represent target candiates for both putative regulators.

``` r
# Simulate proxy regulon data frame:
proxy_regulon <- data.frame(
  regulator=rep(paste0("Feature", 1:2), each=N),
  target_candidate=rep(paste0("Feature", 1:N), 2),
  ortholog_module_status=rbinom((2*N), 1, 0.5),
  ME=rbinom((2*N), 1, 0.5),
  PE=rbinom((2*N), 1, 0.5))
proxy_regulon$PE <- ifelse(proxy_regulon$ortholog_module_status == 0, 0, proxy_regulon$PE)
head(proxy_regulon)
#>   regulator target_candidate ortholog_module_status ME PE
#> 1  Feature1         Feature1                      0  0  0
#> 2  Feature1         Feature2                      1  0  1
#> 3  Feature1         Feature3                      0  1  0
#> 4  Feature1         Feature4                      1  0  1
#> 5  Feature1         Feature5                      1  0  0
#> 6  Feature1         Feature6                      0  0  0
```

`expression` comprises an `N*M` matrix where each row corresponds to a
putative target feature of interest and each column corresponds to a
sample/experimental condition; the `n,m`th element corresponds to the
expression value for the `n`th feature of interest under the `m`th
sample/experimental condition. Each row should be named by its
associated feature.

Let’s simulate an `expression` structure for each of the `N` features of
interest across the `M` samples/experimental conditions:

``` r
# Simulate expression matrix:
expression <- matrix(rgamma((N*M), 1, 1), nrow=N) # Simulate expression data.
rownames(expression) <- paste0("Feature", 1:N) # Names of features of interest.
colnames(expression) <- paste0("Sample", 1:M) # Names of samples/experimental conditions.
expression[1:5, 1:5]
#>            Sample1    Sample2    Sample3   Sample4   Sample5
#> Feature1 0.9419525 1.86626376 0.83296264 1.2543381 1.6411305
#> Feature2 0.6543574 1.26581845 0.19091606 0.5420611 1.4454896
#> Feature3 4.0205033 0.04930548 2.37388265 0.5261245 1.8328811
#> Feature4 1.0204354 0.43533193 0.80722164 0.2643896 0.4373572
#> Feature5 1.0671127 0.65192464 0.09135318 1.1526399 2.5746012
```

### Run BINDER

Given these two data structures, the BINDER model implementation can be
invoked; let’s run the `BINDER::binder` function on the simulated
`proxy_regulon` and `expression` data
    structures:

    results <- BINDER::binder(proxy_regulon=proxy_regulon, expression=expression)

The binder function argument `is.coexpression` defaults to `FALSE`; this
means that BINDER treats the expression matrix supplied to the function
as an `N*M` matrix where each row corresponds to a feature of interest
and each column corresponds to a sample/experimental condition.

(If a coexpression matrix is supplied to the BINDER function, the
`is.coexpression` function argument should be set to `TRUE`. A
coexpression matrix can be computed from an expression matrix using the
`BINDER::compute_coexpression` function (see
`?BINDER::compute_coexpression`).)

### Access BINDER results

The object returned by the call to `BINDER::binder` comprises a list of
lists; each element of the top-level list corresponds to a regulator,
itself comprising a list of model results as pertains to that regulator.
Various results derived from the BINDER modelling process can be
accessed from the list of lists returned by the call to
`BINDER::binder`.

Posterior mean for `theta` for the first 5 target candidates pertaining
to the putative regulator Feature1:

``` r
results[["Feature1"]]$mean_theta[1:5]
#>  Feature1  Feature2  Feature3  Feature4  Feature5 
#> 0.2416932 0.2584999 0.2533846 0.2581473 0.2569062
```

Posterior 0.025%-0.975% credible interval for `zeta` for the putative
regulator Feature2:

``` r
results$Feature2$zeta_interval
#>     0.025     0.975 
#> -8.944301 -3.077080
```

Similarly, the S4 class stanfit objects (see `?rstan::sampling`)
representing the fitted results for each regulator model can also be
inspected from via the object returned by the `BINDER::binder` function
call; here are the first 3 `Rhat` values derived from `stanfit` objetct
for Feature1:

``` r
rstan::summary(results[["Feature1"]]$model_object)$summary[1:3, "Rhat"]
#>     zeta   tau[1]   tau[2] 
#> 1.001776 1.000745 1.000233
```
