
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
#> 1  Feature1         Feature1                      1  1  0
#> 2  Feature1         Feature2                      1  0  0
#> 3  Feature1         Feature3                      1  1  0
#> 4  Feature1         Feature4                      0  1  0
#> 5  Feature1         Feature5                      0  0  0
#> 6  Feature1         Feature6                      1  0  1
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
#>            Sample1   Sample2   Sample3    Sample4    Sample5
#> Feature1 0.2604576 0.1690360 1.0320920 1.27054736 0.78125372
#> Feature2 3.5840574 0.5195996 0.2093035 0.01030329 0.40737647
#> Feature3 2.4036786 0.5116417 2.7342324 0.77338070 1.83508608
#> Feature4 2.7776564 0.5551269 1.0838316 1.25154840 0.38639540
#> Feature5 2.7821946 0.6652486 1.1292559 1.53459183 0.08226841
```

### Run BINDER

Given these two data structures, the BINDER model implementation can be
invoked.

If you want to run `BINDER::binder` in parallel, you can alter the
`mc.cores` global option; for example, to run `BINDER::binder` on 3
cores, you can specify the following (see `?options` and `?parallel`):

``` r
options(mc.cores=3)
```

Let’s run the `BINDER::binder` function on the simulated `proxy_regulon`
and `expression` data
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
#> 0.2471730 0.2235725 0.2388696 0.2385041 0.2304853
```

Posterior 0.025%-0.975% credible interval for `zeta` for the putative
regulator Feature2:

``` r
results$Feature2$zeta_interval
#>     0.025     0.975 
#> -9.638953 -3.631587
```

Similarly, the S4 class stanfit objects (see `?rstan::sampling`)
representing the fitted results for each regulator model can also be
inspected from via the object returned by the `BINDER::binder` function
call; here are the first 3 `Rhat` values derived from `stanfit` objetct
for Feature1:

``` r
rstan::summary(results[["Feature1"]]$model_object)$summary[1:3, "Rhat"]
#>      zeta    tau[1]    tau[2] 
#> 1.0007600 0.9996117 1.0003351
```
