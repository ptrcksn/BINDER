
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

For more information, see corresponding
[publication](https://doi.org/10.1186/s12859-019-3042-8).

## BINDER R package

### Install BINDER

Perhaps the easiest way to install the `BINDER` R package from GitHub is
via the `install_github` function from the `devtools` (see `?devtools`
and `?devtools::install_github`) R
    package:

    if(!require("devtools")){ install.packages("devtools") } # Install devtools R package.
    devtools::install_github("ptrcksn/BINDER") # Install BINDER R package from GitHub.

This will install the most recent version of the `BINDER` R package.

(NOTE: The most recent iteration of `BINDER` uses a Gibbs sampling
approach (as opposed to the original HMC approach facilitated by Stan)
to estimate the posterior distribution of interest; in order to
facilitate the Gibbs sampling procedure, the scale parameters in the
current iteration assume independent inverse gamma prior distributions
(as opposed to independent truncated normal prior distributions as per
original approach) and the slope parameters in the auxiliary stratum
assume independent normal prior distributions (as opposed to independent
truncated normal prior distributions as per original approach). We find
little distinction between the implementations in terms of the estimated
posterior distribution; however, we find the current iteration is
considerably faster to run. If you desire the original HMC approach
facilitated by Stan, you can install this version by specifying this in
the call to
    `devtools::installgithub`:

    devtools::install_github("ptrcksn/BINDER", ref="bcfbec3") # Install HMC implementation of BINDER R package from GitHub.

If you are unsure as to which implementation to install, we recommend
ignoring this note and installing the most recent implementation as per
the initial installation instructions.)

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
represent target candiates for both putative regulators:

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
#> 1  Feature1         Feature1                      1  0  1
#> 2  Feature1         Feature2                      0  0  0
#> 3  Feature1         Feature3                      1  0  0
#> 4  Feature1         Feature4                      0  0  0
#> 5  Feature1         Feature5                      1  1  1
#> 6  Feature1         Feature6                      0  1  0
```

`expression` comprises an `N*M` matrix where each row corresponds to a
putative target feature of interest and each column corresponds to a
sample/experimental condition; element `n,m` corresponds to the
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
#>            Sample1    Sample2   Sample3   Sample4   Sample5
#> Feature1 3.5349622 3.14105065 0.3371845 0.5599447 2.9905482
#> Feature2 1.3354500 0.14327524 2.0137719 0.1678082 0.7546942
#> Feature3 1.0236223 0.08811973 0.7159417 0.3331889 0.3986848
#> Feature4 1.8661619 1.85259789 0.5513087 1.2355948 0.1033417
#> Feature5 0.2733524 0.36622549 2.6259289 0.2850141 0.4128588
```

### Run BINDER

Given these two data structures, the BINDER model implementation can be
invoked.

Although only the arguments `proxy_regulon` and `expression` require
user-specification, the `BINDER::binder` function also comprises several
function arguments in its function header that pertain to
hyperparameters, initial values for the MCMC sampler, thresholding
values, number of MCMC chains to run, number of draws from each MCMC
chain, extent of thinning and burn-in on the resulting MCMC chains,
parallelisation options and more. All function arguments except for
`proxy_regulon` and `expression` exhibit default behaviour in the
absence of user-specifcation; for details see `?BINDER::binder`.

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

The object returned by the call to `BINDER::binder` comprises a
two-dimensional list; each element of the top-level list corresponds to
a regulator, itself comprising a list of MCMC chains (of length equal to
`n_chains`); for a given regulator `r` and chain `c`, each chain
comprises a two-dimensional structure containing the MCMC draws (rows)
for each model parameter (columns) for chain `c` and regulator `r`.
Various results derived from the BINDER modelling process can be
accessed from the structure returned by the call to `BINDER::binder`.

The first 10 posterior draws for `theta` for {Feature1, Feature1},
{Feature1, Feature2} and {Feature1, Feature3} on chain 1 can be accessed
by:

``` r
results[["Feature1"]][[1]][1:10, 1:3]
#>       theta[Feature1] theta[Feature2] theta[Feature3]
#>  [1,]       0.1681374       0.3477076       0.6200443
#>  [2,]       0.1893283       0.3218594       0.4150221
#>  [3,]       0.2744142       0.4866501       0.3690196
#>  [4,]       0.5386341       0.4681985       0.4341110
#>  [5,]       0.5489026       0.4662690       0.4696044
#>  [6,]       0.5723213       0.5326793       0.4855688
#>  [7,]       0.5631092       0.6393360       0.5263637
#>  [8,]       0.5980026       0.5883900       0.4762252
#>  [9,]       0.5602296       0.5303741       0.5674310
#> [10,]       0.5826446       0.6010076       0.5746433
```

The 5th, 9th and 11th posterior draws for `theta` for {Feature2,
Feature11}, logit(`theta`) for {Feature2, Feature9}, `zeta` and `tau_ME`
on chain 1 can be accessed
by:

``` r
results$Feature2[[1]][c(5, 9, 11), c("theta[Feature11]", "logit(theta[Feature9])", "zeta", "tau_ME")]
#>      theta[Feature11] logit(theta[Feature9])        zeta     tau_ME
#> [1,]        0.4965259              0.1295014 -0.08838271 0.15016206
#> [2,]        0.5802706              0.2583725  0.25203249 0.03674674
#> [3,]        0.5652968              0.2808059  0.28476981 0.01532989
```

You can also access summaries (mean, standard deviation and 0%, 25%,
50%, 75% and 100% quantiles) pertaining to the posterior distribution of
`theta` for each regulator-target candidate pair of interest using the
custom `summary` method with an instance of `posteriorDraws` (such as
that returned by `BINDER::binder`):

``` r
head(summary(results))
#>   regulator target_candidate      mean   std.dev.        0%       25%
#> 1  Feature1         Feature1 0.5928781 0.03178056 0.1681374 0.5788284
#> 2  Feature1         Feature2 0.5778794 0.02750820 0.3218594 0.5625920
#> 3  Feature1         Feature3 0.5670677 0.02586060 0.3690196 0.5512209
#> 4  Feature1         Feature4 0.5815275 0.02680061 0.3458604 0.5660052
#> 5  Feature1         Feature5 0.5841839 0.03266487 0.1612005 0.5695483
#> 6  Feature1         Feature6 0.5655742 0.02804707 0.2737798 0.5503719
#>         50%       75%      100%
#> 1 0.5946083 0.6099758 0.6667352
#> 2 0.5786746 0.5952952 0.6580838
#> 3 0.5670795 0.5831421 0.6405395
#> 4 0.5814419 0.5979161 0.6492646
#> 5 0.5855834 0.6015050 0.6622406
#> 6 0.5675410 0.5827345 0.6378946
```

The default behaviour of `summary` is to provide posterior summaries for
all regulator-target candidate pairs; however, if you are interested in
a particular set of summaries, you can specify the the regulator-target
canididates through the `target_candidates` function
argument:

``` r
relevant_target_candidates <- list("Feature1"=c("Feature1", "Feature7"), "Feature2"=c("Feature11", "Feature15"))
summary(results, target_candidates=relevant_target_candidates)
#>   regulator target_candidate      mean   std.dev.        0%       25%
#> 1  Feature1         Feature1 0.5928781 0.03178056 0.1681374 0.5788284
#> 2  Feature1         Feature7 0.5624138 0.02499361 0.3595835 0.5462711
#> 3  Feature2        Feature11 0.5610081 0.02414697 0.2610603 0.5469466
#> 4  Feature2        Feature15 0.5703521 0.02483849 0.2536209 0.5560276
#>         50%       75%      100%
#> 1 0.5946083 0.6099758 0.6667352
#> 2 0.5628532 0.5776084 0.6499368
#> 3 0.5611947 0.5755981 0.6524184
#> 4 0.5704219 0.5845328 0.6950577
```
