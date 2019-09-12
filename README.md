
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
#> 1  Feature1         Feature1                      0  1  0
#> 2  Feature1         Feature2                      0  0  0
#> 3  Feature1         Feature3                      1  0  0
#> 4  Feature1         Feature4                      1  0  1
#> 5  Feature1         Feature5                      0  1  0
#> 6  Feature1         Feature6                      0  0  0
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
#>            Sample1    Sample2   Sample3   Sample4    Sample5
#> Feature1 2.3489117 0.35507667 0.1508703 0.1016637 1.62519006
#> Feature2 0.1918560 1.25832293 0.3032065 0.5291326 2.87581919
#> Feature3 0.1183431 0.85321436 0.5866373 2.1859207 0.05352903
#> Feature4 1.0078704 0.06905414 0.1479612 1.0205423 0.13920528
#> Feature5 0.2869065 0.05668388 1.2943560 0.1759438 0.53792276
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
`is.coexpression` function argument should be set to `TRUE`. If desired,
a coexpression matrix can be computed from an expression matrix using
the `BINDER::compute_coexpression` function (see
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
#>  [1,]       0.3444320       0.3474237      0.73338267
#>  [2,]       0.2932088       0.2418452      0.38029164
#>  [3,]       0.2947132       0.3212269      0.21580408
#>  [4,]       0.5383204       0.2216253      0.20809058
#>  [5,]       0.4237041       0.1628518      0.18740706
#>  [6,]       0.3749246       0.1871047      0.14907236
#>  [7,]       0.2768919       0.3940304      0.16958468
#>  [8,]       0.3602780       0.2560915      0.08207835
#>  [9,]       0.2170996       0.1190923      0.22465461
#> [10,]       0.2528819       0.2751432      0.24524074
```

The 5th, 9th and 11th posterior draws for `theta` for {Feature2,
Feature11}, logit(`theta`) for {Feature2, Feature9}, `zeta` and `tau_ME`
on chain 1 can be accessed
by:

``` r
results$Feature2[[1]][c(5, 9, 11), c("theta[Feature11]", "logit(theta[Feature9])", "zeta", "tau_ME")]
#>      theta[Feature11] logit(theta[Feature9])      zeta     tau_ME
#> [1,]        0.2576879             -0.7543109 -1.071112 0.08909659
#> [2,]        0.3375324             -1.2019174 -1.187885 0.10538952
#> [3,]        0.2748118             -1.1587853 -1.267864 0.14843953
```

You can also access summaries (mean, standard deviation and 0%, 25%,
50%, 75% and 100% quantiles) pertaining to the posterior distribution of
`theta` for each regulator-target candidate pair of interest using the
custom `summary` method with an instance of `posteriorDraws` (such as
that returned by `BINDER::binder`):

``` r
head(summary(results))
#>   regulator target_candidate      mean   std.dev.         0%       25%
#> 1  Feature1         Feature1 0.2749750 0.06839296 0.06246811 0.2270924
#> 2  Feature1         Feature2 0.2162535 0.06132870 0.06998652 0.1738338
#> 3  Feature1         Feature3 0.2257045 0.06387967 0.08207835 0.1801842
#> 4  Feature1         Feature4 0.2377410 0.06420592 0.09289226 0.1929313
#> 5  Feature1         Feature5 0.1732752 0.05070386 0.05180642 0.1362150
#> 6  Feature1         Feature6 0.2344679 0.06260424 0.08615899 0.1885697
#>         50%       75%      100%
#> 1 0.2704985 0.3169951 0.5383204
#> 2 0.2091073 0.2544163 0.4599965
#> 3 0.2180709 0.2620565 0.7333827
#> 4 0.2302281 0.2777345 0.4568322
#> 5 0.1679466 0.2038884 0.4005127
#> 6 0.2308776 0.2734669 0.4674663
```

The default behaviour of `summary` is to provide posterior summaries for
all regulator-target candidate pairs; however, if you are interested in
a particular set of summaries, you can specify the the regulator-target
canididates through the `target_candidates` function
argument:

``` r
relevant_target_candidates <- list("Feature1"=c("Feature1", "Feature7"), "Feature2"=c("Feature11", "Feature15"))
summary(results, target_candidates=relevant_target_candidates)
#>   regulator target_candidate      mean   std.dev.         0%       25%
#> 1  Feature1         Feature1 0.2749750 0.06839296 0.06246811 0.2270924
#> 2  Feature1         Feature7 0.1850818 0.05277211 0.07245370 0.1474260
#> 3  Feature2        Feature11 0.2422060 0.06548180 0.08439559 0.1964373
#> 4  Feature2        Feature15 0.1923132 0.05621942 0.06748995 0.1541940
#>         50%       75%      100%
#> 1 0.2704985 0.3169951 0.5383204
#> 2 0.1803393 0.2140034 0.4640942
#> 3 0.2354611 0.2797723 0.5416652
#> 4 0.1850068 0.2237455 0.5165300
```
