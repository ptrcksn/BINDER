
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
#> 2  Feature1         Feature2                      1  1  1
#> 3  Feature1         Feature3                      1  1  1
#> 4  Feature1         Feature4                      0  0  0
#> 5  Feature1         Feature5                      0  1  0
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
#>            Sample1   Sample2   Sample3    Sample4    Sample5
#> Feature1 0.2985589 0.4683038 1.2125910 0.20447312 0.37017539
#> Feature2 0.2278294 0.1124240 0.4461277 0.24682153 0.08901655
#> Feature3 2.4278190 0.6272579 0.5405478 2.47196959 0.63504147
#> Feature4 1.1062624 0.3100832 0.6648737 0.06052497 0.28580735
#> Feature5 0.6507402 0.1423443 0.9094857 0.44851433 1.75663189
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
#>  [1,]       0.3590215       0.2173603       0.4595265
#>  [2,]       0.3092595       0.2184844       0.3053142
#>  [3,]       0.3811403       0.4171059       0.3186562
#>  [4,]       0.6072865       0.4564324       0.4384101
#>  [5,]       0.5659700       0.4628380       0.4865641
#>  [6,]       0.5659389       0.5187316       0.4967011
#>  [7,]       0.5379224       0.6207171       0.5423335
#>  [8,]       0.5679212       0.5721556       0.4961128
#>  [9,]       0.5262889       0.5239661       0.5797066
#> [10,]       0.5369010       0.5814917       0.5845185
```

The 5th, 9th and 11th posterior draws for `theta` for {Feature2,
Feature11}, logit(`theta`) for {Feature2, Feature9}, `zeta` and `tau_ME`
on chain 1 can be accessed
by:

``` r
results$Feature2[[1]][c(5, 9, 11), c("theta[Feature11]", "logit(theta[Feature9])", "zeta", "tau_ME")]
#>      theta[Feature11] logit(theta[Feature9])        zeta     tau_ME
#> [1,]        0.4829353             0.03733617 -0.08612575 0.09369568
#> [2,]        0.6120259             0.30166177  0.23387067 0.05668823
#> [3,]        0.5959564             0.35355580  0.28216806 0.01364373
```

You can also access summaries (mean, standard deviation and 0%, 25%,
50%, 75% and 100% quantiles) pertaining to the posterior distribution of
`theta` for each regulator-target candidate pair of interest using the
custom `summary` method with an instance of `posteriorDraws` (such as
that returned by `BINDER::binder`):

``` r
head(summary(results))
#>   regulator target_candidate      mean   std.dev.        0%       25%
#> 1  Feature1         Feature1 0.5517206 0.02377518 0.3092595 0.5386721
#> 2  Feature1         Feature2 0.5557324 0.02730216 0.2173603 0.5421914
#> 3  Feature1         Feature3 0.5704321 0.02479461 0.3053142 0.5571993
#> 4  Feature1         Feature4 0.5798860 0.02460916 0.3272047 0.5665950
#> 5  Feature1         Feature5 0.5493962 0.02507074 0.2869857 0.5358284
#> 6  Feature1         Feature6 0.5980025 0.02627707 0.2737798 0.5847863
#>         50%       75%      100%
#> 1 0.5534270 0.5661216 0.6147910
#> 2 0.5557368 0.5715291 0.6286938
#> 3 0.5712254 0.5845883 0.6346843
#> 4 0.5798895 0.5946995 0.6424778
#> 5 0.5504214 0.5638226 0.6290217
#> 6 0.6000122 0.6132843 0.6595397
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
#> 1  Feature1         Feature1 0.5517206 0.02377518 0.3092595 0.5386721
#> 2  Feature1         Feature7 0.5541087 0.02233644 0.3595835 0.5402106
#> 3  Feature2        Feature11 0.5892315 0.02900208 0.1328382 0.5753144
#> 4  Feature2        Feature15 0.5715759 0.02823166 0.1463963 0.5582078
#>         50%       75%      100%
#> 1 0.5534270 0.5661216 0.6147910
#> 2 0.5545264 0.5677760 0.6328728
#> 3 0.5901010 0.6048001 0.6572572
#> 4 0.5719999 0.5868540 0.6419489
```
