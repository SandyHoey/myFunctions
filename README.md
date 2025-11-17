# myFunctions
a package that contains all the random functions I write for myself to save time

## Installation

You can install the development version of pets from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("SandyHoey/myFunctions")
```

## List of functions
#### boots_param_CI()
a function to bootstrap confidence intervals for model parameters from lme4 and glmmTMB

currently works with
- Poisson
- negative binomial
- binomial
- zero inflated Poisson
- zero inflated negative binomial 
