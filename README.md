
# ecrins

<!-- badges: start -->
[![R-CMD-check](https://github.com/petedodd/ecrins/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/petedodd/ecrins/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of ecrins is to simulate PPD flows and TB progression with ordinary differential equations.

## Installation

If your R set up is able to compile packages, you can install the development version of ecrins from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("petedodd/ecrins")
```


## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ecrins)

data(parms)                       #some default parameters
tt <- seq(from=0, to=10, by=0.1)  #time frame to run over
y <- runmodel(tt,parms)           #run model
head(y)                           #look at answer

```

