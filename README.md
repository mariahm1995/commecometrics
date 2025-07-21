
<!-- README.md is generated from README.Rmd. Please edit that file -->

# commecometrics

<!-- badges: start -->

<!-- badges: end -->

The commecometrics package provides tools for conducting ecometric
analyses using trait data, species distributions, and environmental
variables. It supports workflows for both continuous (e.g.,
precipitation) and categorical (e.g., vegetation) environmental
variables, and includes utilities for fossil data integration,
paleoenvironmental reconstruction, and sensitivity analyses.

## Installation

You can install the development version of commecometrics from
[GitHub](https://github.com/) with:

``` r
# Install devtools if you don't have it
install.packages("devtools")

# Install commecometrics from GitHub
devtools::install_github("mariahm1995/commecometrics")
```

You can also get the official release version from CRAN

``` r
install.packages("commecometrics")
```

## Learn More

Two vignettes demonstrate the full workflow and capabilities of the
package:

- **Introduction to commecometrics**  
  Overview of trait–environment modeling, prediction, and visualization
  using modern data.

- **Reconstructing Paleoenvironments from Fossil Traits**  
  Extends the framework to fossil communities for estimating past
  environmental conditions.

You can also access them from within R:

``` r
# List all vignettes
browseVignettes("commecometrics")
#> No vignettes found by browseVignettes("commecometrics")

# Open a specific vignette
vignette("introduction-to-commecometrics")
#> Warning: vignette 'introduction-to-commecometrics' not found
```
