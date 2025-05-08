
<!-- README.md is generated from README.Rmd. Please edit that file -->

# commecometrics

<!-- badges: start -->
<!-- badges: end -->

commecometrics provides tools for conducting trait-based ecometric
modeling with modern and fossil data. It supports workflows for both
continuous and categorical environmental variables, as well as
sensitivity and transferability analysis.

## Installation

You can install the development version of commecometrics from
[GitHub](https://github.com/) with:

``` r
# Install devtools if you don't have it
install.packages("devtools")
#> Installing package into '/private/var/folders/4g/l0bfkncj5vvgbtq2lccc64900000gn/T/RtmpJy4czF/temp_libpathfabc29e803e1'
#> (as 'lib' is unspecified)
#> 
#> The downloaded binary packages are in
#>  /var/folders/4g/l0bfkncj5vvgbtq2lccc64900000gn/T//RtmpZQxLsB/downloaded_packages

# Install commecometrics from GitHub
devtools::install_github("mariahm1995/commecometrics")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo mariahm1995/commecometrics@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/4g/l0bfkncj5vvgbtq2lccc64900000gn/T/RtmpZQxLsB/remotes1015476d4631a/mariahm1995-commecometrics-615454284f9163f3f744c2e94ce6e9522838e5ce/DESCRIPTION’ ... OK
#> * preparing ‘commecometrics’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘commecometrics_0.0.0.9000.tar.gz’
#> Installing package into '/private/var/folders/4g/l0bfkncj5vvgbtq2lccc64900000gn/T/RtmpJy4czF/temp_libpathfabc29e803e1'
#> (as 'lib' is unspecified)
```

## Example

This is a basic example of how use commecometrics

``` r
library(commecometrics)

# Load example data
data("points", package = "commecometrics")
data("traits", package = "commecometrics")
data("polygons", package = "commecometrics")

# Summarize trait values at points
traitsByPoint <- summarize_traits_by_point(
  points_df = points,
  trait_df = traits,
  species_polygons = polygons,
  trait_column = "RBL",
  species_name_col = "sci_name"
)

# Run ecometric model
ecoModel <- ecometric_model(
  points_df = traitsByPoint$points,
  env_var = "BIO12",
  transform_fun = function(x) log(x + 1),
  inv_transform_fun = function(x) exp(x) - 1,
  min_species = 3
)

# Plot trait–environment surface
ecoPlot <- ecometric_space(ecoModel, env_name = "Precipitation (loge mm)")
print(ecoPlot)
```

## Learn More

See the vignette for a full walkthrough:

``` r
vignette("Introduction to commecometrics")
#> Warning: vignette 'Introduction to commecometrics' not found
```
