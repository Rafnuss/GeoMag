
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GeoMag

<!-- badges: start -->

[![R-CMD-check](https://github.com/Rafnuss/GeoMag/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Rafnuss/GeoMag/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Rafnuss/GeoMag/graph/badge.svg)](https://app.codecov.io/gh/Rafnuss/GeoMag)
[![lint](https://github.com/Rafnuss/GeoMag/actions/workflows/lint.yaml/badge.svg)](https://github.com/Rafnuss/GeoMag/actions/workflows/lint.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

**GeoMag** is an R package to estimate animal geolocation using magnetic
field data.  
It provides tools for robust calibration, outlier removal, likelihood
map computation, and interactive 3D visualization of magnetic and
acceleration data.  
GeoMag is designed to work seamlessly with
[GeoPressureR](https://raphaelnussbaumer.com/GeoPressureR), enabling
high-resolution migratory track reconstruction using multi-sensor
archival tags.

------------------------------------------------------------------------

## 🚀 Overview

GeoMag implements a robust workflow for geolocation based on magnetic
field measurements, including:

- **Calibration:** Correct sensor biases, scale factors, and rotations
  using in situ or external datasets.
- **Outlier removal:** Automatically detect and exclude spurious
  measurements.
- **Likelihood mapping:** Compute spatial likelihood maps using the
  World Magnetic Model (WMM), comparing observed and predicted intensity
  and inclination.
- **Interactive visualization:** Explore raw and calibrated data,
  ellipsoid fits, and 3D scatterplots.

------------------------------------------------------------------------

## 🦅 Main Features

- **Magnetic sensor calibration** (including ellipsoid and sphere
  fitting)
- **Likelihood mapping** of stationary periods (`geomag_map`)
- **Integration with
  [GeoPressureR](https://raphaelnussbaumer.com/GeoPressureR) tag
  analysis workflow**
- **Interactive 3D plotting** of sensor data (`plot_mag`)
- **Support for WMM-based reference maps**

------------------------------------------------------------------------

## 📦 Installation

To install the latest version from GitHub:

``` r
# install.packages("pak")
pak::pkg_install("Rafnuss/GeoMag")
```

------------------------------------------------------------------------

## 📖 Learn More

- The
  [GeoPressureManual](https://raphaelnussbaumer.com/GeoPressureManual/)
  includes a dedicated section on magnetic-based geolocation and
  demonstrates the GeoMag workflow with real data.
- See also the [GeoPressureR
  reference](https://raphaelnussbaumer.com/GeoPressureR/reference/index.html)
  for tag object structure and integration.

------------------------------------------------------------------------

## 🛠️ Example Usage

``` r
library(GeoMag)
library(GeoPressureR)

library(GeoPressureR)
withr::with_dir(system.file("extdata", package = "GeoMag"), {
  # Create a GeoPressureR tag object (see GeoPressureR documentation)
  tag <- tag_create("14DM")

  # Label the tag
  tag <- tag_label(tag)
})

# Calibrate the tag's magnetic data
tag <- geomag_calib(tag)

# Interactive 3D plot of calibrated magnetic field
plot_mag(tag, type = "acceleration")
plot_mag(tag, type = "magnetic")
plot_mag(tag, type = "calib")

# Compute the spatial likelihood map for each stationary period
tag <- geomag_map(tag)

plot(tag, "map_magnetic")
```

------------------------------------------------------------------------

## 📚 Citation

If you use GeoMag in your research, please cite:

> Nussbaumer, R. (2025). GeoMag: Magnetic Field-Based Geolocation in R.
> <https://github.com/Rafnuss/GeoMag>

For citation information in R:

``` r
citation("GeoMag")
```

------------------------------------------------------------------------

## 🔗 Related Projects

- [GeoPressureR](https://raphaelnussbaumer.com/GeoPressureR):
  Pressure-based geolocation
- [GeoLocatoR](https://github.com/Rafnuss/GeoLocatoR): Data packaging
  for geolocator datasets
- [GeoPressureTemplate](https://github.com/Rafnuss/GeoPressureTemplate):
  Template repository for geolocator analysis

------------------------------------------------------------------------

## 📝 License

This package is licensed under [CC BY
4.0](https://creativecommons.org/licenses/by/4.0/).
