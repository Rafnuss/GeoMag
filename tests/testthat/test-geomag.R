library(GeoPressureR)
library(GeoMag)
library(testthat)

setwd(system.file("extdata", package = "GeoMag"))

tag <- GeoPressureR::tag_create("14DM", quiet = TRUE)

plot_mag(tag, "acceleration")
expect_error(plot_mag(tag, "acceleration_p"))
plot_mag(tag, "magnetic")
expect_error(plot_mag(tag, "calib"))
expect_error(plot_mag(tag, "timeseries"))
expect_error(plot_mag(tag, "histogram"))


# Calib without stap: ellipse
tag2 <- geomag_calib(tag, quiet = TRUE)
expect_true(all(c("magnetic_xc", "magnetic_yc", "magnetic_zc") %in% names(tag2$magnetic)))
expect_true("is_static" %in% names(tag2$magnetic))
expect_true(tag2$param$geomag_calib$calib_method == "ellipse")
plot_mag(tag2, "magnetic")
plot_mag(tag2, "calib")
plot_mag(tag2, "acceleration_p")
plot_mag(tag2, "timeseries")
expect_error(plot_mag(tag2, "histogram"))

# Calib with stap
tag <- tag_label(tag, quiet = TRUE)
# tag2 <- geomag_calib(tag, quiet = TRUE)
# plot_mag(tag2,"calib")

# Calib with known
tag <- tag_set_map(tag,
  extent = c(-18, 23, 0, 50),
  scale = 2,
  known = data.frame(
    stap_id = c(1, -1),
    known_lon = 7.27,
    known_lat = 46.19
  )
)
tag <- geomag_calib(tag, quiet = TRUE)
plot_mag(tag, "magnetic")
plot_mag(tag, "calib")
plot_mag(tag, "acceleration_p")
plot_mag(tag, "timeseries")
expect_error(plot_mag(tag2, "histogram"))

# GeoMag
tag <- geomag_map(tag, quiet = TRUE)

plot(tag, "map_magnetic")

# Plot
path = tag2path(tag)
plot_mag(tag, "timeseries", path=path)
plot_mag(tag2, "histogram", path=path)
