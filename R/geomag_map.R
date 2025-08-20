#' Compute Magnetic Likelihood Maps
#'
#' @description
#' Estimates a likelihood map for each stationary period based on observed magnetic data
#' (intensity and inclination) compared to the expected values from the World Magnetic Model (WMM).
#'
#' @param tag A [GeoPressureR tag object
#'  ](https://raphaelnussbaumer.com/GeoPressureR/reference/tag_create.html).
#' @param compute_known Logical. If TRUE, computes likelihood maps for known stationary periods;
#'   if FALSE, fixes the likelihood at the known location.
#' @param sd_e_f Numeric. Standard deviation of observation noise of intensity error (single value
#' or one per stap).
#' @param sd_e_i Numeric. Standard deviation of observation noise of inclination error (single
#' value or one per stap).
#' @param sd_m_f Numeric. Standard deviation of stationary-period-specific noise of intensity error
#' (single value or per stap).
#' @param sd_m_i Numeric. Standard deviation of stationary-period-specific noise of inclination
#' error (single value or per stap).
#' @param ref_map Raster stack or list. Reference magnetic maps (intensity and inclination)
#' computed by default using `geomag_map_ref()`. Provide it if you've already computed it to save
#' computational time.
#' @param quiet Logical. If TRUE, suppresses progress messages.
#'
#' @return A [GeoPressureR tag object
#' ](https://raphaelnussbaumer.com/GeoPressureR/reference/tag_create.html) with likelihood maps
#' added as:
#'   - `tag$map_magnetic_intensity`
#'   - `tag$map_magnetic_inclination`
#'   - `tag$map_magnetic`
#'   and updated parameters in `tag$param`.
#'
#' @examples
#' library(GeoPressureR)
#' withr::with_dir(system.file("extdata", package = "GeoMag"), {
#'   tag <- tag_create("14DM", quiet = TRUE)
#'   tag <- tag_label(tag, quiet = TRUE)
#'   tag <- tag_set_map(tag,
#'     extent = c(-18, 23, 0, 50),
#'     scale = 2,
#'     known = data.frame(
#'       stap_id = c(1, -1),
#'       known_lon = 7.27,
#'       known_lat = 46.19
#'     )
#'   )
#'   tag <- geomag_calib(tag, quiet = TRUE)
#'   tag <- geomag_map(tag, quiet = TRUE)
#'   plot(tag, type = "map_magnetic_intensity")
#'   plot(tag, type = "map_magnetic_inclination")
#' })
#' @export
geomag_map <- function(tag,
                       compute_known = FALSE,
                       sd_e_f = .007,
                       sd_e_i = 5,
                       sd_m_f = .01,
                       sd_m_i = 3.2,
                       ref_map = geomag_map_ref(tag),
                       quiet = FALSE) {
  GeoPressureR:::tag_assert(tag, "setmap")
  assertthat::assert_that(is.logical(compute_known))
  if (!inherits(ref_map, "SpatRaster")) {
    cli::cli_abort(c(
      "{.var ref_map} must be a {.cls SpatRaster} with columns {.var lon}, {.var lat},
      {.var intensity}, and {.var inclinaison}."
    ))
  }
  g <- GeoPressureR:::map_expand(tag$param$tag_set_map$extent, tag$param$tag_set_map$scale)
  # Check ref_map is coheren with dimension of g
  if (!all(dim(ref_map)[c(1, 2)] == c(g$dim[1], g$dim[2]))) {
    cli::cli_abort(c(
      "{.var ref_map} must have dimensions {.val {g$dim[1]}} x {.val {g$dim[2]}} (same as
      {.var tag$param$tag_set_map}).",
      ">" = "Current dimensions are {.val {dim(ref_map)[1]}} x {.val {dim(ref_map)[2]}}."
    ))
  }

  check_sd <- function(x, name) {
    assertthat::assert_that(is.numeric(x))
    if (length(x) == 1) {
      x <- rep(x, times = nrow(tag$stap))
    } else if (length(x) != nrow(tag$stap)) {
      cli::cli_abort(c(
        "x" = "{.var {name}} is of length {.val {length(x)}}.",
        ">" = "{.var {name}} must be length {.val 1} or {.val {nrow(tag$stap)}}."
      ))
    }
    assertthat::assert_that(all(x >= 0))
    x
  }
  sd_e_i <- check_sd(sd_e_i, "sd_e_i")
  sd_e_f <- check_sd(sd_e_f, "sd_e_f")
  sd_m_i <- check_sd(sd_m_i, "sd_m_i")
  sd_m_f <- check_sd(sd_m_f, "sd_m_f")

  mag <- geomag_clean(tag)
  likelihood <- function(v, map, sd_m, sd_e) {
    n <- length(v)
    if (sd_m == 0) {
      mse <- sapply(map, \(x) mean((x - v)^2))
      l <- (1 / (2 * pi * sd_e^2))^(n / 2) * exp(-n / (2 * sd_e^2) * mse)
    } else {
      sum_Y_minus_X <- sapply(map, \(x) sum(x - v))
      sum_Y_minus_X_squared <- sapply(map, \(x) sum((x - v)^2))
      a <- n / sd_e^2 + 1 / sd_m^2
      b <- sum_Y_minus_X / sd_e^2
      ll <- -0.5 * n * log(2 * pi * sd_e^2) - 0.5 * log(2 * pi * sd_m^2)
      ll <- ll + 0.5 * log(2 * pi / a)
      ll <- ll - (0.5 / sd_e^2) * sum_Y_minus_X_squared
      ll <- ll + (b^2 / (2 * a))
      l <- exp(ll - max(ll))
    }
    l <- matrix(l, nrow = dim(map)[1], ncol = dim(map)[2])
    l[tag$map_pressure_mse$mask_water] <- NA
    l
  }

  map_mag <- vector("list", nrow(tag$stap))
  map_mag_f <- vector("list", nrow(tag$stap))
  map_mag_i <- vector("list", nrow(tag$stap))

  for (istap in seq_len(nrow(tag$stap))) {
    mag_istap <- mag[istap == mag$stap_id, ]
    map_mag_f[[istap]] <- likelihood(
      mag_istap$F,
      as.matrix(ref_map$intensity, wide = TRUE),
      sd_m_f[istap], sd_e_f[istap]
    )
    map_mag_i[[istap]] <- likelihood(
      mag_istap$I * 180 / pi,
      as.matrix(ref_map$inclinaison, wide = TRUE),
      sd_m_i[istap], sd_e_i[istap]
    )
    map_mag[[istap]] <- map_mag_f[[istap]] * map_mag_i[[istap]]
  }

  # Override likelihood at known location if requested
  if (!compute_known) {
    for (stap_id in tag$stap$stap_id[!is.na(tag$stap$known_lat)]) {
      zero_matrix <- matrix(0, nrow = nrow(ref_map), ncol = ncol(ref_map))
      cells <- terra::cellFromXY(
        ref_map,
        cbind(tag$stap$known_lon[stap_id], tag$stap$known_lat[stap_id])
      )
      rc <- terra::rowColFromCell(ref_map, cells)
      zero_matrix[cbind(rc[, 1], rc[, 2])] <- 1
      map_mag_f[[stap_id]] <- zero_matrix
      map_mag_i[[stap_id]] <- zero_matrix
      map_mag[[stap_id]] <- zero_matrix
    }
  }

  tag$map_magnetic_intensity <- GeoPressureR::map_create(
    data   = map_mag_f,
    extent = tag$param$tag_set_map$extent,
    scale  = tag$param$tag_set_map$scale,
    stap   = tag$stap,
    id     = tag$param$id,
    type   = "magnetic_intensity"
  )
  tag$map_magnetic_inclination <- GeoPressureR::map_create(
    data   = map_mag_i,
    extent = tag$param$tag_set_map$extent,
    scale  = tag$param$tag_set_map$scale,
    stap   = tag$stap,
    id     = tag$param$id,
    type   = "magnetic_inclination"
  )
  tag$map_magnetic <- GeoPressureR::map_create(
    data   = map_mag,
    extent = tag$param$tag_set_map$extent,
    scale  = tag$param$tag_set_map$scale,
    stap   = tag$stap,
    id     = tag$param$id,
    type   = "magnetic"
  )

  tag$param$geomag_map$sd_e_f <- sd_e_f
  tag$param$geomag_map$sd_e_i <- sd_e_i
  tag$param$geomag_map$sd_m_f <- sd_m_f
  tag$param$geomag_map$sd_m_i <- sd_m_i
  tag
}

#' Compute Reference Magnetic Map from WMM
#'
#' @description
#' Computes raster maps of expected magnetic field intensity and inclination over the region
#' defined in the tag object, using the World Magnetic Model (WMM).
#'
#' @param tag A GeoPressureR tag object.
#' @param quiet Logical. If TRUE, suppresses progress messages.
#' @return A list of two terra raster layers: `intensity` (Gauss) and `inclinaison` (degrees).
#' @noRd
geomag_map_ref <- function(tag, quiet = FALSE) {
  height <- -(288.15 / -0.0065) * (1 - ((stats::median(tag$pressure$value) / 1013.25)^(1 / 5.2561)))
  time <- stats::median(tag$magnetic$date)
  g <- GeoPressureR:::map_expand(tag$param$tag_set_map$extent, tag$param$tag_set_map$scale)
  m <- expand.grid(lat = g$lat, lon = g$lon)

  if (!quiet) {
    pb <- cli::cli_progress_bar("Calculating Magnetic Field map", total = nrow(m))
    map_mag <- mapply(\(loni, lati) {
      tmp <- wmm::GetMagneticFieldWMM(loni, lati, height, time)
      cli::cli_progress_update(id = pb)
      c(tmp$f, tmp$i)
    }, m$lon, m$lat)
    cli::cli_progress_done(id = pb)
  } else {
    map_mag <- mapply(\(loni, lati) {
      tmp <- wmm::GetMagneticFieldWMM(loni, lati, height, time)
      c(tmp$f, tmp$i)
    }, m$lon, m$lat)
  }

  map_f <- terra::rast(
    matrix(map_mag[1, ] / 100000, nrow = g$dim[1], ncol = g$dim[2]),
    extent = tag$param$tag_set_map$extent,
    crs = "epsg:4326"
  )
  names(map_f) <- "intensity"

  map_i <- terra::rast(
    matrix(map_mag[2, ], nrow = g$dim[1], ncol = g$dim[2]),
    extent = tag$param$tag_set_map$extent,
    crs = "epsg:4326"
  )
  names(map_i) <- "inclinaison"

  c(map_f, map_i)
}

#' Clean and Filter Magnetic Data for Likelihood Calculation
#'
#' @description
#' Cleans and filters the magnetic data prior to likelihood map computation:
#' - Assigns stap_id if missing.
#' - Removes points outside valid magnetic field range.
#' - Removes outliers per stationary period for both intensity and inclination.
#'
#' @param tag A GeoPressureR tag object.
#' @return A filtered data.frame of magnetic observations.
#' @noRd
geomag_clean <- function(tag) {
  mag <- tag$magnetic
  if (!("stap_id" %in% names(mag))) {
    mag$stap_id <- GeoPressureR:::find_stap(tag$stap, mag$date)
  }
  mag <- mag[mag$F > .25 & mag$F < .65, ]
  G <- round(mag$stap_id)
  outlier_F <- unsplit(lapply(split(mag$F, G), is_outlier), G)
  outlier_I <- unsplit(lapply(split(mag$I, G), is_outlier), G)
  mag <- mag[!outlier_F & !outlier_I, ]
  mag
}
