#' Calibrate Magnetic and Acceleration Data for 3-Axis Sensors
#'
#' @description
#' Performs tilt compensation and magnetic calibration for 3-axis sensor data, with support for
#' in situ and in vitro calibration routines, outlier removal, and computation of orientation and
#' field parameters. This function is designed for use with `GeoPressureR` tag objects and can
#' utilize calibration datasets if available.
#'
#' ## Workflow
#'
#' 1. **Determine static/movement states:**
#'    Classify each data point as static or moving by evaluating acceleration signals.
#' 2. **Magnetic Data Calibration:**
#'    - Select calibration data source (raw or from calibration dataset).
#'    - Optionally remove extreme or outlier values.
#'    - Fit and apply a calibration model (sphere/ellipse or their stap variants).
#' 3. **Tilt Compensation:**
#'    - Compute pitch and roll from acceleration.
#'    - Project gravity and calibrated magnetic data into the horizontal plane of the Earth frame.
#' 4. **Orientation and Field Parameters:**
#'    - Calculate heading, field intensity, and inclination.
#'    - Store calibration metadata and processed data in the tag object.
#'
#' @param tag A `GeoPressureR` tag object containing magnetic and acceleration data.
#' @param calib_data Logical, character, or `NULL`.
#'        If `TRUE`, uses calibration data from `magCalib/` subfolder.
#'        If `FALSE`, calibrates using in situ data.
#'        If a character path, uses calibration data from the specified directory.
#'        If `NULL`, auto-detects calibration folder.
#' @param calib_method Character. Calibration method, one of `"sphere"`, `"ellipse"`,
#'        `"near-sphere"`, `"sphere_stap"`, or `"ellipse_stap"`. If `NULL`, chosen automatically.
#' @param rm_outlier Logical. If `TRUE`, removes outliers from calibration data (recommended).
#' @param quiet Logical. If `TRUE`, suppresses progress messages.
#'
#' @return Modified `GeoPressureR` tag object. The `$magnetic` data frame contains:
#'   - `date`: Timestamp (POSIXct or numeric)
#'   - `acceleration_x`, `acceleration_y`, `acceleration_z`: Raw acceleration data
#'   - `magnetic_x`, `magnetic_y`, `magnetic_z`: Raw magnetic data
#'   - `is_static`: Scaled MAD of acceleration (0 = static, >1 = movement)
#'   - `pitch`, `roll`: Orientation angles (radian)
#'   - `acceleration_xp`, `acceleration_yp`, `acceleration_zp`: Gravity projected in NED frame
#'   - `is_outlier`: Logical, marks outliers in magnetic data
#'   - `magnetic_xc`, `magnetic_yc`, `magnetic_zc`: Calibrated magnetic data
#'   - `magnetic_xcp`, `magnetic_ycp`, `magnetic_zcp`: Calibrated magnetic data projected in NED
#'   frame
#'   - `H`: Heading (radian, North=0)
#'   - `F`: Magnetic field intensity (Gauss)
#'   - `I`: Inclination (radian)
#' Also returns (invisibly) the calibration dataset used (`tag$mag_calib`) and calibration
#' parameters (`tag$param$geomag_calib`).
#'
#'
#' @section Details:
#' This function is part of the `GeoMag` package and is intended for use with animal-attached tags
#' that record magnetic and acceleration data. It supports several calibration workflows and robust
#' outlier detection.
#'
#' @examples
#' library(GeoPressureR)
#' withr::with_dir(system.file("extdata", package = "GeoMag"), {
#'   tag <- tag_create("14DM", quiet = TRUE)
#'   tag <- tag_label(tag, quiet = TRUE)
#'   tag <- geomag_calib(tag, quiet = TRUE)
#'   tag$param$geomag_calib
#'   head(tag$magnetic)
#' })
#' @export
geomag_calib <- function(tag,
                         calib_data = NULL,
                         calib_method = NULL,
                         rm_outlier = TRUE,
                         quiet = FALSE) {
  GeoPressureR::tag_assert(tag, "magnetic")

  # Step 1: Mark static/movement states using acceleration
  tag <- is_static(tag)
  mag <- tag$magnetic

  # Step 2: Determine calibration data source (in situ or calibration folder)
  if (is.null(calib_data)) {
    directory <- glue::glue("{tag$param$tag_create$directory}/magCalib")
    if (dir.exists(directory)) {
      calib_data <- TRUE
    } else {
      calib_data <- FALSE
    }
  } else if (is.character(calib_data) && dir.exists(calib_data)) {
    directory <- calib_data
    calib_data <- TRUE
  } else if (identical(calib_data, TRUE)) {
    directory <- glue::glue("{tag$param$tag_create$directory}/magCalib")
  }

  # Load calibration data if requested/available
  if (identical(calib_data, TRUE)) {
    tag_calib <- GeoPressureR:::tag_create_soi("",
      directory = directory, quiet = quiet
    )
    assertthat::assert_that(
      !is.null(tag_calib$magnetic),
      msg = "No magnetic data found in the calibration data."
    )
    tag_calib <- is_static(tag_calib)
    mag_calib <- tag_calib$magnetic
  } else if (is.data.frame(calib_data)) {
    mag_calib <- calib_data
    cli::cli_alert_info("Using provided data.frame {.arg calib_data}")
  } else {
    # Use in situ data
    mag_calib <- mag
    cli::cli_alert_info("Using raw magnetic data for calibrationd data")
  }

  # Step 3: Remove outliers/extreme values from calibration data
  if (rm_outlier) {
    mag_calib <- geomag_calib_rm(mag_calib, tag)
  }

  # Step 4: Select calibration method if not specified
  if (is.null(calib_method)) {
    if ("stap_id" %in% names(mag_calib)) {
      calib_method <- "ellipse_stap"
    } else {
      calib_method <- "ellipse"
    }
  }

  # Step 5: Calibrate magnetic data using selected method
  mag <- geomag_calib_fit(
    mag = mag,
    mag_calib = mag_calib,
    method = calib_method,
    stap = tag$stap
  )

  # Step 6: Tilt compensation and orientation calculation
  # Compute acceleration magnitude
  gn <- sqrt(mag$acceleration_x^2 + mag$acceleration_y^2 + mag$acceleration_z^2)

  # Pitch: rotation around Y' axis; Roll: rotation around X axis
  mag$pitch <- asin(-mag$acceleration_x / gn)
  mag$roll <- atan2(mag$acceleration_y / gn, mag$acceleration_z / gn)

  # Project gravity into horizontal plane of the Earth frame
  gr <- geomag_calib_rotate(
    matrix(
      c(mag$acceleration_x, mag$acceleration_y, mag$acceleration_z),
      nrow(mag), 3
    ),
    roll = mag$roll, pitch = mag$pitch
  )
  mag$acceleration_xp <- gr[, 1]
  mag$acceleration_yp <- gr[, 2]
  mag$acceleration_zp <- gr[, 3]

  # Project calibrated magnetic data into the horizontal plane of the Earth frame
  mr <- geomag_calib_rotate(
    matrix(
      c(mag$magnetic_xc, mag$magnetic_yc, mag$magnetic_zc),
      nrow(mag), 3
    ),
    roll = mag$roll, pitch = mag$pitch
  )
  mag$magnetic_xcp <- mr[, 1]
  mag$magnetic_ycp <- mr[, 2]
  mag$magnetic_zcp <- mr[, 3]

  # Step 7: Compute heading (H), intensity (F), and inclination (I)
  mag$H <- atan2(-mag$magnetic_ycp, mag$magnetic_xcp)
  mag$F <- sqrt(rowSums(mr^2))
  mag$I <- -asin(mag$magnetic_zcp / mag$F)

  # Store calibration metadata and used calibration data
  tag$param$geomag_calib <- attr(mag, "geomag_calib")
  tag$param$geomag_calib$calib_data <- ifelse(identical(calib_data, TRUE), directory, FALSE)
  tag$param$geomag_calib$rm_outlier <- rm_outlier
  attr(mag, "geomag_calib") <- NULL

  # Save calibration data and processed magnetic data to tag
  tag$mag_calib <- mag_calib
  tag$magnetic <- mag

  return(tag)
}


#' @noRd
geomag_calib_rm <- function(mag_calib, tag) {

  # Remove values with excessive magnetic field intensity
  mn <- sqrt(mag_calib$magnetic_x^2 + mag_calib$magnetic_y^2 +
    mag_calib$magnetic_z^2)
  mag_calib <- mag_calib[mn < 1, ]
  if (nrow(mag_calib) == 0) {
    cli::cli_abort(c(
      "x" = "No calibration data left after removing extreme values.",
      ">" = "Check the magnetic sensor unit for issues."
    ))
    return(mag_calib)
  }

  # Initial offset estimation using spherical model
  mag_tmp <- geomag_calib_fit(
    mag = tag$magnetic,
    mag_calib = mag_calib,
    method = "sphere",
    stap = tag$stap
  )
  offset <- attr(mag_tmp, "geomag_calib")$offset

  # Remove values with physically implausible field strengths
  mn <- sqrt((mag_calib$magnetic_x - offset[1])^2 +
    (mag_calib$magnetic_y - offset[2])^2 +
    (mag_calib$magnetic_z - offset[3])^2)
  mag_calib <- mag_calib[mn > 0.25, ]
  mag_calib <- mag_calib[mn < 0.65, ]

  # Outlier detection within stap groups if available
  mn <- sqrt((mag_calib$magnetic_x - offset[1])^2 +
    (mag_calib$magnetic_y - offset[2])^2 +
    (mag_calib$magnetic_z - offset[3])^2)
  if ("stap_id" %in% names(mag_calib)) {
    G <- round(mag_calib$stap_id)
    split_df <- split(mn, G)
    result <- lapply(split_df, is_outlier)
    mag_calib <- mag_calib[!unsplit(result, G), ]
  } else {
    mag_calib <- mag_calib[!is_outlier(mn), ]
  }

  # Step 4: Ensure sufficient calibration data
  if (nrow(mag_calib) < 10) {
    cli::cli_abort(c(
      "x" = "Removal of outliar lead to remove all data."
    ))
  }

  mag_calib
}
