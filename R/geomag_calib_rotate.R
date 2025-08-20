#' Apply Rotation Matrix to Sensor Data for Frame Alignment
#'
#' @description
#' Rotates each row of a 3-column matrix of sensor data from the sensor's frame to the
#' world (NED) frame, using per-row roll, pitch, and yaw angles. This is typically used for tilt
#' compensation of acceleration or magnetic data.
#'
#' @param a Numeric matrix (N x 3), each row is a 3-vector in the sensor frame.
#' @param roll Numeric vector of roll angles (radians), length 1 or nrow(a).
#' @param pitch Numeric vector of pitch angles (radians), length 1 or nrow(a).
#' @param yaw Numeric vector of yaw angles (radians), length 1 or nrow(a). Default: 0.
#'
#' @return A numeric matrix (N x 3), where each row is the rotated vector in the world/NED frame.
#'
#' @examples
#' # Example: Rotate a single vector by roll = 0.1 rad, pitch = 0.2 rad, yaw = 0
#' v <- matrix(c(1, 0, 0), ncol = 3)
#' geomag_calib_rotate(v, roll = 0.1, pitch = 0.2)
#'
#' @noRd
geomag_calib_rotate <- function(a, roll, pitch, yaw = 0) {
  # Input validation
  if (is.null(roll) || length(roll) == 0) {
    cli::cli_abort("{.arg roll} cannot be {.val NULL} or empty")
  }
  if (is.null(pitch) || length(pitch) == 0) {
    cli::cli_abort("{.arg pitch} cannot be {.val NULL} or empty")
  }
  if (is.null(yaw) || length(yaw) == 0) {
    cli::cli_abort("{.arg yaw} cannot be {.val NULL} or empty")
  }
  assertthat::assert_that(assertthat::are_equal(dim(a)[2], 3))

  # Replicate angles if single value provided
  if (length(roll) == 1) roll <- rep(roll, nrow(a))
  if (length(pitch) == 1) pitch <- rep(pitch, nrow(a))
  if (length(yaw) == 1) yaw <- rep(yaw, nrow(a))

  assertthat::assert_that(assertthat::are_equal(nrow(a), length(roll)))
  assertthat::assert_that(assertthat::are_equal(nrow(a), length(pitch)))
  assertthat::assert_that(assertthat::are_equal(nrow(a), length(yaw)))

  ar <- matrix(NA, nrow(a), ncol(a))
  for (i in seq_len(nrow(a))) {
    ar[i, ] <- rot(roll[i], pitch[i], yaw[i]) %*% a[i, ]
  }
  return(ar)
}
