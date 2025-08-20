#' Compute Rotation Matrix from Roll, Pitch, and Yaw Angles
#'
#' @description
#' Constructs a 3D rotation matrix given roll, pitch, and yaw angles (in radians).
#' Can take a vector of three angles (`c(roll, pitch, yaw)`) or separate arguments.
#'
#' @param x Either roll angle (radian) or a vector of length 3: c(roll, pitch, yaw).
#' @param pitch Pitch angle (radian), ignored if x is a vector.
#' @param yaw Yaw angle (radian), ignored if x is a vector.
#'
#' @return 3x3 rotation matrix.
#' @seealso [geomag_calib_rotate()]
#' @examples
#' rot(0.1, 0.2, 0.3)
#' rot(c(0.1, 0.2, 0.3))
#' @noRd
rot <- function(x, pitch, yaw) {
  if (length(x) == 3) {
    roll <- x[1]
    pitch <- x[2]
    yaw <- x[3]
  } else {
    roll <- x
  }
  cc <- cos(roll)
  ss <- sin(roll)
  rx <- matrix(c(1, 0, 0, 0, cc, ss, 0, -ss, cc), nrow = 3, ncol = 3, byrow = FALSE)
  cc <- cos(pitch)
  ss <- sin(pitch)
  ry <- matrix(c(cc, 0, -ss, 0, 1, 0, ss, 0, cc), nrow = 3, ncol = 3, byrow = FALSE)
  cc <- cos(yaw)
  ss <- sin(yaw)
  rz <- matrix(c(cc, 0, ss, 0, 1, 0, -ss, 0, cc), nrow = 3, ncol = 3, byrow = FALSE)
  rz %*% ry %*% rx
}

#' Refine Rotation Matrix to Align Major Axis with Diagonal
#' @param rotM 3x3 rotation matrix.
#' @param radius_shape Numeric vector, ellipsoid radii.
#' @return List with elements: `rotM`, `radius_shape`.
#' @noRd
rot_align <- function(rotM, radius_shape) {
  # Align largest component to diagonal
  if (max(abs(rotM)) > 0) {
    rm <- which(abs(rotM) == max(abs(rotM)), arr.ind = TRUE)[1]
    cm <- which(abs(rotM) == max(abs(rotM)), arr.ind = TRUE)[2]
  }
  if (rm != cm) {
    t <- rotM[, cm]
    rotM[, cm] <- rotM[, rm]
    rotM[, rm] <- t
    t <- radius_shape[cm]
    radius_shape[cm] <- radius_shape[rm]
    radius_shape[rm] <- t
  }
  if (rm == 1) k <- c(2, 3)
  if (rm == 2) k <- c(1, 3)
  if (rm == 3) k <- c(1, 2)
  rm <- 0
  cm <- 0
  m <- 0
  r <- c(1:2)
  c <- c(1:2)
  for (i in seq_along(r)) {
    for (j in seq_along(c)) {
      if (abs(rotM[k[i], k[j]]) > m) {
        m <- abs(rotM[k[i], k[j]])
        rm <- k[i]
        cm <- k[j]
      }
    }
  }
  if (rm != cm) {
    t <- rotM[, cm]
    rotM[, cm] <- rotM[, rm]
    rotM[, rm] <- t
    t <- radius_shape[cm]
    radius_shape[cm] <- radius_shape[rm]
    radius_shape[rm] <- t
  }
  if (rotM[1, 1] < 0) rotM[, 1] <- -rotM[, 1]
  if (rotM[2, 2] < 0) rotM[, 2] <- -rotM[, 2]
  if (rotM[3, 3] < 0) rotM[, 3] <- -rotM[, 3]
  list(rotM = rotM, radius_shape = radius_shape)
}

#' Identify Outliers in a Numeric Vector (MAD Method)
#' @param x Numeric vector.
#' @param thresh Threshold (default 3).
#' @return Logical vector: TRUE if outlier.
#' @noRd
is_outlier <- function(x, thresh = 3) {
  med <- stats::median(x, na.rm = TRUE)
  mad_val <- stats::mad(x, na.rm = TRUE)
  abs(x - med) / (mad_val + .Machine$double.eps) > thresh
}

#' Identify Outliers in 3D Data (Any Axis)
#' @param x Numeric matrix/data.frame with 3 columns.
#' @param thresh Numeric, MAD threshold (default 3).
#' @return Logical vector: TRUE if any axis is an outlier for a row.
#' @noRd
is_outlier_3d <- function(x, thresh = 3) {
  stopifnot(ncol(x) == 3)
  med <- apply(x, 2, stats::median, na.rm = TRUE)
  mad_val <- apply(x, 2, stats::mad, na.rm = TRUE)
  scaled <- abs(sweep(x, 2, med)) / (mad_val + .Machine$double.eps)
  apply(scaled, 1, function(row) any(row > thresh))
}

#' Classify Static/Motion State Based on Acceleration
#' @param tag GeoPressureR tag object, must contain `$magnetic` (and optionally `$acceleration`).
#' @param thresh_hard Numeric, absolute threshold for 1g (default 0.1).
#' @param thresh_act Numeric, activity threshold (default 10).
#' @param thresh Numeric, MAD threshold for outlier removal (default 3).
#' @return Tag object with updated `magnetic$is_static` logical column.
#' @noRd
is_static <- function(tag,
                      thresh_hard = .1, # g
                      thresh_act = 10,
                      thresh = 3) {
  if (is.null(tag$magnetic)) {
    cli::cli_abort("No magnetic data found in the tag object.")
  }
  mag <- tag$magnetic
  gn <- sqrt(mag$acceleration_x^2 + mag$acceleration_y^2 + mag$acceleration_z^2)
  is_static <- abs(gn - 1) < thresh_hard
  if (!is.null(tag$acceleration)) {
    acc <- tag$acceleration
    acc$act <- acc$value
    mag <- merge(mag, acc[, c("date", "act")], by = "date", all.x = TRUE)
    is_static[mag$act > thresh_act] <- FALSE
  }
  is_static[is_static] <- !is_outlier(gn[is_static], thresh = thresh)
  mag$is_static <- is_static
  tag$magnetic <- mag
  tag
}
