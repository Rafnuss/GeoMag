utils::globalVariables(c(".data"))

#' Interactive Visualization of Magnetic and Acceleration Data
#'
#' @description
#' Provides interactive 3D and 2D plots for exploring [GeoPressureR tag object
#' ](https://raphaelnussbaumer.com/GeoPressureR/reference/tag_create.html) sensor data, including
#' raw and calibrated magnetic data, fitted calibration ellipsoids, and (projected) acceleration
#' data. Color-coding by stationary period or movement state is supported. Also supports time
#' series and histogram error plots against reference paths.
#'
#' Plot types:
#' \itemize{
#'   \item \strong{"magnetic"}: Raw magnetic data in sensor frame, colored by stationary period.
#'   \item \strong{"calib"}: Calibration points and fitted ellipsoids for selected periods.
#'   \item \strong{"acceleration"}: Raw acceleration, colored by static/moving classification.
#'   \item \strong{"acceleration_p"}: Projected acceleration (NED frame), with reference gravity.
#'   \item \strong{"timeseries"}: Time series of inclination and intensity, with optional reference
#'    path overlay.
#'   \item \strong{"histogram"}: Histogram of errors (sample and mean per stationary period)
#'   against reference path.
#' }
#'
#' @param tag A [GeoPressureR tag object
#' ](https://raphaelnussbaumer.com/GeoPressureR/reference/tag_create.html) containing magnetic (and
#'  optionally calibration) data.
#' @param type Character, plot type. One of "magnetic", "calib", "acceleration", "acceleration_p",
#'  "timeseries", or "histogram".
#' @param stap_id Integer or vector, stationary period(s) to plot calibration fit for
#' (type="calib").
#' @param path Optional, a data frame with columns `start`, `end`, `stap_id`, `lon`, `lat` for
#' plotting reference paths in "timeseries" or "histogram" types.
#'
#' @return A plotly `plot_ly` object (interactive 3D/2D plot) or a ggplotly object for time
#' series/histogram types.
#'
#' @details
#' - Uses the `scico` or `viridisLite` palettes for clear color separation.
#' - For type "calib", if `stap_id` is missing, periods with min/max radii are shown.
#' - For acceleration plots, static/moving state is estimated if not present.
#' - For "timeseries" and "histogram" types, a reference path is required for error analysis.
#'
#' @examples
#' library(GeoPressureR)
#' withr::with_dir(system.file("extdata", package = "GeoMag"), {
#'   tag <- tag_create("14DM", quiet = TRUE)
#'   tag <- tag_label(tag, quiet = TRUE)
#'   tag <- geomag_calib(tag, quiet = TRUE)
#' })
#' plot_mag(tag, type = "acceleration")
#' plot_mag(tag, type = "magnetic")
#' plot_mag(tag, type = "calib")
#' @export
plot_mag <- function(tag, type = "magnetic", stap_id = NULL, path = NULL) {
  # Ensure stap_id exists for all points
  if (!"stap_id" %in% names(tag$magnetic)) {
    tag$magnetic$stap_id <- 1
    stap_id <- 1
    tag$stap <- data.frame(stap_id = 1)
  }
  cols <- get_stap_palette(nrow(tag$stap))

  if (type == "magnetic") {
    p <- plotly::plot_ly() |>
      add_3d_scatter(tag$magnetic, "magnetic_x", "magnetic_y", "magnetic_z", "stap_id", cols)
  } else if (type == "calib") {
    if (!"I" %in% names(tag$magnetic)) {
      cli::cli_abort(c(
        x = "Magnetic data has not yet be calibrated.",
        ">" = "Please run {.fun geomag_calib} first."
      ))
    }

    if (!"stap_id" %in% names(tag$mag_calib)) {
      tag$mag_calib$stap_id <- 1
      stap_id <- 1
    }
    if (is.null(stap_id)) {
      tstap <- table(tag$mag_calib$stap_id)
      poss_stap <- as.numeric(names(tstap[tstap > 5]))
      stap_id <- poss_stap[c(
        which.min(tag$param$geomag_calib$radius_amplitude[poss_stap]),
        which.max(tag$param$geomag_calib$radius_amplitude[poss_stap])
      )]
    }
    stap_id[stap_id < 0] <- max(tag$mag_calib$stap_id) - 1 - stap_id[stap_id < 0]


    p <- plotly::plot_ly() |>
      add_3d_scatter(tag$mag_calib, "magnetic_x", "magnetic_y", "magnetic_z", "stap_id", cols)

    for (i_stap in stap_id) {
      radius <- tag$param$geomag_calib$radius_shape *
        tag$param$geomag_calib$radius_amplitude[i_stap]
      p <- p |>
        add_ellipsoid_mesh(
          radius = radius,
          offset = tag$param$geomag_calib$offset,
          color = cols[i_stap]
        )
    }
  } else if (type %in% c("acceleration", "acceleration_p")) {
    # Combine acceleration and acceleration_p logic
    if (!"is_static" %in% names(tag$magnetic)) {
      tag <- is_static(tag)
    }
    tag$magnetic$is_static_label <- factor(
      ifelse(tag$magnetic$is_static, "Static", "Moving"), c("Moving", "Static")
    )
    if (type == "acceleration") {
      cols_acc <- c("acceleration_x", "acceleration_y", "acceleration_z")
      axis_titles <- list(
        title = "Raw Acceleration (sensor frame)",
        x = "X (forward)", y = "Y (right)", z = "Z (down)"
      )
    } else {
      if (!"acceleration_xp" %in% names(tag$magnetic)) {
        cli::cli_abort(c(
          x = "Projected acceleration data (acceleration_xp) is missing.",
          ">" = "Please run {.fun geomag_calib} first."
        ))
      }
      cols_acc <- c("acceleration_xp", "acceleration_yp", "acceleration_zp")
      axis_titles <- list(
        title = "Projected Acceleration (Horizontal plane of Earth)",
        x = "X_p (forward)", y = "Y_p (right)", z = "Z_p (down)"
      )
    }
    p <- plotly::plot_ly() |>
      add_3d_scatter(
        data = tag$magnetic,
        xcol = cols_acc[1],
        ycol = cols_acc[2],
        zcol = cols_acc[3],
        colorcol = "is_static_label",
        colors = viridisLite::viridis(5)
      ) |>
      add_ellipsoid_mesh(color = "lightblue") |>
      plotly::layout(
        title = axis_titles$title,
        scene = list(
          xaxis = list(title = axis_titles$x),
          yaxis = list(title = axis_titles$y),
          zaxis = list(title = axis_titles$z)
        )
      )

    if (type == "acceleration_p") {
      # Plot the center
      p <- p |> plotly::add_trace(
        x = c(0, 0), y = c(0, 0), z = c(0, 1),
        type = "scatter3d",
        mode = "lines",
        line = list(color = "red", width = 4),
        name = "Reference gravity field"
      )
    }
  } else if (type == "timeseries" || type == "histogram") {
    if (!"I" %in% names(tag$magnetic)) {
      cli::cli_abort(c(
        x = "Magnetic data has not yet be calibrated.",
        ">" = "Please run {.fun geomag_calib} first."
      ))
    }

    mag <- geomag_clean(tag) |>
      dplyr::filter(stap_id == round(stap_id)) |>
      dplyr::mutate(stap_id = factor(stap_id)) |>
      dplyr::mutate(I = .data$I * 180 / pi) |>
      tidyr::pivot_longer(cols = c(.data$I, .data$F), names_to = "variable", values_to = "value") |>
      dplyr::select("date", "stap_id", "variable", "value")

    segments <- mag |>
      dplyr::group_by(.data$stap_id, .data$variable) |>
      dplyr::summarise(
        start = min(date),
        end = max(date),
        mean_val = mean(.data$value, na.rm = TRUE),
        .groups = "drop"
      )

    if (!is.null(path)) {
      time <- as.POSIXct(rowMeans(cbind(path$start, path$end)))
      path[c("F", "I")] <- t(vapply(
        seq_len(nrow(path)),
        \(i) {
          out <- wmm::GetMagneticFieldWMM(path$lon[i], path$lat[i], 0, time[i])
          c(out$f / 100000, out$i)
        },
        numeric(2)
      ))
      path_long <- path |>
        dplyr::select("start", "end", "stap_id", "F", "I") |>
        tidyr::pivot_longer(cols = c(.data$F, .data$I), names_to = "variable", values_to = "val") |>
        dplyr::mutate(stap_id = factor(.data$stap_id))
    }

    if (type == "timeseries") {
      p <- ggplot2::ggplot(mag, ggplot2::aes(x = .data$date, y = .data$value)) +
        ggplot2::geom_point(size = 1.5, alpha = 0.5) +
        ggplot2::geom_segment(
          data = segments,
          ggplot2::aes(
            x = .data$start, xend = .data$end, y = .data$mean_val, yend = .data$mean_val,
            color = .data$stap_id
          ),
          linewidth = 1.2,
          show.legend = FALSE
        )

      if (!is.null(path)) {
        p <- p +
          ggplot2::geom_segment(
            data = path_long,
            ggplot2::aes(x = .data$start, xend = .data$end, y = .data$val, yend = .data$val),
            color = "red",
            linewidth = 1.2,
            show.legend = FALSE
          )
      }

      p <- p +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::facet_wrap(
          ~variable,
          ncol = 1,
          scales = "free_y",
          labeller = ggplot2::labeller(variable = c(
            I = "Inclinaison (\u00B0)",
            F = "Intensity (nT)"
          ))
        ) +
        ggplot2::theme_minimal()

      p <- plotly::ggplotly(p) |> plotly::layout(showlegend = FALSE)
    } else {
      if (is.null(path)) {
        cli::cli_abort("{.arg path} is required for histogram plots.")
      }

      err <- mag |>
        dplyr::left_join(path_long, by = c("stap_id", "variable")) |>
        dplyr::mutate(
          err = .data$value - .data$val
        ) |>
        dplyr::select(-c("date", "start", "end", "value", "val"))

      err_stap <- err |>
        dplyr::group_by(.data$stap_id, .data$variable) |>
        dplyr::summarise(
          err = mean(err, na.rm = TRUE),
          .groups = "drop"
        )

      df <- rbind(
        err |> dplyr::mutate(type = "observation (\\u03c3_e)"),
        err_stap |> dplyr::mutate(type = "stap (\\u03c3_m)")
      )

      sds <- df |>
        dplyr::group_by(.data$type, .data$variable) |>
        dplyr::summarise(sd = stats::sd(err, na.rm = TRUE), .groups = "drop")

      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$err)) +
        ggplot2::geom_histogram(bins = 40) +
        ggplot2::facet_grid(type ~ variable, scales = "free") +
        ggplot2::geom_text(
          data = sds,
          ggplot2::aes(
            x = Inf, y = Inf,
            label = paste0("SD=", round(sd, 4))
          ),
          hjust = 1.1, vjust = 1.5, inherit.aes = FALSE
        ) +
        ggplot2::labs(x = "Error obs - WMM", y = NULL)
    }
  }
  p
}


# ---- Utility functions for plot_mag ----
# Generate color palette for stap_id
get_stap_palette <- function(n) {
  cols <- scico::scico(n, palette = "romaO")
  n2 <- length(cols) %/% 2
  cols <- c(cols[(n2 + 1):length(cols)], cols[1:n2])
  stats::setNames(cols, seq_len(n))
}

# Add ellipsoid mesh overlay to a plotly object
add_ellipsoid_mesh <- function(p, radius = c(1, 1, 1), offset = c(0, 0, 0), color = "lightblue") {
  # Mesh grid for ellipsoid overlays
  phi <- seq(0, 360, length.out = 100) / 180 * pi
  theta <- seq(-90, 90, length.out = 100) / 180 * pi
  x <- c(sin(phi) %*% t(cos(theta)))
  y <- c(sin(phi) %*% t(sin(theta)))
  z <- c(cos(phi) %*% t(rep(1, length(theta))))

  suppressWarnings(
    p |> plotly::add_mesh(
      x = x * radius[1] + offset[1],
      y = y * radius[2] + offset[2],
      z = z * radius[3] + offset[3],
      alphahull = 0,
      opacity = 0.2,
      showscale = FALSE,
      intensity = rep(1, length(x)),
      colorscale = list(c(0, 1), c(color, color)),
      color = I(color),
      showlegend = FALSE
    )
  )
}
# Helper for 3D scatter plot
add_3d_scatter <- function(p, data, xcol, ycol, zcol, colorcol, colors, ...) {
  data <- data[stats::complete.cases(data[, c(xcol, ycol, zcol, colorcol)]), ]
  suppressWarnings(
    p |> plotly::add_markers(
      data = data,
      x = data[[xcol]],
      y = data[[ycol]],
      z = data[[zcol]],
      color = data[[colorcol]],
      colors = colors,
      text = data$stap_id,
      ...
    )
  )
}
