#' run_xeno
#' 
#' Run Xeno Mouse Fraction Analysis
#' 
#' @return a \code{data.frame} with the estimated mouse content
#' @param idat_path Path to the idat file for which the mouse fraction is to be computed
#' @param n_cores The number of cores to be used
#' @export
#' @author Nima Esmaeelpour
rnb.run_xeno <- function(idat_path, n_cores = 1) {
  rnb.require('sesame')
  rnb.require('BiocParallel')
  rnb.require('dplyr')
    if (!is.character(idat_path) || length(idat_path) != 1) {
    stop("`idat_path` must be a single character string.", call. = FALSE)
  }
  if (!dir.exists(idat_path)) {
    stop("Directory does not exist: ", idat_path, call. = FALSE)
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) {
    stop("`n_cores` must be a positive integer.", call. = FALSE)
  }
  
  # ---- LOAD DATA ----
  
  idat_files <- list.files(idat_path, "\\.idat$")
  if (length(idat_files) == 0) {
    stop("No .idat files found in: ", idat_path, call. = FALSE)
  }

  # Read .idat files
  logger.info(sprintf("Started Reading %s directory...", idat_path))
  sdfs <- openSesame(
    idat_path,
    func = NULL,
    BPPARAM = MulticoreParam(n_cores)
  )
  logger.info(sprintf("Finished Reading %s directory...", idat_path))

  # Convert to list if only one sample (i.e. two .idat files)
  if (!is.list(sdfs) || is(sdfs, "SigDF")) {
    sdfs <- list(sample = sdfs)
  }

  # Add xeno content
  results <- vector("list", length(sdfs))
  i <- 1
  for (sample_name in names(sdfs)) {
    results[[i]] <- add_xeno(sdfs[[sample_name]], sample_name)
    i <- i + 1
  }

  # Return final df
  dplyr::bind_rows(results)
}

# ---- HELPER FUNCTIONS ----

standard_error <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  stats::sd(x) / sqrt(length(x))
}

filter_type_I_probes <- function(df, manifest) {
  dplyr::filter(df, Probe_ID %in% manifest$IlmnID)
}

calculate_wrong_color_fraction <- function(df, manifest) {
  channel_map <- dplyr::select(manifest, IlmnID, Color_Channel)

  df |>
    dplyr::left_join(channel_map, by = c("Probe_ID" = "IlmnID")) |>
    dplyr::mutate(
      total_signal = MR + UR + MG + UG,
      wrong_color_signal = dplyr::if_else(
        Color_Channel == "Grn",
        MR + UR,
        MG + UG
      ),
      wrong_color_fraction = wrong_color_signal / total_signal
    )
}

# ---- CORE LOGIC ----

MIN_TOTAL_SIGNAL <- 1000
LOW_SIGNAL_CUTOFF <- 0.4

#' Compute xenograft (interspecies) contamination metrics for one sample
#'
#' @param df          Raw probe-level data frame (columns: Probe_ID, MR, UR, MG, UG)
#' @param sample_name Character scalar identifying the sample
#' @return Single-row data.frame with contamination metrics
#' @noRd
add_xeno <- function(df, sample_name) {
  platform <- resolve_platform(df, sample_name)
  manifest <- platform$manifest
  interspecies_probes <- platform$interspecies_probes

  df <- df |>
    filter_type_I_probes(manifest) |>
    calculate_wrong_color_fraction(manifest) |>
    dplyr::filter(MR + UR + MG + UG > MIN_TOTAL_SIGNAL)

  background <- stats::median(df$wrong_color_fraction, na.rm = TRUE)
  background_se <- standard_error(df$wrong_color_fraction, na.rm = TRUE)
  background_mad <- stats::mad(df$wrong_color_fraction, na.rm = TRUE)

  xeno_signal_stats <- df |>
    dplyr::filter(Probe_ID %in% interspecies_probes) |>
    dplyr::summarize(
      n_interspecies_probes = dplyr::n(),
      median_signal = stats::median(wrong_color_fraction, na.rm = TRUE),
      se_signal = standard_error(wrong_color_fraction, na.rm = TRUE),
      mad_signal = stats::mad(wrong_color_fraction, na.rm = TRUE)
    )

  mouse_fraction <- compute_mouse_fraction(
    xeno_signal_stats$median_signal,
    background
  )

  data.frame(
    sample = sample_name,
    probes_used = xeno_signal_stats$n_interspecies_probes,
    mouse_signal = xeno_signal_stats$median_signal,
    sem_mouse_signal = xeno_signal_stats$se_signal,
    mad_mouse_signal = xeno_signal_stats$mad_signal,
    background_signal = background,
    sem_background = background_se,
    mad_background = background_mad,
    mouse_fraction = mouse_fraction,
    stringsAsFactors = FALSE
  )
}

#' Adjust raw mouse signal by background, clamped to `[0, 1]`
#' Returns 1 (worst case) when no interspecies probes are available.
#' @noRd
compute_mouse_fraction <- function(mouse_med, background) {
  dplyr::case_when(
    is.na(mouse_med) ~ 1,
    mouse_med < LOW_SIGNAL_CUTOFF ~ pmax(0, mouse_med - background),
    TRUE ~ mouse_med
  )
}

#' Resolve platform and return manifest and interspecies probes
#'
#' @param df Raw probe-level data frame
#' @param sample_name Character scalar identifying the sample
#' @return List with manifest and interspecies_probes
#' @noRd
resolve_platform <- function(df, sample_name) {
  load(system.file("data/xenomix.RData", package = "RnBeads"))
  n <- nrow(df)

  if (n == 937690) {
    list(
      manifest = epic_v2_manifest,
      interspecies_probes = interspecies_probes_v2
    )
  } else {
    list(
      manifest = epic_v1_manifest,
      interspecies_probes = interspecies_probes_v1
    )
  }
}
          