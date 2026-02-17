#' run_xeno
#' 
#' Run Xeno Mouse Fraction Analysis
#' 
#' @return a \code{data.frame} with the estimated mouse content
#' @param idat_path Path to the idat file for which the mouse fraction is to be computed
#' @param n_cores The number of cores to be used
#' @export
#' @author Nima Esmaeelpour
#' @importFrom dplyr bind_rows
#' @importFrom sesame openSesame
#' @importFrom methods is
rnb.run_xeno <- function(idat_path, n_cores = 1) {
  rnb.require('sesame')
  rnb.require('BiocParallel')
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
  bind_rows(results)
}

#' @noRd
#' @importFrom dplyr filter
#' @importFrom dplyr summarize
#' @importFrom dplyr pull
#' @importFrom stats median
add_xeno <- function(df, sample_name) {
  platform <- resolve_platform(df, sample_name)
  manifest <- platform$manifest
  interspecies_probes <- platform$interspecies_probes

  df <- filter_type_I_probes(df, manifest)
  df <- calculate_wrong_color_fraction(df, manifest)
  df <- df %>% filter(rowSums(.[2:5]) > 1000)
  bg <- median(df$wrong_color_fraction, na.rm = TRUE)

  mouse <- df |>
    filter(Probe_ID %in% interspecies_probes) |>
    summarize(med = median(wrong_color_fraction)) |>
    pull(med)

  frac <- ifelse(is.na(mouse), 1,
    ifelse(mouse < 0.33, pmax(0, mouse - bg), mouse)
  )

  data.frame(
    sample = sample_name,
    background_signal = bg,
    mouse_signal = mouse,
    mouse_fraction = frac,
    stringsAsFactors = FALSE
  )
}

#' @noRd
#' @importFrom dplyr filter
filter_type_I_probes <- function(df, manifest) {
  df <- df |> filter(Probe_ID %in% manifest$IlmnID)
}

#' @noRd
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr select
calculate_wrong_color_fraction <- function(df, manifest) {
  df |>
    left_join(
      manifest |> select(IlmnID, Color_Channel),
      by = c("Probe_ID" = "IlmnID")
    ) |>
    mutate(
      wrong_color_fraction =
        ifelse(Color_Channel == "Grn",
          (MR + UR) / (MR + UR + MG + UG),
          (MG + UG) / (MR + UR + MG + UG)
        )
    )
}

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

          