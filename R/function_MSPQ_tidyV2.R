#' Tidying and preparing the MultispeQ raw dataset.
#'
#' @param df data.frame. The raw data set obtained from the PhotosynQ network.
#' @param genotype string with the name of the genotype column on df.
#' @param time.diff logical. Were the measurements repeated twice in a day? Morning & Afternoon. TRUE if so, FALSE otherwise
#' @param data_name string with the name of the experiment.
#' @param plotIm logical. To plot the imputed variables. FALSE by default.
#'
#' @return A list with 8 objects if plotIm = FALSE. 9 elements otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#' # MSPQ_csv_data: any csv dataset obtained from the PhotosynQ Network.
#'
#' tidy_data <- MSPQ_tidy(
#'   df = MSPQ_csv_data,
#'   genotype = "Genotype",
#'   time.diff = FALSE,
#'   data_name = "My_Project",
#'   plotIm = TRUE
#' )
#' }
MSPQ_tidy <- function(df, genotype, time.diff, data_name = NULL, plotIm = FALSE) {
  

  if ("Leaf.Temperature.Differenial" %in% names(df)) df %<>% dplyr::select(-Leaf.Temperature.Differenial)

  namescol <- c(
    "ID", "Series", "Repeat", "Ambient.Humidity",
    "Ambient.Temperature", "ECSt.mAU", "gH.", "Leaf.Angle",
    "Leaf.Temperature.Differential", "LEF", "Light.Intensity..PAR.", "NPQt",
    "Phi2", "PhiNO", "PhiNPQ", "PS1.Active.Centers",
    "PS1.Open.Centers", "PS1.Over.Reduced.Centers", "PS1.Oxidized.Centers", "Relative.Chlorophyll",
    "Thickness", "vH.", "time", "note",
    "absorbance_420", "absorbance_530", "absorbance_605", "status",
    "absorbance_650", "absorbance_730", "absorbance_850", "absorbance_880",
    "absorbance_940", "air_temp_kinetics", "Ambient.Pressure", "B",
    "data_raw_PAM", "ECS_averaged_trace", "ECS_tau", "fitinput",
    "FmPrime", "FoPrime", "Fs", "FvP_over_FmP",
    "G", "humidity2_K", "humidity2_K_T", "humidity_K",
    "humidity_K_T", "kP700", "labels", "Leaf.Temperature",
    "LEAF_temp", "LEAF_temp_T", "LEFd_trace", "lights_length",
    "n_sets", "outdata", "P700_DIRK_ampl", "P700_DIRK_averaged_trace",
    "P700_fitinput", "P700_outdata", "PSI_data_absorbance", "pump",
    "qL", "R", "SPAD_420", "SPAD_420_intensity",
    "SPAD_530_intensity", "SPAD_605_intensity", "SPAD_650", "SPAD_650_intensity",
    "SPAD_730_intensity", "SPAD_850_intensity", "SPAD_880_intensity", "t",
    "test_data_raw_PAM", "thick2", "thick3", "thick4",
    "Time.of.Day", "tP700", "tt1", "tt2",
    "tt3", "tt4", "v_initial_P700", "value1",
    "value2", "value3", "YII", "User",
    "Device.ID", "Latitude", "Longitude", "Issues"
  )


  SoV_names <- names(df)[!tolower(names(df)) %in% tolower(namescol)]

  if ("leaf_thickness" %in% SoV_names & "SPAD" %in% SoV_names) {
    names(df)[names(df) == "leaf_thickness"] <- "Thickness"
    names(df)[names(df) == "SPAD"] <- "Relative.Chlorophyll"
    SoV_names <- SoV_names[-which(SoV_names %in% c("leaf_thickness", "SPAD"))] # SoV
  }


  ind <- grep(pattern = "^col|^row|^fil", x = tolower(SoV_names))
  if (length(ind) != 0) {
    SoV_names <- SoV_names[-ind]
  }
  rm(ind)

  #-----delete Plot.1 column (repeated)-----
  r_full <- nrow(df)
  c_full <- ncol(df)
  delete <- which(names(df) == "Plot.1")
  if (!installr::is.empty(delete)) {
    df <- df[, -delete]
  }
  rm(delete)
  #-----Replacing null for NA-----
  cat("Replacing null for NA\n")

  for (i in 1:ncol(df)) {
    if (any(df[, i] == "null", na.rm = T)) {
      df[, i] <- as.numeric(gsub(x = df[, i], pattern = "null", replacement = NA))
    }
  }
  rm(i)

  #-----------------Calculating Phi index--------------------------------
  cat("Calculating Phi Index\n")

  df %<>%
    dplyr::mutate(phi_index = Phi2 / (PhiNPQ + PhiNO))

  #-----Formating dates and creating time (morning/afternoon) column-----
  if (time.diff) {
    cat("Formatting dates and creating the time (morning/afternoon) variable\n")
  }

  df$time %<>% as.character

  if (grepl(x = df$time[1], pattern = ":.*AM|PM&AM|PM")) {
    x <- list()

    for (i in 1:length(df$time)) {
      x[[i]] <- unlist(strsplit(df$time[i], " "))[-2]
      x[[i]] <- paste0(x[[i]][1], " ", x[[i]][2])
    }

    x <- unlist(x)
    df$time <- x
    rm(x, i)

    names(df)[names(df) == "time"] <- "date"

    df <- tidyr::separate(df, date, into = c("date", "time"), sep = " ")

    df$date <- as.Date(df$date, format = "%m/%d/%Y")
    df$time <- as.factor(df$time)
  } else {
    stop("There is not AM/PM indicator in column time and/or hour is missing, check out first\n")
  }

  if (!time.diff) {
    df$time <- factor("")
  }

  #-----Discarding rows with issues-----

  cat("Discarding rows with issues\n")

  issues_index <- which(!is.na(df$Issues))
  r_issue <- length(issues_index)
  if (r_issue == 0) {
    issues_index <- c()
  }

  #-----Separating factors and character columns out from df (NOT to be included into final df)-----

  factor_index <- c()

  for (i in 1:ncol(df)) {
    if (is.factor(df[, i])) {
      factor_index[i] <- i
    } else if (is.character(df[, i])) {
      factor_index[i] <- i
    } else if (length(which(is.na(df[, i]))) == nrow(df) || length(unique(df[, i])) == 1 || length(unique(df[, i])) == 2) {
      factor_index[i] <- i
    } else if (length(unique(df[, i])) == 3 || length(which(is.na(df[, i]))) >= round(nrow(df) * 0.5, digits = 0)) {
      factor_index[i] <- i
    } else if (lubridate::is.Date(df[, i])) {
      factor_index[i] <- i
    }
  }
  rm(i)
  factor_index <- factor_index[!is.na(factor_index)]
  factor_index <- c(factor_index, which(names(df) == "ID"))

  #-----Applying control structures-----

  cat("Applying control structures\n")

  SoV_index <- which(names(df) %in% c("date", "time", SoV_names))
  factor_index <- factor_index[!factor_index %in% SoV_index]
  factor_index <- factor_index[!factor_index %in% which(names(df) == "Device.ID")]

  num_df <- df[, -factor_index]
  factor_df <- df[, factor_index]
  c_factor <- length(factor_index)
  rm(factor_index)

  #----Imputation from num_df ---------------

  cat("Data Imputation\n")

  nas <- function(var, val) {
    c(which(num_df[, var] > val), which(num_df[, var] < (val * -1)))
  }

  if ("PS1.Active.Centers" %in% namescol) {
    delete <- nas("PS1.Active.Centers", 50)
    if (!installr::is.empty(delete)) {
      num_df[delete, "PS1.Active.Centers"] <- NA_real_
    }
  }

  if ("PS1.Open.Centers" %in% namescol) {
    delete <- nas("PS1.Open.Centers", 50)
    if (!installr::is.empty(delete)) {
      num_df[delete, "PS1.Open.Centers"] <- NA_real_
    }
  }

  if ("PS1.Over.Reduced.Centers" %in% namescol) {
    delete <- nas("PS1.Over.Reduced.Centers", 50)
    if (!installr::is.empty(delete)) {
      num_df[delete, "PS1.Over.Reduced.Centers"] <- NA_real_
    }
  }

  if ("PS1.Oxidized.Centers" %in% namescol) {
    delete <- nas("PS1.Oxidized.Centers", 50)
    if (!installr::is.empty(delete)) {
      num_df[delete, "PS1.Oxidized.Centers"] <- NA_real_
    }
  }

  if ("NPQt" %in% namescol) {
    delete <- nas("NPQt", 15)
    if (!installr::is.empty(delete)) {
      num_df[delete, "NPQt"] <- NA_real_
    }
  }


  if ("vH." %in% namescol) {
    delete <- nas("vH.", 10)
    if (!installr::is.empty(delete)) {
      num_df[delete, "vH."] <- NA_real_
    }
  }

  if ("ECSt.mAU" %in% namescol) {
    delete <- nas("ECSt.mAU", 10)
    if (!installr::is.empty(delete)) {
      num_df[delete, "ECSt.mAU"] <- NA_real_
    }
  }

  if ("gH." %in% namescol) {
    delete <- nas("gH.", 1100)
    if (!installr::is.empty(delete)) {
      num_df[delete, "gH."] <- NA_real_
    }
  }

  if ("ECS_tau" %in% namescol) {
    delete <- nas("ECS_tau", 100)
    if (!installr::is.empty(delete)) {
      num_df[delete, "ECS_tau"] <- NA_real_
    }
  }
  rm(delete)

  d <- num_df

  d <- d[, -which(names(d) %in% c("Latitude", "Longitude"))]

  d
  miss <- VIM::aggr(d,
    col = c("navyblue", "yellow"), plot = FALSE,
    numbers = TRUE, sortVars = TRUE,
    labels = names(d), cex.axis = .7,
    gap = 3, ylab = c("Missing data", "Pattern")
  )$missings

  miss %<>%
    dplyr::filter(Count != 0) ####### OUTPUT

  if (nrow(miss) > 0) {
    variables <- miss$Variable

    to_sort <- names(num_df)


    to_impute <- num_df[, c("date", variables)]

    new <- list()

    j <- ncol(to_impute) - 1
    pbar <- txtProgressBar(
      min = 0,
      max = j,
      style = 3,
      width = 50,
      char = "="
    )

    for (j in 2:length(to_impute)) {
      var <- to_impute[, c(1, j)] #### c(1,j)


      dates <- lapply(unique(d$date), function(x) {
        dplyr::filter(var, date == x)
      })

      imputes <- list()

      for (i in 1:length(dates)) {
        if (any(is.na(dates[[i]][2]))) {
          imputes[[i]] <- suppressWarnings(mice::mice(dates[[i]], m = 5, maxit = 50, method = "sample", seed = 500, printFlag = FALSE))
        }
      }

      rm(i)

      for (i in 1:length(imputes)) {
        if (!is.null(imputes[[i]])) {
          imputes[[i]] <- mice::complete(imputes[[i]], 5, include = FALSE)
        }
      }

      if (length(imputes) < length(dates)) {
        dif <- length(dates) - length(imputes)
        for (k in 1:dif) {
          imputes[length(imputes) + k] <- list(NULL)
        }
        rm(k)
      }
      rm(i)

      for (i in 1:length(dates)) {
        if (is.null(imputes[[i]])) {
          imputes[[i]] <- dates[[i]]
        }
      }

      dates <- do.call(rbind, imputes)

      new[[j]] <- data.frame(dates[, 2])
      colnames(new[[j]]) <- names(var)[2]
      rm(dates, var, imputes)
      setTxtProgressBar(pbar, j)
    }
    close(pbar)
    cat("\n")

    rm(i, j)

    new[[1]] <- data.frame(date = to_impute$date)

    to_impute <- to_impute[, -1]

    miss <- collapsibleTree::collapsibleTree(miss,
      hierarchy = names(miss),
      collapsed = FALSE, linkLength = 120,
      fontSize = 15
    )

    if (plotIm) {
      plots <- list()

      for (i in 1:length(to_impute)) {
        nas <- which(is.na(to_impute[, i]))
        plots[[i]] <- ggplot2::ggplot(
          cbind(new[[1]], new[[i + 1]]) %>% dplyr::mutate(value = "Measured"),
          ggplot2::aes_string(x = "date", group = "date", y = names(to_impute)[i])
        ) +
          ggplot2::geom_boxplot(ggplot2::aes(fill = value)) +
          ggplot2::scale_fill_manual(values = "grey50") +
          ggplot2::geom_point(
            data = cbind(new[[1]], new[[i + 1]])[nas, ] %>% dplyr::mutate(value = "Imputed"),
            ggplot2::aes_string(x = "date", group = "date", y = names(to_impute)[i], color = "value")
          ) +
          ggplot2::ylim(c(median(new[[i + 1]][, 1]) - (median(new[[i + 1]][, 1]) * 4), median(new[[i + 1]][, 1]) + (median(new[[i + 1]][, 1]) * 4))) +
          ggplot2::labs(title = paste0("Imputed data for ", names(to_impute)[i])) +
          ggplot2::theme_bw()
        names(plots)[i] <- names(to_impute)[i]
      }
      rm(i)
    }

    to_impute <- do.call(cbind, new)

    # rm(new)

    num_df <- suppressWarnings(dplyr::select_(num_df, .dots = names(num_df)[!to_sort %in% variables]))
    num_df <- data.frame(num_df, to_impute[, 2:length(to_impute)])
    num_df <- suppressWarnings(dplyr::select_(num_df, .dots = to_sort))

    rm(d, to_sort, variables)
  } else {
    miss <- "There was not any missing data in raw file."
    if (plotIm) {
      plots <- miss
    }
  }

  #----finding absorbances with NA's (NOT to be included into final df)----

  if (any(grepl(x = names(num_df), pattern = "absorbance_"))) {
    absorbance_index <- which(grepl(x = names(num_df), pattern = "absorbance_"))
    absorbances <- unlist(strsplit(names(num_df)[absorbance_index], "absorbance_"))
    absorbances <- absorbances[-which(absorbances == "")]

    absorbance <- c()
    for (i in 1:length(absorbance_index)) {
      if (length(which(is.na(num_df[, paste0("absorbance_", absorbances[i])]))) != 0) {
        absorbance[i] <- absorbances[i]
      }
    }
    rm(i)

    absorbance <- absorbance[!is.na(absorbance)]

    absorbance_index <- c()

    for (i in 1:length(absorbance)) {
      absorbance[i] <- paste0("absorbance_", absorbance[i])
      absorbance_index[i] <- which(names(num_df) == absorbance[i])
    }
    rm(i)
    rm(absorbance)
    rm(absorbances)
    c_absorbance <- length(absorbance_index)
  } else {
    absorbance_index <- NULL
    c_absorbance <- 0
  }

  if (is.null(absorbance_index)) {
    absorbance_index <- NA_integer_
    c_absorbance <- 0
  }

  #----Comparing Rel_clo with SPAD_605, if identical (NOT to be included into final df)----

  if (identical(df$Relative.Chlorophyll, df$SPAD_650)) {
    clo_index <- as.numeric(which(names(num_df) == "SPAD_650"))
    clo_rem <- 1
  } else {
    clo_index <- NA
    clo_rem <- 0
  }

  #----if NA's in spatial coordinates, removing the columns (NOT to be included into final df)----

  if (any(names(num_df) == "Latitude")) {
    if (any(is.na(num_df$Latitude))) {
      location_index <- c(which(names(num_df) == "Latitude"), which(names(num_df) == "Longitude"))
      c_lon_lat_a <- 2
      c_lon_lat <- 2
    } else {
      location_index <- NULL
      c_lon_lat <- 0
      c_lon_lat_a <- 0
    }
  } else {
    location_index <- NULL
    c_lon_lat <- 2
    c_lon_lat_a <- 0
  }

  f_indices <- unique(c(absorbance_index, clo_index, location_index))
  f_indices <- f_indices[!is.na(f_indices)]

  if (!installr::is.empty(f_indices)) {
    factor_df <- cbind(factor_df, num_df[, f_indices])
    num_df <- num_df[, -f_indices]
  }

  rm(absorbance_index, clo_index, f_indices)

  #-----Finding and removing rows with NA's-----
  cat("Finding and removing rows with NA's\n")
  row_NA <- list()


  for (i in 1:ncol(num_df)) {
    row_NA[[i]] <- num_df[, i] %>%
      is.na() %>%
      which()
  }
  rm(i)

  row_NA <- unique(unlist(row_NA))
  r_rem <- length(row_NA)

  #-----Removing PAR Phi2, PhiNPQ, PhiNO and Relative.Chlorophyll outliers-----

  cat("Removing PAR, Phi2 and Relative.Chlorophyll outliers\n")

  par_out <- which(num_df$Light.Intensity..PAR. > 2500 | num_df$Light.Intensity..PAR. < 1)
  phi2_out <- which(num_df$Phi2 > 0.85 | num_df$Phi2 < 0.03)
  # phinpq_out <- which(num_df$PhiNPQ > 0.85 | num_df$PhiNPQ < 0)
  # phino_out <- which(num_df$PhiNO > 0.5 | num_df$PhiNO < 0)

  if (clo_rem == 1) {
    rel_clo_out <- which(num_df$Relative.Chlorophyll > 75 | num_df$Relative.Chlorophyll < 0)
  } else if (clo_rem == 0) {
    rel_clo_out <- which(num_df$SPAD_650 > 75 | num_df$SPAD_650 < 0)
  }


  rows_out <- unique(c(issues_index, row_NA))
  out_param <- unique(c(par_out, phi2_out, rel_clo_out))
  final_removals <- unique(c(rows_out, out_param))
  r_out_param <- length(out_param)

  #-----Making final df-----

  cat("Making final df\n")
  if (!installr::is.empty(final_removals)) {
    num_df <- num_df[-final_removals, ]
    factor_df <- factor_df[-final_removals, ]
    rm_data <- df[final_removals, ]
  } else {
    rm_data <- NA
  }


  #-----Making summary table-----
  cat("Making summary table\n")
  summary_table <- as.data.frame(matrix(NA, 1))
  names(summary_table) <- "file_name"
  summary_table$file_name <- data_name ###### index with list.files()
  summary_table$total_obs <- r_full
  summary_table$total_variables <- c_full
  summary_table$non_numeric_variables_removed <- c_factor
  summary_table$absorbance_variables_removed_by_NA <- c_absorbance
  summary_table$LAT_LON_removed <- c_lon_lat
  summary_table$SPAD_650_removed <- clo_rem
  summary_table$obs_removed_by_Issues <- r_issue
  summary_table$obs_removed_by_NA <- r_rem
  summary_table$obs_removed_by_param_out <- r_out_param
  summary_table$Total_obs_removed <- length(final_removals)
  summary_table$perc_removed_obs <- round(((length(final_removals)) * 100 / r_full), digits = 2)
  summary_table$perc_removed_var <- round(((c_lon_lat_a + clo_rem + c_factor + c_absorbance) * 100 / c_full), digits = 1)


  #----Dropping empty levels in factors----
  num_df[, genotype] <- as.factor(num_df[, genotype])
  num_df <- droplevels(num_df)
  factor_df <- droplevels(factor_df)
  if (is.data.frame(rm_data)) {
    rm_data <- droplevels(rm_data)
  }

  #----output list----
  if (summary_table$perc_removed_obs != 0) {
    rm_table <- rm_data %>%
      dplyr::select_(.dots = c("date", SoV_names)) %>%
      table() %>%
      as.data.frame()
    if (any(rm_table$Freq == 0)) {
      rm_table <- rm_table[-which(rm_table[, "Freq"] == 0), ]
    }
    rm_table$date <- lubridate::ymd(rm_table$date)
    rm_table %<>% dplyr::arrange(desc(Freq))
    rm_table <- collapsibleTree::collapsibleTree(rm_table,
      hierarchy = names(rm_table),
      collapsed = TRUE, linkLength = 120,
      fontSize = 15
    )
  } else {
    rm_table <- "No observations were removed"
  }

  summary_table %<>% t %>% as.data.frame()
  names(summary_table) <- "value"
  summary_table$factor <- rownames(summary_table)
  rownames(summary_table) <- NULL
  summary_table %<>%
    dplyr::select(factor, value)
  summary_table <- collapsibleTree::collapsibleTree(summary_table,
    hierarchy = names(summary_table),
    collapsed = FALSE, zoomable = F,
    linkLength = 220,
    fontSize = 15
  )

  if (plotIm) {
    out <- list(num_df, factor_df, summary_table, rm_data, rm_table, SoV_names, genotype, miss, plots)
    names(out) <- c(
      "numeric_dataset", "non_numeric_dataset", "summary", "removed_observations", "removed_freq",
      "Sources_of_Variation", "Genotype", "imputed_variables", "imputed_plots"
    )
  } else {
    out <- list(num_df, factor_df, summary_table, rm_data, rm_table, SoV_names, genotype, miss)
    names(out) <- c(
      "numeric_dataset", "non_numeric_dataset", "summary", "removed_observations", "removed_freq", "Sources_of_Variation",
      "Genotype", "imputed_variables"
    )
  }

  if (!genotype %in% names(num_df)) {
    warning("The provided genotype argument does not match with the database names, the function MSPQ_ranks() might return an error", immediate. = T)
  }


  cat("Done!!!\n")

  return(out)
}
