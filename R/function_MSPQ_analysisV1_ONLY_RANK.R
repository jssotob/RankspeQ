#' MSPQ_ranks
#'
#' @param out The list obtained by MSPQ_tidy.
#' @param perIter Integer. Number of iterations to conduct on the Phi Index permutational analysis. 100 by default.
#' @param PerSeed Integer. Randomization seed for the Phi Index permutational analysis. 123 by default.
#' @param spats Logical. If spatial analysis is needed and row & column coordinates are provided. FALSE by default.
#' @param row Character. Name of the row coordinate variable if spats = TRUE. NULL by default.
#' @param column Character. Name of the column coordinate variable if spats = TRUE. NULL by default.
#' @param pl.date Logical. Will a sampling date be provided? FALSE by default.
#'
#' @return A list with 11 elements if spats = TRUE. 9 otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#' MSPQ_ranks()
#' }
MSPQ_ranks <- function(out, perIter = 100, PerSeed = 123,
                          spats = FALSE, row = NULL, column = NULL,
                          pl.date = FALSE){

  # Begining ----------------------------------------------------------------

  df_a <- out[["numeric_dataset"]]

  SoV <- c("date", "time", out$Sources_of_Variation)

  if(spats){
    if(any(is.null(row), is.null(column))){
      stop("Both row and column names must be provided. They are case sensitive")
    }
  }

  if(!spats & !is.null(c(row, column))){
    stop("row and column arguments must be removed from MSPQ_analysis() due to spats = FALSE")
  }

  if(!all(c(row,column) %in% names(df_a))){
    stop("The row/column names provided not found in the dataset. Please check spelling or change the arguments")
  }


  if(!spats){
    ind <- grep(pattern = "^col|^row|^fil", x = tolower(SoV))
    if(!installr::is.empty(ind)){
      SoV <- SoV[-ind]
      rm(ind)
    }
  }


    arr_SoV <- SoV[-c(grep("date", SoV),grep("time", SoV), grep(out$Genotype, SoV))]

  if(pl.date){
    multi <- (readline("Are there multiple planting dates? Y/N: "))
    multi <- toupper(multi)

    while(!any(c("Y", "YES", "N", "NO") %in% multi)){
      multi <- (readline("Wrong answer, please try again. Are there multiple planting dates? Y/N: "))
      multi <- toupper(multi)
    }

    if(any(c("N", "NO") %in% multi)){
      p.date <- suppressWarnings(lubridate::mdy((readline("Please type the planting date in format mm/dd/yyyy: "))))

      while(is.na(p.date)){
        p.date <- suppressWarnings(lubridate::mdy((readline("Wrong answer, check the date format. Please type the planting date in format mm/dd/yyyy: "))))
      }

      cat(paste0("The planting date provided is: ", p.date, '\n'))

      df_a %<>%
        mutate(planting_date = p.date,
               DAS = paste0(date-planting_date, " DAS"))

    } else{

      p.dates.factor <- (readline(paste0("Which of the '",
                                         paste0(arr_SoV, collapse = ", "),
                                         "' '", out$Genotype, "' was sown in multiple dates: ")))

      while(!any(c(arr_SoV,out$Genotype) %in% p.dates.factor)){
        p.dates.factor <- (readline(paste0("Wrong answer, please try again. Which of the '",
                                           paste0(arr_SoV, collapse = ", "), "' '",
                                           out$Genotype, "' was sown in multiple dates: ")))
      }

      p.dates <- list()

      for(a in 1:length(unique(df_a[,p.dates.factor]))){
        p.dates[[a]] <- (readline(paste0("Please type the planting date for ", unique(df_a[,p.dates.factor])[a], " in format mm/dd/yyyy: ")))
        p.dates[[a]] <- suppressWarnings(lubridate::mdy(p.dates[[a]]))

        while(is.na(p.dates[[a]])){
          p.dates[[a]] <- suppressWarnings(lubridate::mdy(readline(paste0("Wrong answer, check the date format. Please type the planting date for ", unique(df_a[,p.dates.factor])[a], " in format mm/dd/yyyy: "))))
        }
        cat(paste0("The planting date provided for ", unique(df_a[,p.dates.factor])[a],  " is: ", p.dates[[a]], '\n'))

        p.dates[[a]] <- data.frame(x = unique(df_a[,p.dates.factor])[a], planting_date = p.dates[[a]])
        names(p.dates[[a]])[1] <- p.dates.factor
      }
      rm(a)

      p.dates <- do.call(rbind, p.dates)
      df_a <- suppressMessages(left_join(df_a, p.dates) %>%
                                 mutate(DAS = paste0(date - planting_date, " DAS")))

    }
    SoV <- c("DAS", SoV)
  }



  cat("Descriptive tables \n")
  means <- df_a %>%              #Outputs
    as_tibble %>%
    group_by(.dots = SoV) %>%
    summarise_if(is.numeric, .funs = mean) %>%
    ungroup

  st_d <- df_a %>%               #Outputs
    as_tibble %>%
    group_by(.dots = SoV) %>%
    summarise_if(is.numeric, .funs = sd) %>%
    ungroup

  var_coef <- df_a %>%           #Outputs
    as_tibble %>%
    group_by(.dots = SoV) %>%
    summarise_if(is.numeric, .funs = EnvStats::cv) %>%
    ungroup

  medians <- df_a %>%            #Outputs
    as_tibble %>%
    group_by(.dots = SoV) %>%
    summarise_if(is.numeric, .funs = median) %>%
    ungroup


  # Phi_index_permutational -------------------------------------------------


  if(length(levels(out$numeric_dataset$time))==2){
  cat("Phi index permutational analysis \n")

  meanRatio <- function(x, jornada) {
    mr = mean(x[jornada == "AM"], trim = 0.025) / mean(x[jornada == "PM"], trim = 0.025)
    return(mr)
  }

  set.seed(PerSeed)

  spl <- df_a %>%
    dplyr::select_(.dots = c(SoV, "phi_index")) %>%
    group_by_(.dots = SoV[-which(SoV == "time")]) %>%
    group_split()

  odd_rat <- list()

  j <- length(spl)
  pbar <- txtProgressBar(min = 0,
                         max = j,
                         style = 3,
                         width = 50,
                         char = "=")

  for(j in 1:length(spl)){
    Dm <- numeric(length = perIter)
    N <- nrow(spl[[j]])

    for(i in seq_len(length(Dm) - 1)) {
      perm <- permute::shuffle(N)
      Dm[i] <- with(spl[[j]], meanRatio(phi_index, time[perm]))
    }


    rm(i)
    Dm[length(Dm)] <- with(spl[[j]], meanRatio(phi_index, time))

    D <- sum(Dm<Dm[length(Dm)])
    pv_perm <- D/length(Dm)
    odd_rat[[j]] <- spl[[j]] %>%
      select_(.dots = SoV[-2]) %>% .[1,] %>%
      mutate(phi_ratio = Dm[length(Dm)],
             Rless = D,
             p.value = pv_perm,
             Eval = pv_perm<0.05)
    setTxtProgressBar(pbar, j)
  }
  close(pbar)
  cat("\n")

  rm(j,spl)
  odd_rat <-do.call(rbind,odd_rat)



  if("DAS" %in% SoV){

    odd_out <- odd_rat %>%
      dplyr::filter(Eval == TRUE) %>%
      group_by_(.dots = SoV[4:length(SoV)]) %>%
      tally %>%
      arrange(desc(n))

  } else {

    odd_out <- odd_rat %>%
      dplyr::filter(Eval == TRUE) %>%
      group_by_(.dots = SoV[3:length(SoV)]) %>%
      tally %>%
      arrange(desc(n))
  }

  odd_out <- collapsibleTree::collapsibleTree(odd_out,
                                              hierarchy = names(odd_out),
                                              linkLength = 120)
  } else {
    odd_out <- "Procedure not done since time.dif = FALSE in MSPQ_tidy()"
  }
  # SPATS module ------------------------------------------------------------

  if(spats){
    cat("Adjusting variables with spatial components \n")
    if(all(c(row,column) %in% names(df_a))){

      group_vars <- c("date", "time", arr_SoV)

      spat_split <- as.list(group_split(medians, .dots = group_vars))

      spat_split <- lapply(spat_split, function(x)as.data.frame(x))

      varia <- medians %>% select_if(is.numeric) %>% names

      x <- c("LEF","NPQt","Phi2","PhiNO","PhiNPQ","PS1.Oxidized.Centers","FmPrime","FvP_over_FmP",
             "phi_index", "gH.","P700_DIRK_ampl","Leaf.Temperature.Differential","PS1.Active.Centers",
             "PS1.Over.Reduced.Centers", "Relative.Chlorophyll","vH.", "FoPrime", "Fs",
             "kP700", "tP700", "v_initial_P700", "value1")


      varia <- varia[varia %in% x]

      blups <- list()

      s <- length(spat_split)
      pbar <- txtProgressBar(min = 0,
                             max = s,
                             style = 3,
                             width = 50,
                             char = "=")

      for(s in 1:length(spat_split)){
        modelos <- lapply(varia, function(z){
          suppressWarnings(RankspeQ:::SpATS_mrbean(data = spat_split[[s]], response = z,genotype = out$Genotype, col = column, row =row, segm=T ,ncols=NULL, nrows=NULL, rep="",
                                        fix_fact = NULL, ran_fact = NULL, gen_ran = T, covariate = c("Leaf.Temperature", "Light.Intensity..PAR."),
                                        clean_out= FALSE, iterations=1, checks = NULL))
        })

        blups[[s]] <- do.call(cbind, lapply(modelos, function(x)data.frame(v = EnvStats::predict(x, which = out$Genotype, predFixed = "marginal")[,"predicted.values"])))
        names(blups[[s]]) <- varia

        blups[[s]] <- cbind(spat_split[[s]] %>% dplyr::select_(.dots = c(SoV, row, column)), blups[[s]])
        setTxtProgressBar(pbar, s)

      }
      close(pbar)
      cat("\n")
      rm(s)

      blups <- do.call(rbind,blups)

    } else {
      warning("The row/column names provided not found in the dataset, the procedure will continue with no spatial adjusted BLUP's")
    }

  }


  # Ranks -------------------------------------------------------------------

  cat("Ranking... \n")

  if(spats& "blups" %in% ls()){
    to_rank <- blups
    rm(blups)
  } else {
    to_rank <- medians
  }


  dates_rank <- lapply(as.character(unique(to_rank$date)), function(x)dplyr::filter(to_rank, date == x))

  x <- c("LEF","NPQt","Phi2","PhiNO","PhiNPQ","PS1.Oxidized.Centers","FmPrime","FvP_over_FmP",
         "phi_index", "gH.", "P700_DIRK_ampl","Leaf.Temperature.Differential","PS1.Active.Centers",
         "PS1.Over.Reduced.Centers", "Relative.Chlorophyll","vH.", "FoPrime", "Fs",
         "kP700", "tP700", "v_initial_P700", "value1")

  x <- names(df_a)[names(df_a) %in% x]


  ranks <- list()

  for(y in 1:length(x)){

    per_date <- list()

    for(i in 1:length(dates_rank)){

      per_date[[i]] <- as.list(dates_rank[[i]] %>%
                                 dplyr::select(all_of(c(SoV,x[y]))) %>%
                                 group_split(.dots = c("time", arr_SoV)))
    }
    rm(i)

    for(i in 1:length(per_date)){
      for(j in 1:length(per_date[[i]])){
        if(any(x[y] %in% c("LEF", "NPQt", "PhiNPQ", "PS1.Oxidized.Centers", "Vh.",
                           "v_initial_P700", "gH.", "P700_DIRK_ampl",
                           "kP700", "Leaf.Temperature.Differential"))){
          per_date[[i]][[j]] <- per_date[[i]][[j]][order(-per_date[[i]][[j]][,x[y]]),]

        } else {
          per_date[[i]][[j]] <- per_date[[i]][[j]][order(per_date[[i]][[j]][,x[y]]),]

        }
        per_date[[i]][[j]] %<>%
          droplevels %>%
          dplyr::select(-x[y]) %>%
          as_tibble %>%
          dplyr::mutate(x = 1:dplyr::n()) %>%
          #arrange_(.dots = SoV[grep("Gen|gen", SoV)]) %>%
          as.data.frame
      }
      per_date[[i]] <- do.call(rbind, per_date[[i]])
    }
    ranks[[y]] <- do.call(rbind, per_date)
  }

  ranks <- plyr::join_all(ranks, by= SoV, type='left')
  scores <- paste0(x,"_score")
  names(ranks)[grep("x", names(ranks))] <- scores
  ranks$final_score <- rowSums(ranks[,grep("_score", names(ranks))])


  # plots -------------------------------------------------------------------

  plots <- as.list(group_split(ranks, .dots = c("time",arr_SoV)))



  if("DAS" %in% SoV){
    plot_rank <- lapply(plots, function(s){
      plotly::ggplotly(s %>%
                 ggplot(aes_string(x = paste0("tidytext::reorder_within(",out$Genotype,",final_score, date)"), y = 'final_score')) +
                 geom_point() +
                 facet_wrap(~DAS, scales = "free_y", nrow = 2) +
                 coord_flip() +
                 tidytext::scale_x_reordered()+
                 theme(axis.text.x = element_text(angle = 90, hjust = 1),
                       axis.text.y = element_text(size = 5))+
                 labs(title = paste0(unique(s[,arr_SoV]), " ", as.character(unique(s$time))),
                      x = out$Genotype, y = "Final Score"))}
    )

    } else {

      plot_rank <- lapply(plots, function(s){
        plotly::ggplotly(s %>%
                   ggplot(aes_string(x = paste0("tidytext::reorder_within(",out$Genotype,",final_score, date)"), y = 'final_score')) +
                   geom_point() +
                   facet_wrap(~date, scales = "free_y", nrow = 2) +
                   coord_flip() +
                   tidytext::scale_x_reordered()+
                   theme(axis.text.x = element_text(angle = 90, hjust = 1),
                         axis.text.y = element_text(size = 5))+
                   labs(title = paste0(unique(s[,arr_SoV]), " ", as.character(unique(s$time)), collapse = "_"),
                        x = out$Genotype, y = "Final Score"))}
      )
    }

  cat("Making return \n")

  conam <- do.call(rbind, plots)

  if(spats){
      SPATS <- c(row, column)

    output <- list(means = means,
                   std = st_d,
                   var_coef = var_coef,
                   medians = medians,
                   BLUP_df = to_rank,
                   rank_table = conam,
                   rank_plots = plot_rank,
                   permutes = odd_out,
                   Sources_of_variation = SoV,
                   Genotype = out$Genotype,
                   SPATS_variables = SPATS)
  } else {
    output <- list(means = means,
                   std = st_d,
                   var_coef = var_coef,
                   medians = medians,
                   rank_table = conam,
                   rank_plots = plot_rank,
                   permutes = odd_out,
                   Sources_of_variation = SoV,
                   Genotype = out$Genotype)
  }
  cat("Done!!! \n")
  return(output)

}
