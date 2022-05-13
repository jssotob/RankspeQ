#' MSPQ_ranks
#'
#' @param out output
#' @param perIter iterations
#' @param PerSeed seed
#' @param spats TRUE or FALSE if spatial analysis is needed
#' @param row row information
#' @param column column information
#' @param pl.date date
#'
#' @return list
#' @export
#'
#' @examples
#' \dontrun{
#' MSPQ_ranks()
#' }
MSPQ_ranks <- function(out, perIter = 100, PerSeed = 123,
                          spats = FALSE, row = NULL, column = NULL,
                          pl.date = FALSE){

  # library(pacman)
  # p_load(tidyverse, magrittr, lmerTest, lme4, permute, vegan,
  #        collapsibleTree, EnvStats, SpATS, lubridate,cvms, broom, tibble, moments, tidytext)


  # SPATS_mrbean ------------------------------------------------------------

  SpATS_mrbean <- function(data, response ,genotype,
                           col, row, segm ,ncols, nrows, rep,
                           fix_fact, ran_fact, gen_ran, covariate,
                           clean_out= FALSE, iterations=1, checks = NULL){
    dt <- data
    dt[ , genotype] <- as.factor(dt[ , genotype])
    dt$col =  dt[,col]
    dt$row = dt[,row]
    dt$col_f = factor( dt[,col])
    dt$row_f = factor( dt[,row])
    ncols = ifelse(is.null(ncols)|isFALSE(segm),length(unique(dt[,col])), ncols )
    nrows = ifelse(is.null(nrows)|isFALSE(segm),length(unique(dt[,row])), nrows )
    for (i in 1:length(fix_fact)) {
      dt[, fix_fact[i]] <- as.factor(dt[, fix_fact[i]])
    }
    for (i in 1:length(ran_fact)) {
      dt[, ran_fact[i]] <- as.factor(dt[, ran_fact[i]])
    }
    if(is.null(fix_fact)&is.null(covariate)) {Fijo <- as.formula(~NULL)
    } else if(!is.null(fix_fact)&is.null(covariate)){
      Fijo <- as.formula(paste("", paste(fix_fact, collapse=" + "), sep=" ~ "))
    } else if(!is.null(fix_fact)&!is.null(covariate)){
      Fijo <- as.formula(paste("", paste(c(fix_fact,covariate), collapse=" + "), sep=" ~ "))
    } else if(is.null(fix_fact)&!is.null(covariate)) {
      Fijo <- as.formula(paste("",paste(covariate, collapse = " + "),sep = " ~ "))
    }

    #--
    if(rep!=""){
      dt$rep <- as.factor(dt[, rep])
      if(nlevels(dt$rep)>1){
        if(is.null(ran_fact)) Random <- as.formula(~ rep:col_f + rep:row_f)
        else Random <-  as.formula(paste("" ,paste(c(ran_fact,"rep:col_f","rep:row_f"), collapse=" + "), sep=" ~ "))
        if(is.null(fix_fact)&is.null(covariate)){
          Fijo <- as.formula(~ rep)
        } else {
          Fijo <-  paste(paste0(as.character(Fijo), collapse = ""), "+ rep" )
        }
      } else {
        if(is.null(ran_fact)) Random <- as.formula(~ col_f+row_f)
        else Random <-  as.formula(paste("" ,paste(c(ran_fact,"col_f","row_f"), collapse=" + "), sep=" ~ "))
      }
    } else {
      if(is.null(ran_fact)) Random <- as.formula(~ col_f+row_f)
      else Random <-  as.formula(paste("" ,paste(c(ran_fact,"col_f","row_f"), collapse=" + "), sep=" ~ "))
    }

    # if(is.null(ran_fact)) Random <- as.formula(~ col_f+row_f)
    # else Random <-  as.formula(paste("" ,paste(c(ran_fact,"col_f","row_f"), collapse=" + "), sep=" ~ "))

    if(!is.null(checks)& gen_ran == T){
      dt <- check_gen_SpATS(gen = genotype, data = dt, check_gen = checks )
      if("checks" %in% names(dt) ){
        if(is.null(fix_fact) & is.null(covariate) & rep==""){
          Fijo <- as.formula(~ checks)
        } else {
          Fijo <-  paste(paste0(as.character(Fijo), collapse = ""), "+ checks" )
        }
      }
    }

    Modelo=try(SpATS(response=response,
                     genotype=genotype, genotype.as.random=gen_ran,
                     fixed= Fijo,
                     spatial = ~ PSANOVA(col, row, nseg = c(ncols,nrows), degree=c(3,3),nest.div=2),
                     random = Random, data=dt,
                     control = list(tolerance=1e-03, monitoring=0)),silent = T)
    # tryCatch(
    #   {
    #     if(class(Modelo)=="try-error") stop("Error in the components of model")
    #   },
    #   error = function(e) {
    #     shinytoastr::toastr_error(title = "Warning:", conditionMessage(e),position =  "bottom-right",progressBar = TRUE)
    #   }
    # )
    # if(class(Modelo)=="try-error") return()

    if(clean_out){
      resum_out <- msa_residuals(Modelo)
      if(resum_out>0){
        tmp_out <- 1
        counter <- 1
        while (tmp_out>0 & counter<=iterations) {
          c_datos <- res_raw_data(Modelo)
          c_datos[, response] <- ifelse(c_datos$Classify=="Outlier", NA, c_datos[, response] )
          c_datos <- c_datos %>% dplyr::select(-weights)
          Modelo = try(SpATS(response=response,
                             genotype=genotype, genotype.as.random=gen_ran,
                             fixed= Fijo,
                             spatial = ~ PSANOVA(col, row, nseg = c(ncols,nrows), degree=c(3,3),nest.div=2),
                             random = Random, data=c_datos,
                             control = list(tolerance=1e-03, monitoring=0)), silent = T)
          tmp_out <- msa_residuals(Modelo)
          if(iterations>1) resum_out <-  resum_out + tmp_out
          counter <- counter + 1
        }
      }
    }

    Modelo
  }


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


  # if(spats){
  #   ind <- which(SoV %in% c(row, column))
  #   SoV <- SoV[-ind]
  #   rm(ind)
  # }

  if(!spats){
    ind <- grep(pattern = "^col|^row|^fil", x = tolower(SoV))
    if(!is.empty(ind)){
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
      p.date <- suppressWarnings(mdy((readline("Please type the planting date in format mm/dd/yyyy: "))))

      while(is.na(p.date)){
        p.date <- suppressWarnings(mdy((readline("Wrong answer, check the date format. Please type the planting date in format mm/dd/yyyy: "))))
      }

      cat(paste0("The planting date provided is: ", p.date, '\n'))

      df_a %<>%
        mutate(planting_date = p.date,
               DAS = paste0(date-planting_date, " DAS"))

    } else{

      p.dates.factor <- (readline(paste0("Which of the '", paste0(arr_SoV, collapse = ", "), "' '", out$Genotype, "' was sown in multiple dates: ")))
      p.dates <- list()

      for(a in 1:length(unique(df_a[,p.dates.factor]))){
        p.dates[[a]] <- (readline(paste0("Please type the planting date for ", unique(df_a[,p.dates.factor])[a], " in format mm/dd/yyyy: ")))
        p.dates[[a]] <- suppressWarnings(mdy(p.dates[[a]]))

        while(is.na(p.dates[[a]])){
          p.dates[[a]] <- suppressWarnings(mdy(readline(paste0("Wrong answer, check the date format. Please type the planting date for ", unique(df_a[,p.dates.factor])[a], " in format mm/dd/yyyy: "))))
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
    summarise_if(.predicate = "is.numeric", .funs = "mean") %>%
    ungroup

  st_d <- df_a %>%               #Outputs
    as_tibble %>%
    group_by(.dots = SoV) %>%
    summarise_if(.predicate = "is.numeric", .funs = "sd") %>%
    ungroup

  var_coef <- df_a %>%           #Outputs
    as_tibble %>%
    group_by(.dots = SoV) %>%
    summarise_if(.predicate = "is.numeric", .funs = "cv") %>%
    ungroup

  medians <- df_a %>%            #Outputs
    as_tibble %>%
    group_by(.dots = SoV) %>%
    summarise_if(.predicate = "is.numeric", .funs = "median") %>%
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
  pbar <- create_progress_bar('text')
  pbar$init(j)

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
    pbar$step()
  }
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
      pbar <- create_progress_bar('text')
      pbar$init(s)

      for(s in 1:length(spat_split)){
        modelos <- lapply(varia, function(z){
          suppressWarnings(SpATS_mrbean(data = spat_split[[s]], response = z,genotype = out$Genotype, col = column, row =row, segm=T ,ncols=NULL, nrows=NULL, rep="",
                                        fix_fact = NULL, ran_fact = NULL, gen_ran = T, covariate = c("Leaf.Temperature", "Light.Intensity..PAR."),
                                        clean_out= FALSE, iterations=1, checks = NULL))
        })

        blups[[s]] <- do.call(cbind, lapply(modelos, function(x)data.frame(v = EnvStats::predict(x, which = out$Genotype, predFixed = "marginal")[,"predicted.values"])))
        names(blups[[s]]) <- varia

        blups[[s]] <- cbind(spat_split[[s]] %>% dplyr::select_(.dots = c(SoV, row, column)), blups[[s]])
        pbar$step()

      }

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

  ranks <- join_all(ranks, by= SoV, type='left')
  scores <- paste0(x,"_score")
  names(ranks)[grep("x", names(ranks))] <- scores
  ranks$final_score <- rowSums(ranks[,grep("_score", names(ranks))])


  # plots -------------------------------------------------------------------

  plots <- as.list(group_split(ranks, .dots = c("time",arr_SoV)))



  if("DAS" %in% SoV){
    plot_rank <- lapply(plots, function(s){
      ggplotly(s %>%
                 ggplot(aes_string(x = paste0("reorder_within(",out$Genotype,",final_score, date)"), y = 'final_score')) +
                 geom_point() +
                 facet_wrap(~DAS, scales = "free_y", nrow = 2) +
                 coord_flip() +
                 scale_x_reordered()+
                 theme(axis.text.x = element_text(angle = 90, hjust = 1),
                       axis.text.y = element_text(size = 5))+
                 labs(title = paste0(unique(s[,arr_SoV]), " ", as.character(unique(s$time))),
                      x = out$Genotype))}
    )

    } else {

      plot_rank <- lapply(plots, function(s){
        ggplotly(s %>%
                   ggplot(aes_string(x = paste0("reorder_within(",out$Genotype,",final_score, date)"), y = 'final_score')) +
                   geom_point() +
                   facet_wrap(~date, scales = "free_y", nrow = 2) +
                   coord_flip() +
                   tidytext::scale_x_reordered()+
                   theme(axis.text.x = element_text(angle = 90, hjust = 1),
                         axis.text.y = element_text(size = 5))+
                   labs(title = paste0(unique(s[,arr_SoV]), " ", as.character(unique(s$time)), collapse = "_"),
                        x = out$Genotype))}
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
