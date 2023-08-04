#' Comparisons between raknings and target trait of interest.
#'
#' @param ranks The list obtained by MSPQ_ranks.
#' @param target.trait.file data.frame in narrow format with the yield component or any other trait to contrast against MultispeQ ranking. An observation per experimental unit is required as well as the row & column coordinates if spats = TRUE in MSPQ_ranks.
#' @param target.trait.name Character. name of the yield variable on target.trait.file to contrast against MultispeQ ranking.
#' @param metadata data.frame in narrow format with extra information about the genotypes i.e. gene pool, known biotic/abiotic traits, etc. Optional process. NULL by default.
#'
#' @return A list with the confusion matrices and summary of predictions.
#' @export
#'
#' @examples
#' \dontrun{
#' # yield_data: any table containing the same evaluated genotypes and the yield or any target to contrast against MSPQ rankings
#'
#' target_trait_comparison(ranks = ranks,
#'                         target.trait.file = yield_data,
#'                         target.trait.name = "kg_ha",
#'                         metadata = NULL)
#' }
target_trait_comparison <- function(ranks, target.trait.file, target.trait.name, metadata = NULL){


# yield calibration -------------------------------------------------------


  if(!is.null(target.trait.file)){
    if(!is.data.frame(target.trait.file)){
      stop(paste0("The argument target.trait.file must be a data frame that contains the '", ranks$Genotype, "' column and the target trait to evaluate."))
    } else if(!ranks$Genotype %in% names(target.trait.file)){
      stop(paste0("The variable '", ranks$Genotype ,"' was not found into the target.trait.file dataset."))
    } else if(!all(unique(target.trait.file[,ranks$Genotype]) %in% unique(unlist(c(ranks$means[,ranks$Genotype]))))){
      gen_index <- unique(target.trait.file[,ranks$Genotype])[which(!unique(target.trait.file[,ranks$Genotype]) %in% unique(unlist(c(ranks$means[,ranks$Genotype]))))]
      gen_lets <- substr(gen_index, 1, 3)

      gen_real <- list()
      for(i in 1:length(gen_lets)){
      gen_real[[i]] <- unique(unlist(c(ranks$means[,ranks$Genotype])))[grep(gen_lets[i],unique(unlist(c(ranks$means[,ranks$Genotype]))))]
      }
      gen_real <- unlist(gen_real)

      if(!length(gen_real) == 0){
      stop(paste0("The genotype(s) '", paste0((gen_index), collapse = "', '"), "' from the target trait file not found into the ranked genotypes.
                  Should it be named as '",paste0((gen_real), collapse = "', '") , "'?"))
      } else{
        stop(paste0("The genotype(s) '", paste0((gen_index), collapse = "', '"), "' from the target trait file not found into the ranked genotypes."))
      }
    }
  }

  if("DAS" %in% ranks$Sources_of_variation){
    DAS_ind <- 4
  } else {
    DAS_ind <- 3
  }


  if(!is.null(ranks[["SPATS_variables"]])){
    if(!all(ranks[["SPATS_variables"]] %in% names(target.trait.file))){
      stop(paste0("The MulstispeQ ranks obtained by the function MSPQ_ranks() were adjusted by SpATS. The spatial variables ", paste0("'",ranks[["SPATS_variables"]],"'", collapse = ", "),
                " were not found into the target trait file"))
    } else if(!all(c(ranks$Sources_of_variation[DAS_ind:length(ranks$Sources_of_variation)], target.trait.name) %in% names(target.trait.file))){
      pos <- which(!c(ranks$Sources_of_variation[DAS_ind:length(ranks$Sources_of_variation)], target.trait.name) %in% names(target.trait.file))
      stop(paste0("The variable(s) ", paste0("'",c(ranks$Sources_of_variation[DAS_ind:length(ranks$Sources_of_variation)], target.trait.name)[pos],"'", collapse = ", "),
                  " do not exist into the target trait file. Make sure that it/they exist and try again"))
      rm(pos, DAS_ind)
    }

    cat("Adjusting target trait with spatial components... \n")


    if("DAS" %in% ranks$Sources_of_variation){
      sov_rm <- which(!ranks$Sources_of_variation %in% c(ranks$Genotype, "DAS", "date", "time"))
    } else {
      sov_rm <- which(!ranks$Sources_of_variation %in% c(ranks$Genotype, "date", "time"))
    }

    arr_SoV <- ranks$Sources_of_variation[sov_rm]
    rm(sov_rm)

    yield.list <- as.list(group_split(target.trait.file, .dots = arr_SoV))
    yield.list <- lapply(yield.list, function(x)as.data.frame(x))

    yield.model <- list()
    yield.blup <- list()

    for(i in 1:length(yield.list)){
      yield.model[[i]] <- suppressWarnings(SpATS_mrbean(data = yield.list[[i]], response = target.trait.name,genotype = ranks$Genotype, col = ranks$SPATS_variables[2], row =ranks$SPATS_variables[1], segm=T ,ncols=NULL, nrows=NULL, rep="",
                                                        fix_fact = NULL, ran_fact = NULL, gen_ran = T, covariate = NULL,
                                                        clean_out= FALSE, iterations=1, checks = NULL))
      yield.blup[[i]] <-data.frame(predict(yield.model[[i]], which = ranks$Genotype, predFixed = "marginal")[,c(ranks$Genotype,"predicted.values")])
      names(yield.blup[[i]])[2] <- target.trait.name
      yield.blup[[i]]$x = yield.list[[i]][,arr_SoV][1]
      for(s in 1:length(arr_SoV)){
        names(yield.blup[[i]])[2+s] <- arr_SoV[s]
      }
    }

    yield.transf <- do.call(rbind,yield.blup)
    yield.list <- yield.blup
    rm(i,s, DAS_ind)

  } else {
    if("DAS" %in% ranks$Sources_of_variation){
      sov_rm <- which(!ranks$Sources_of_variation %in% c(ranks$Genotype, "DAS", "date", "time"))
    } else {
      sov_rm <- which(!ranks$Sources_of_variation %in% c(ranks$Genotype, "date", "time"))
    }

    arr_SoV <- ranks$Sources_of_variation[sov_rm]
    rm(sov_rm)

    yield.transf <- target.trait.file %>%
      group_by_(.dots = c(ranks$Genotype, arr_SoV)) %>%
      summarise_if(.predicate = "is.numeric", .funs = "mean", na.rm=T) %>%
      ungroup %>% as.data.frame
}


# Confusion matrices ------------------------------------------------------


cat("Plotting confusion matrices... \n")

  conam <- ranks$rank_table

  if("DAS" %in% ranks$Sources_of_variation){
  conam %<>%
    dplyr::select_(.dots = c("date", "time", ranks$Genotype, arr_SoV, "cumulative_trait_score", "DAS")) %>%
    left_join(., yield.transf, by = c(ranks$Genotype, arr_SoV))

  conam <- as.list(group_split(conam,date, DAS, time, .dots = arr_SoV))
  } else {
    conam %<>%
      dplyr::select_(.dots = c("date", "time", ranks$Genotype, arr_SoV, "cumulative_trait_score")) %>%
      left_join(., yield.transf, by = c(ranks$Genotype, arr_SoV))

    conam <- as.list(group_split(conam,date, time, .dots = arr_SoV))
  }


  fin <- list()


    grupos <- 10
  falses <- round(grupos*0.3)

  for(i in 1:length(conam)){

    conam[[i]]$y.rank <- as.numeric(rank(conam[[i]][,target.trait.name]))
    conam[[i]]$m.rank <- as.numeric(rank(conam[[i]][,"cumulative_trait_score"]))
    conam[[i]] %<>% arrange(y.rank)
    clus <- round(nrow(conam[[i]])*0.1)
    conam[[i]]$clus.y <- rep(1:grupos,times = clus, length.out = nrow(conam[[i]])) %>% sort

    conam[[i]] %<>% arrange(m.rank)
    conam[[i]]$clus.m <- rep(1:grupos,times = clus, length.out = nrow(conam[[i]])) %>% sort


    conam[[i]] %<>%
      arrange(clus.y) %>%
      mutate(eval = case_when(clus.m== clus.y~TRUE,
                              clus.y+1 == clus.m ~ TRUE,
                              clus.y-1 == clus.m~TRUE,
                              # clus.y+2 == clus.m ~ TRUE,
                              # clus.y-2 == clus.m~TRUE,
                              TRUE~FALSE),
             col = case_when(eval == TRUE~"Predicted",
                             clus.m >= (max(clus.m)-falses)+1&clus.y <= falses~"False Positive",
                             clus.y >= (max(clus.m)-falses)+1&clus.m <= falses~"False Negative",
                             TRUE~"Low Prediction"))
  }
  rm(i)


  eval <- list()

  # if("DAS" %in% ranks$Sources_of_variation){
  #   for(i in 1:length(conam)){
  #     eval[[i]] <- DT::datatable(conam[[i]] %>%
  #       group_by(col) %>%
  #       dplyr::summarise(Count = dplyr::n()) %>%
  #       mutate(conf_matrix = paste0(unique(conam[[i]][,arr_SoV]), " ", conam[[i]]$DAS[1]," ", as.character(conam[[i]]$time[1]))) %>%
  #       dplyr::select(conf_matrix, col, Count) %>%
  #       dplyr::rename(Variable = col))
  #     names(eval)[[i]] <- eval[[i]]$x$data$conf_matrix[1]
  #   }
  # } else {
  #   for(i in 1:length(conam)){
  #     eval[[i]] <- DT::datatable(conam[[i]] %>%
  #                                  group_by(col) %>%
  #                                  dplyr::summarise(Count = dplyr::n()) %>%
  #                                  mutate(conf_matrix = paste0(unique(conam[[i]][,arr_SoV]), " ", conam[[i]]$date[1]," ", as.character(conam[[i]]$time[1]))) %>%
  #                                  dplyr::select(conf_matrix, col, Count) %>%
  #                                  dplyr::rename(Variable = col))
  #     names(eval)[[i]] <- eval[[i]]$x$data$conf_matrix[1]
  #   }
  # }

  if("DAS" %in% ranks$Sources_of_variation){
    for(i in 1:length(conam)){
      eval[[i]] <- conam[[i]] %>%
                         group_by(col) %>%
                         dplyr::summarise(Count = dplyr::n()) %>%
                         mutate(conf_matrix = paste0(unique(conam[[i]][,arr_SoV]), " ", conam[[i]]$DAS[1]," ", as.character(conam[[i]]$time[1]))) %>%
                         dplyr::select(conf_matrix, col, Count) %>%
                         dplyr::rename(Variable = col)

    }
    eval <- DT::datatable(do.call(rbind,eval),rownames = F)
  } else {
    for(i in 1:length(conam)){
      eval[[i]] <- conam[[i]] %>%
                        group_by(col) %>%
                        dplyr::summarise(Count = dplyr::n()) %>%
                        mutate(conf_matrix = paste0(unique(conam[[i]][,arr_SoV]), " ", conam[[i]]$date[1]," ", as.character(conam[[i]]$time[1]))) %>%
                        dplyr::select(conf_matrix, col, Count) %>%
                        dplyr::rename(Variable = col)
    }
    eval <- DT::datatable(do.call(rbind,eval),rownames = F)
  }

  conf_matrices <- list()

  for(i in 1:length(conam)){
    data_loop <- conam[[i]]
    aaa <- as.list(group_split(data_loop, clus.y, clus.m))
    aaa <- lapply(aaa, function(x){
      x$name <- paste0(unlist(c(x[,ranks$Genotype])), collapse = "\n")
      return(x)})
    aaa <- do.call(rbind, aaa) %>%
      dplyr::select(clus.y,clus.m, name, col) %>%
      .[!duplicated(.),]



    if("DAS" %in% ranks$Sources_of_variation){
    conf_matrices[[i]] <- plotly::ggplotly(ggplot2::ggplot(aaa,
                                                           ggplot2::aes_string(x = "factor(clus.y)", y = "factor(clus.m)", fill = "factor(col)",
                                                         text = "name"))+
                                           ggplot2::geom_raster(alpha = 0.95)+
                                           ggplot2::theme_classic()+
                                           ggplot2::labs(title = paste0(data_loop$DAS[1]," ", as.character(data_loop$time[1]), " ", unique(data_loop[,arr_SoV])),
                                          x = "Yield cluster", y = "MultispeQ cluster", fill = ""), tooltip="text")
    names(conf_matrices)[i] <- paste0(data_loop$DAS[1]," ", as.character(data_loop$time[1]), " ", unique(data_loop[,arr_SoV]))
    } else {
      conf_matrices[[i]] <- plotly::ggplotly(ggplot2::ggplot(aaa,
                                                             ggplot2::aes_string(x = "factor(clus.y)", y = "factor(clus.m)", fill = "factor(col)",
                                                           text = "name"))+
                                            ggplot2::geom_raster(alpha = 0.95)+
                                            ggplot2::theme_classic()+
                                            ggplot2::labs(title = paste0(data_loop$date[1]," ", as.character(data_loop$time[1]), " ", unique(data_loop[,arr_SoV])),
                                            x = "Yield cluster", y = "MultispeQ cluster", fill = ""), tooltip="text")
      names(conf_matrices)[i] <- paste0(data_loop$date[1]," ", as.character(data_loop$time[1]), " ", unique(data_loop[,arr_SoV]))
    }
  }
  rm(data_loop,i,aaa)


# metadata ----------------------------------------------------------------

   if(!is.null(metadata)){
    if(!is.data.frame(metadata)){
      warning(paste0("The argument metadata must be a data frame that contains the '", ranks$Genotype, "' column and the traits to evaluate. This procedure is not done"), immediate. = T)
    } else if(!ranks$Genotype %in% names(metadata)){
      warning(paste0("The variable '", ranks$Genotype ,"' was not found into the metadata dataset."), immediate. = T)
    } else if(!all(unique(metadata[,ranks$Genotype]) %in% unique(unlist(c(ranks$means[,ranks$Genotype]))))){
      gen_index <- unique(metadata[,ranks$Genotype])[which(!unique(metadata[,ranks$Genotype]) %in% unique(unlist(c(ranks$means[,ranks$Genotype]))))]
      gen_lets <- substr(gen_index, 1, 3)

      gen_real <- list()
      for(i in 1:length(gen_lets)){
        gen_real[[i]] <- unique(unlist(c(ranks$means[,ranks$Genotype])))[grep(gen_lets[i],unique(unlist(c(ranks$means[,ranks$Genotype]))))]
      }
      gen_real <- unlist(gen_real)

      if(!length(gen_real) == 0){
        warning(paste0("The genotype(s) '", paste0((gen_index), collapse = "', '"), "' from the metadata file not found into the ranked genotypes.
                  Should it be named as '",paste0((gen_real), collapse = "', '") , "'?. This procedure is not done"), immediate. = T)
      } else{
        warning(paste0("The genotype(s) '", paste0((gen_index), collapse = "', '"), "' from the metadata file not found into the ranked genotypes. This procedure is not done"), immediate. = T)
      }
    } else{

      metadatos <- suppressMessages(lapply(1:length(conam),function(x)left_join(conam[[x]], metadata)))

      cols <- names(metadatos[[1]])[which(!names(metadatos[[1]])%in% names(conam[[1]]))]

      cat(paste0("The genotypes metadata to analyse is: ", paste0(cols, collapse = " "), "\n"))

      metadatos <- do.call(rbind, metadatos)


      descriptors <- list()

      for(i in 1:length(cols)){
      descriptors[[i]] <- metadatos %>%
        group_by_(.dots = c(cols[i], arr_SoV)) %>%
        group_split %>% as.list
      }
      rm(i)

      if("DAS" %in% ranks$Sources_of_variation){

      for(i in 1:length(descriptors)){
        for(j in 1:length(descriptors[[i]])){
          suppressMessages(descriptors[[i]][[j]] %<>%
            group_by_(.dots = c("DAS", "date", arr_SoV)) %>%
            dplyr::summarise(Predicted = sum(col == "Predicted", na.rm = T),
                             No_Predicted = sum(col == "Low Prediction", na.rm = T),
                             False_Positive = sum(col == "False Positive", na.rm = T),
                             False_Negative = sum(col == "False Negative", na.rm = T)) %>%
            ungroup %>%
            tidyr::gather(., key = "x", value = "Number_of_Genotypes",Predicted, No_Predicted, False_Positive, False_Negative) %>%
            mutate(.,s = as.data.frame(descriptors[[i]][[j]])[1, cols[i]],
                   date = lubridate::ymd(date)))
          names(descriptors[[i]])[[j]] <- paste0(as.data.frame(descriptors[[i]][[j]])[1, "s"], " ", paste0(descriptors[[i]][[j]][1, arr_SoV], collapse = " "))
        }
        names(descriptors)[[i]] <- cols[i]
      }
      } else {
        for(i in 1:length(descriptors)){
          for(j in 1:length(descriptors[[i]])){
            suppressMessages(descriptors[[i]][[j]] %<>%
                               group_by_(.dots = c("date", arr_SoV)) %>%
                               dplyr::summarise(Predicted = sum(col == "Predicted", na.rm = T),
                                                No_Predicted = sum(col == "Low Prediction", na.rm = T),
                                                False_Positive = sum(col == "False Positive", na.rm = T),
                                                False_Negative = sum(col == "False Negative", na.rm = T)) %>%
                               ungroup %>%
                               tidyr::gather(., key = "x", value = "Number_of_Genotypes",Predicted, No_Predicted, False_Positive, False_Negative) %>%
                               mutate(.,s = as.data.frame(descriptors[[i]][[j]])[1, cols[i]],
                                      date = lubridate::ymd(date)))
            names(descriptors[[i]])[[j]] <- paste0(as.data.frame(descriptors[[i]][[j]])[1, "s"], " ", paste0(descriptors[[i]][[j]][1, arr_SoV], collapse = " "))
          }
          names(descriptors)[[i]] <- cols[i]
        }
      }



      meta_tables <- descriptors


      if("DAS" %in% ranks$Sources_of_variation){

      for(i in 1:length(descriptors)){
        for(j in 1:length(descriptors[[i]])){
      title <- paste0(as.data.frame(descriptors[[i]][[j]])[1, "s"], " ", paste0(descriptors[[i]][[j]][1, arr_SoV], collapse = " "))

      descriptors[[i]][[j]] <-descriptors[[i]][[j]] %>%
        ggplot2::ggplot(ggplot2::aes(x = DAS, y = Number_of_Genotypes,
                                     fill = factor(x, levels = c("Predicted","No_Predicted","False_Positive","False_Negative"))))+
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(1))+
        ggplot2::scale_fill_discrete(labels = c("Predicted","Low Prediction","False Positive","False Negative"))+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.title = ggplot2::element_blank())+
        ggplot2::labs(title = title, x = "Days After Sowing", y = "Occurrences (#)")

        }
      }
      } else {
        for(i in 1:length(descriptors)){
          for(j in 1:length(descriptors[[i]])){
            title <- paste0(as.data.frame(descriptors[[i]][[j]])[1, "s"], " ", paste0(descriptors[[i]][[j]][1, arr_SoV], collapse = " "))

            descriptors[[i]][[j]] <- plotly::ggplotly(descriptors[[i]][[j]] %>%
                                                      ggplot2::ggplot(ggplot2::aes(x = date,
                                                                                   y = Number_of_Genotypes,
                                                                                   fill = factor(x, levels = c("Predicted","No_Predicted","False_Positive","False_Negative"))))+
                                                        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(4.9))+
                                                        ggplot2::theme_bw()+
                                                        ggplot2::theme(legend.title = ggplot2::element_blank())+
                                                        ggplot2::labs(title = title, x = "Date", y = "Occurrences (#)"), tooltip = c("x", "y"))

          }
        }
    }



      ssss <- suppressMessages(do.call(rbind, conam) %>% left_join(.,metadata))

      # temporal, no commit
      aa <- ssss

      suppressMessages(ssss %<>%
        group_by_(.dots = c(ranks$Genotype, arr_SoV, cols)) %>%
        dplyr::summarise(Predicted = sum(col == "Predicted", na.rm = T),
                         Low_Prediction = sum(col == "Low Prediction", na.rm = T),
                         False_Positive = sum(col == "False Positive", na.rm = T),
                         False_Negative = sum(col == "False Negative", na.rm = T)) %>%
        ungroup())

      suppressMessages(dates <-  do.call(rbind, conam) %>%
        group_by_(.dots = ranks$Genotype, arr_SoV) %>%
        dplyr::summarise(Times = dplyr::n(),
                         Dates = Times/2) %>%
        ungroup())



      summary_table <- suppressMessages(DT::datatable(left_join(ssss, dates)))



    }
  } else{
    ssss <- do.call(rbind, conam)

    # temporal, no commit
    aa <- ssss

    suppressMessages(ssss %<>%
      group_by_(.dots = ranks$Genotype, arr_SoV) %>%
      dplyr::summarise(Predicted = sum(col == "Predicted", na.rm = T),
                       Low_Prediction = sum(col == "Low Prediction", na.rm = T),
                       False_Positive = sum(col == "False Positive", na.rm = T),
                       False_Negative = sum(col == "False Negative", na.rm = T)) %>%
      ungroup())

    suppressMessages(dates <-  do.call(rbind, conam) %>%
      group_by_(.dots = ranks$Genotype, arr_SoV) %>%
      dplyr::summarise(Times = dplyr::n(),
                       Dates = Times/2) %>%
      ungroup())



    summary_table <- suppressMessages(DT::datatable(left_join(ssss, dates)))
  }

  cat("Making return\n")

  if(!is.null(metadata)){
    output <- list(Conf_matrices = conf_matrices,
                   Summary_of_predictions = eval,
                   Summary_metadata = summary_table,
                   Metadata_tables = meta_tables,
                   Metadata_plots = descriptors,
                   tabla = conam) #temporal)
    if(!is.null(ranks[["SPATS_variables"]])){
      output$yield.BLUP <- yield.transf
    }
    cat("Done!!!\n")
  return(output)

  } else{
    output <- list(Conf_matrices = conf_matrices,
                   Summary_of_predictions = eval,
                   Summary_metadata = summary_table,
                   tabla = conam) #temporal)
    if(!is.null(ranks[["SPATS_variables"]])){
      output$yield.BLUP <- yield.transf
    }
    cat("Done!!!\n")
    return(output)
  }

 }



