#' yield_calibration
#'
#' @param ranks ranks
#' @param yield.file file
#' @param yield.name name
#' @param metadata metadata
#'
#' @return list
#' @export
#'
#' @examples
#' \dontrun{
#' yield_calibration()
#' }
yield_calibration <- function(ranks, yield.file, yield.name, metadata = NULL){


# yield calibration -------------------------------------------------------
# p_load(reshape, DT, plotly)

  if(!is.null(yield.file)){
    if(!is.data.frame(yield.file)){
      stop(paste0("The argument yield.file must be a data frame that contains the '", ranks$Genotype, "' column and the yield component to evaluate."))
    } else if(!ranks$Genotype %in% names(yield.file)){
      stop(paste0("The variable '", ranks$Genotype ,"' was not found into the yield.file dataset."))
    } else if(!all(unique(yield.file[,ranks$Genotype]) %in% unique(unlist(c(ranks$means[,ranks$Genotype]))))){
      gen_index <- unique(yield.file[,ranks$Genotype])[which(!unique(yield.file[,ranks$Genotype]) %in% unique(unlist(c(ranks$means[,ranks$Genotype]))))]
      gen_lets <- substr(gen_index, 1, 3)

      gen_real <- list()
      for(i in 1:length(gen_lets)){
      gen_real[[i]] <- unique(unlist(c(ranks$means[,ranks$Genotype])))[grep(gen_lets[i],unique(unlist(c(ranks$means[,ranks$Genotype]))))]
      }
      gen_real <- unlist(gen_real)

      if(!length(gen_real) == 0){
      stop(paste0("The genotype(s) '", paste0((gen_index), collapse = "', '"), "' from the yield file not found into the ranked genotypes.
                  Should it be named as '",paste0((gen_real), collapse = "', '") , "'?"))
      } else{
        stop(paste0("The genotype(s) '", paste0((gen_index), collapse = "', '"), "' from the yield file not found into the ranked genotypes."))
      }
    }
  }

  if("DAS" %in% ranks$Sources_of_variation){
    DAS_ind <- 4
  } else {
    DAS_ind <- 3
  }


  if(!is.null(ranks[["SPATS_variables"]])){
    if(!all(ranks[["SPATS_variables"]] %in% names(yield.file))){
      stop(paste0("The MulstispeQ ranks obtained by the function MSPQ_ranks() were adjusted by SpATS. The spatial variables ", paste0("'",ranks[["SPATS_variables"]],"'", collapse = ", "),
                " were not found into the yield file."))
    } else if(!all(c(ranks$Sources_of_variation[DAS_ind:length(ranks$Sources_of_variation)], yield.name) %in% names(yield.file))){
      pos <- which(!c(ranks$Sources_of_variation[DAS_ind:length(ranks$Sources_of_variation)], yield.name) %in% names(yield.file))
      stop(paste0("The variable(s) ", paste0("'",c(ranks$Sources_of_variation[DAS_ind:length(ranks$Sources_of_variation)], yield.name)[pos],"'", collapse = ", "),
                  " do not exist into the yield file. Make sure that it/they exist and try again"))
      rm(pos, DAS_ind)
    }

    cat("Adjusting yield with spatial components... \n")


    if("DAS" %in% ranks$Sources_of_variation){
      sov_rm <- which(!ranks$Sources_of_variation %in% c(ranks$Genotype, "DAS", "date", "time"))
    } else {
      sov_rm <- which(!ranks$Sources_of_variation %in% c(ranks$Genotype, "date", "time"))
    }

    arr_SoV <- ranks$Sources_of_variation[sov_rm]
    rm(sov_rm)

    yield.list <- as.list(group_split(yield.file, .dots = arr_SoV))
    yield.list <- lapply(yield.list, function(x)as.data.frame(x))

    yield.model <- list()
    yield.blup <- list()

    for(i in 1:length(yield.list)){
      yield.model[[i]] <- suppressWarnings(SpATS_mrbean(data = yield.list[[i]], response = yield.name,genotype = ranks$Genotype, col = ranks$SPATS_variables[2], row =ranks$SPATS_variables[1], segm=T ,ncols=NULL, nrows=NULL, rep="",
                                                        fix_fact = NULL, ran_fact = NULL, gen_ran = T, covariate = NULL,
                                                        clean_out= FALSE, iterations=1, checks = NULL))
      yield.blup[[i]] <-data.frame(predict(yield.model[[i]], which = ranks$Genotype, predFixed = "marginal")[,c(ranks$Genotype,"predicted.values")])
      names(yield.blup[[i]])[2] <- yield.name
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

    yield.transf <- yield.file %>%
      group_by_(.dots = c(ranks$Genotype, arr_SoV)) %>%
      summarise_if(.predicate = "is.numeric", .funs = "mean", na.rm=T) %>%
      ungroup %>% as.data.frame
  }


# Confusion matrices ------------------------------------------------------


cat("Plotting confusion matrices... \n")

  conam <- ranks$rank_table

  if("DAS" %in% ranks$Sources_of_variation){
  conam %<>%
    dplyr::select_(.dots = c("date", "time", ranks$Genotype, arr_SoV, "final_score", "DAS")) %>%
    left_join(., yield.transf, by = c(ranks$Genotype, arr_SoV))

  conam <- as.list(group_split(conam,date, DAS, time, .dots = arr_SoV))
  } else {
    conam %<>%
      dplyr::select_(.dots = c("date", "time", ranks$Genotype, arr_SoV, "final_score")) %>%
      left_join(., yield.transf, by = c(ranks$Genotype, arr_SoV))

    conam <- as.list(group_split(conam,date, time, .dots = arr_SoV))
  }


  fin <- list()


    grupos <- 10
  falses <- round(grupos*0.3)

  for(i in 1:length(conam)){

    conam[[i]]$y.rank <- as.numeric(rank(conam[[i]][,yield.name]))
    conam[[i]]$m.rank <- as.numeric(rank(conam[[i]][,"final_score"]))
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
                              TRUE~FALSE),
             col = case_when(eval == TRUE~"Predicted",
                             clus.m >= (max(clus.m)-falses)+1&clus.y <= falses~"False positive",
                             clus.y >= (max(clus.m)-falses)+1&clus.m <= falses~"False negative",
                             TRUE~"Low prediction"))
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
    conf_matrices[[i]] <- ggplotly(ggplot(aaa,aes_string(x = "factor(clus.y)", y = "factor(clus.m)", fill = "factor(col)",
                                                         text = "name"))+
                                     geom_raster(alpha = 0.95)+
                                     theme_classic()+
                                     labs(title = paste0(data_loop$DAS[1]," ", as.character(data_loop$time[1]), " ", unique(data_loop[,arr_SoV])),
                                          x = "Yield cluster", y = "MultispeQ cluster", fill = ""), tooltip="text")
    names(conf_matrices)[i] <- paste0(data_loop$DAS[1]," ", as.character(data_loop$time[1]), " ", unique(data_loop[,arr_SoV]))
    } else {
      conf_matrices[[i]] <- ggplotly(ggplot(aaa,aes_string(x = "factor(clus.y)", y = "factor(clus.m)", fill = "factor(col)",
                                                           text = "name"))+
                                       geom_raster(alpha = 0.95)+
                                       theme_classic()+
                                       labs(title = paste0(data_loop$date[1]," ", as.character(data_loop$time[1]), " ", unique(data_loop[,arr_SoV])),
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
                             No_Predicted = sum(col == "Low prediction", na.rm = T),
                             False_Positive = sum(col == "False positive", na.rm = T),
                             False_Negative = sum(col == "False negative", na.rm = T)) %>%
            ungroup %>%
            gather(., key = "x", value = "Number_of_Genotypes",Predicted, No_Predicted, False_Positive, False_Negative) %>%
            mutate(.,s = as.data.frame(descriptors[[i]][[j]])[1, cols[i]],
                   date = ymd(date)))
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
                                                No_Predicted = sum(col == "Low prediction", na.rm = T),
                                                False_Positive = sum(col == "False positive", na.rm = T),
                                                False_Negative = sum(col == "False negative", na.rm = T)) %>%
                               ungroup %>%
                               gather(., key = "x", value = "Number_of_Genotypes",Predicted, No_Predicted, False_Positive, False_Negative) %>%
                               mutate(.,s = as.data.frame(descriptors[[i]][[j]])[1, cols[i]],
                                      date = ymd(date)))
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
        ggplot(aes(x = DAS, y = Number_of_Genotypes, fill = factor(x, levels = c("Predicted","No_Predicted","False_Positive","False_Negative"))))+
        geom_bar(stat = "identity", position = position_dodge(1))+
        theme_bw()+
        theme(legend.title = element_blank())+
        labs(title = title, x = "Days After Sowing", y = "Ocurrences (#)")

        }
      }
      } else {
        for(i in 1:length(descriptors)){
          for(j in 1:length(descriptors[[i]])){
            title <- paste0(as.data.frame(descriptors[[i]][[j]])[1, "s"], " ", paste0(descriptors[[i]][[j]][1, arr_SoV], collapse = " "))

            descriptors[[i]][[j]] <- ggplotly(descriptors[[i]][[j]] %>%
                                                ggplot(aes(x = date, y = Number_of_Genotypes, fill = factor(x, levels = c("Predicted","No_Predicted","False_Positive","False_Negative"))))+
                                                geom_bar(stat = "identity", position = position_dodge(4.9))+
                                                theme_bw()+
                                                theme(legend.title = element_blank())+
                                                labs(title = title, x = "Date", y = "Ocurrences (#)"), tooltip = c("x", "y"))

          }
        }
    }



      ssss <- suppressMessages(do.call(rbind, conam) %>% left_join(.,metadata))

      suppressMessages(ssss %<>%
        group_by_(.dots = c(ranks$Genotype, arr_SoV, cols)) %>%
        dplyr::summarise(Predicted = sum(col == "Predicted", na.rm = T),
                         No_Predicted = sum(col == "Low prediction", na.rm = T),
                         False_Positive = sum(col == "False positive", na.rm = T),
                         False_Negative = sum(col == "False negative", na.rm = T)) %>%
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

    suppressMessages(ssss %<>%
      group_by_(.dots = ranks$Genotype, arr_SoV) %>%
      dplyr::summarise(Predicted = sum(col == "Predicted", na.rm = T),
                       No_Predicted = sum(col == "Low prediction", na.rm = T),
                       False_Positive = sum(col == "False positive", na.rm = T),
                       False_Negative = sum(col == "False negative", na.rm = T)) %>%
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
    output <- list(Conf_matrices = conf_matrices, Summary_of_predictions = eval, Summary_metadata = summary_table,
                Metadata_tables = meta_tables, Metadata_plots = descriptors)
    if(!is.null(ranks[["SPATS_variables"]])){
      output$yield.BLUP <- yield.transf
    }
    cat("Done!!!\n")
  return(output)

  } else{
    output <- list(Conf_matrices = conf_matrices, Summary_of_predictions = eval, Summary_metadata = summary_table)
    if(!is.null(ranks[["SPATS_variables"]])){
      output$yield.BLUP <- yield.transf
    }
    cat("Done!!!\n")
    return(output)
  }

 }



