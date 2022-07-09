#' SpATS_mrbean
#'
#' A hidden function to run spatial analysis if the MultispeQ dataset contain row and column components
#'
#' @param data Raw data to be transformed
#' @param response Response variable. Automatically taken from data
#' @param genotype Character vector with the genotype column name. Automatically taken from Tidy output list
#' @param col Character vector with column information.
#' @param row Character vector with row information.
#' @param column column information
#' @param segm TRUE by default
#' @param ncols NULL by default
#' @param nrows NULL by default
#' @param rep empty by default
#' @param fix_fact NULL by default
#' @param ran_fact NULL by default
#' @param gen_ran TRUE by default
#' @param covariate c("Leaf.Temperature", "Light.Intensity..PAR.") by default
#' @param clean_out FALSE by default
#' @param iterations 1 by default
#' @param checks NULL by default
#'
#' @return list
#'
#' @examples
#' \dontrun{
#' SpATS_mrbean()
#' }

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

    Modelo=try(SpATS::SpATS(response=response,
                     genotype=genotype, genotype.as.random=gen_ran,
                     fixed= Fijo,
                     spatial = ~ SpATS::PSANOVA(col, row, nseg = c(ncols,nrows), degree=c(3,3),nest.div=2),
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

