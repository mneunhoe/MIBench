#' Get the results from the output of repeat_MIBench_experiment
#'
#' @param obj The output of a run of `repeat_MIBench_experiment`
#' @return A list of class "MIbench_imputations"
#' @export
#'
#' @examples
#' df <- cars
#' regression <- lm(log(dist)~speed, data = df)
#' # Specifiying no scenario to simulate at the mean of speed.
#' simloglm(regression)
#' # Explicitily specifying a scenario.
#' simloglm(regression, scenario = list(speed = c(5, 10, 20)))
#'
get_MIBench_results <- function(obj) {
  mi_uncongenial_combining_rules <- lapply(obj, function(x)
    tryCatch(
      MIBench:::analyze_mi(x$imputations, x$analysis_model, congenial = FALSE),
      error = function(e)
        NA
    ))

  mi_congenial_combining_rules <- lapply(obj, function(x)
    tryCatch(
      MIBench:::analyze_mi(x$imputations, x$analysis_model, congenial = TRUE),
      error = function(e)
        NA
    ))

  lwd <- lapply(obj, function(x)
    tryCatch({
    tmp <- x$analysis_model(x$D_mis)

    tmp_ci <- confint(tmp)

    if(length(dim(tmp_ci)>2)){
      tmp_ci1 <- NULL
      for(i in 1:dim(tmp_ci)[3]){
        tmp_ci1 <- rbind(tmp_ci1, tmp_ci[,,i])
      }

      tmp_ci <- tmp_ci1[c(matrix(1:dim(tmp_ci1)[1], nrow = dim(tmp_ci)[3], byrow = T)),]

    }

    res <- cbind(as.numeric(coef(tmp)), tmp_ci)

    colnames(res) <- c("estimate", "2.5 %", "97.5 %")
    rownames(res) <- NULL
    return(res)
  },
  error = function(e) NA))

  infeasible <- lapply(obj, function(x)
    tryCatch({
    tmp <- x$analysis_model(x$D)

    tmp_ci <- confint(tmp)

    if(length(dim(tmp_ci)>2)){
      tmp_ci1 <- NULL
      for(i in 1:dim(tmp_ci)[3]){
        tmp_ci1 <- rbind(tmp_ci1, tmp_ci[,,i])
      }

      tmp_ci <- tmp_ci1[c(matrix(1:dim(tmp_ci1)[1], nrow = dim(tmp_ci)[3], byrow = T)),]

    }

    res <- cbind(as.numeric(coef(tmp)), tmp_ci)

    colnames(res) <- c("estimate", "2.5 %", "97.5 %")
    rownames(res) <- NULL
    return(res)
  },
  error = function(e) NA))


  results_uncongenial <-
    MIBench:::summarize_mi_analysis(Filter(Negate(anyNA), mi_uncongenial_combining_rules), as.numeric(obj[[1]]$true_values))

  results_congenial <-
    MIBench:::summarize_mi_analysis(Filter(Negate(anyNA), mi_congenial_combining_rules), as.numeric(obj[[1]]$true_values))

  results_lwd <-
    MIBench:::summarize_mi_analysis(Filter(Negate(anyNA), lwd), as.numeric(obj[[1]]$true_values))

  results_infeasible <-
    MIBench:::summarize_mi_analysis(Filter(Negate(anyNA), infeasible), as.numeric(obj[[1]]$true_values))

  res <- list(
    congenial = results_congenial,
    uncongenial = results_uncongenial,
    lwd = results_lwd,
    infeasible = results_infeasible,
    number_of_runs = length(obj),
    number_of_failed_runs = sum(sapply(mi_uncongenial_combining_rules, anyNA))
  )

  return(res)
}
