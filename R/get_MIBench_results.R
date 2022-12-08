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

    res <- cbind(coef(tmp), confint(tmp))

    colnames(res) <- c("estimate", "2.5 %", "97.5 %")
    rownames(res) <- NULL
    return(res)
  },
  error = function(e) NA))

  infeasible <- lapply(obj, function(x)
    tryCatch({
    tmp <- x$analysis_model(x$D)

    res <- cbind(coef(tmp), confint(tmp))

    colnames(res) <- c("estimate", "2.5 %", "97.5 %")
    rownames(res) <- NULL
    return(res)
  },
  error = function(e) NA))


  results_uncongenial <-
    MIBench:::summarize_mi_analysis(Filter(Negate(anyNA), mi_uncongenial_combining_rules), obj[[1]]$true_values)

  results_congenial <-
    MIBench:::summarize_mi_analysis(Filter(Negate(anyNA), mi_congenial_combining_rules), obj[[1]]$true_values)

  results_lwd <-
    MIBench:::summarize_mi_analysis(Filter(Negate(anyNA), lwd), obj[[1]]$true_values)

  results_infeasible <-
    MIBench:::summarize_mi_analysis(Filter(Negate(anyNA), infeasible), obj[[1]]$true_values)

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
