#' Run an experiment with MIBench
#'
#' @param dgp A dgp function in the format of MIBench
#' @param MIalgorithm A multiple imputation algorithm in the format of MIBench
#' @param m The number of imputations
#' @param store_runs Set to TRUE if you want to store the output of the experiments. Default is FALSE.
#' @param suffix A suffix for the files, especially useful if you run multiple experiments in a loop.
#' @param ... Additional arguments for the dgp function
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
MIBench_experiment <-
  function(dgp = NULL,
            MIalgorithm = NULL,
            m = 10,
            store_runs = FALSE,
            suffix = 1,
            ...) {


    res <- dgp(...)
    store_seed <- attr(res, "seed")
    MI_res <- tryCatch(MIalgorithm(res$D_mis, m = m), error = function(e) list(imputations = NA, MI_name = MIalgorithm()$MI_name))

    res <- append(res, MI_res)
    attr(res, "seed") <- store_seed

    if (store_runs) {
      path <- here::here(paste0("experiments/", res$dgp_name, "/", MIalgorithm()$MI_name,
                     "/"))
      if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
      }

      saveRDS(res, paste0(path, "imputations_",
                          suffix, ".RDS"))
    }

    return(res)


  }


