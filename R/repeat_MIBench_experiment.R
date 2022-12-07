#' Repeat a MIBench_experiment multiple times
#'
#' @param dgp A dgp function in the format of MIBench
#' @param MIalgorithm A multiple imputation algorithm in the format of MIBench
#' @param m The number of imputations
#' @param store_runs Set to TRUE if you want to store the output of the experiments. Default is FALSE.
#' @param n_repetitions The number of repetitions (default is 1000) with fresh draws from the dgp function.
#' @param n_cores The number of cores for parallel processing of the experiments.
#' @param seed A random seed for the experiments. Note that a exact replication depends on the seed and the number of cores `n_cores`.
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

repeat_MIbench_experiment <-
  function(dgp,
           MIalgorithm,
           m = 10,
           store_runs = FALSE,
           n_repetitions = 1000,
           n_cores = 4,
           seed = NULL,
           ...) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)

    bettermc::mclapply(
      X = 1:n_repetitions,
      FUN = function(x)
        MIBench_experiment(
          dgp = dgp,
          MIalgorithm = MIalgorithm,
          m = m,
          suffix = x,
          store_runs = store_runs,
          ...
        ),
      mc.cores = n_cores,
      mc.allow.error = TRUE,
      mc.allow.fatal = TRUE,
      mc.cleanup = TRUE
    )
  }
