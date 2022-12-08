#' Run all experiments for one multiple imputation algorithm
#'
#' @param MIalgorithm A multiple imputation algorithm in the format of MIBench
#' @param m The number of imputations
#' @param store_runs Set to TRUE if you want to store the imputations and all data of the experiments on disk. Default is FALSE.
#' @param store_results Set to TRUE if you want to store the summary of the results of the experiments on disk. Default is FALSE.
#' @param n_repetitions The number of repetitions (default is 1000) with fresh draws from the dgp function.
#' @param n_cores The number of cores for parallel processing of the experiments.
#' @param seed A random seed for the experiments. Note that a exact replication depends on the seed and the number of cores `n_cores`.
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

run_MIBench_all_dgps <-
  function(MIalgorithm,
           m = 10,
           store_runs = FALSE,
           store_results = FALSE,
           n_repetitions = 1000,
           n_cores = 4,
           seed = NULL) {
    # Update this list if new dgps are added
    all_dgps <- c(amelia_data,
                  hd_data,
                  marbach_data,
                  mixed_data,
                  tbm_data)

    res_list <- list()

    for (dgp in all_dgps) {
      for (mis in dgp()$missingness_patterns) {
        tmp <-
          repeat_MIbench_experiment(
            dgp = dgp,
            MIalgorithm = MIalgorithm,
            m = m,
            store_runs = store_runs,
            n_repetitions = n_repetitions,
            n_cores = n_cores,
            seed = seed,
            missingness = mis
          )

        res <- get_MIBench_results(tmp)

        res_list[[paste0(dgp(missingness = mis)$dgp_name, "_", MIalgorithm()$MI_name)]] <-
          res

        if (store_results) {
          path <- here::here(paste0("results/"))
          if (!dir.exists(path)) {
            dir.create(path, recursive = TRUE)
          }
          saveRDS(res,
                  paste0(path, dgp(missingness = mis)$dgp_name, "_", MIalgorithm()$MI_name, ".RDS"))
        }
      }
    }
    return(res_list)
  }
