#' Run a experiment with MIBench
#'
#' @param input_obj Either a list of class "lm" (the output from a call to lm) or a user provided list with the following entries beta_hat (the estimated regression coefficients),
#' varcov_hat (the estimated variance covariance matrix), sigma_hat (the estimated residual standard error), n (the number of observations) and k (the number of regression coefficients).
#' The list can be provided for more flexibility. Most users will call simulate directly on the output from a call to lm.
#' @param nsim_est Number of simulations to simulate estimation uncertainty (defaults to 1000).
#' @param nsim_fund Number of simulations to simulate fundamental uncertainty (defaults to 1000).
#' @param scenario Named list with values (scalar or vector, the vectors need to be of the same length) for each independent variable (names must match the names of the variables in the regression model) of scenarios for which predictions of the regression model should be calculated. (Default is at the mean for all variables only works when using a lm object.)
#' @param X Sometimes it is easier to directly pass the scenario as a model matrix (e.g. when you want to set the average of factors), this overrides any scenario specified under scenario. (Default is NULL)
#' @param observed_value_approach Do you want to use the observed value approach? (Default is FALSE.)
#' @param logged_dv Is the dependent variable in the linear model logged? (Default is TRUE.)
#' @param predicted_values Do you also want to output predicted values? (Default is FALSE. The resulting object can be very big.)
#' @param fast Do you want to speed up computation by not explicitly calculating predicted values? (Default is FALSE.)
#' @return A list of class "simloglm" with two matrices (geometric_mean and arithmetic_mean) with each nsim_est rows and number of scenarios columns and a matrix of predicted values if desired.
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
MIBench_light <-
  function (dgp = NULL,
            MIalgorithm = NULL,
            congenial = FALSE,
            n_iter = 10,
            compare = TRUE,
            n_data = 500,
            m = 10,
            seed = NULL,
            start_i = 1,
            store_runs = FALSE,
            load_runs = FALSE,
            algorithm_prefix = "none")
  {
    set.seed(seed)
    seed_list <- sample(2 ^ 20, size = n_iter)
    res_list <- vector("list", length = n_iter)
    if (store_runs) {
      path <- paste0("experiments/", dgp, "/", algorithm_prefix,
                     "/")
      if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
      }
    }
    if (compare)
      res_list_complete <- vector("list", length = n_iter)
    if (!load_runs) {
      cli::cli_progress_bar("Repeating the imputations", total = length(start_i:n_iter))
      for (i in start_i:n_iter) {
        if (dgp == "khjs-mar") {
          df <- MIBench:::dgp_fct(
            n_data = n_data,
            missingness = "mar1",
            dgp_name = "amelia",
            seed = seed_list[i]
          )
          analysis_model <- function(x) {
            lm(x[, 1] ~ x[, 2] + x[, 3])
          }
          true_values <- c(0,-0.11,-0.089)
        }
        if (dgp == "mixed-mcar") {
          df <- MIBench:::dgp_fct(
            n_data = n_data,
            missingness = "mcar1",
            dgp_name = "mixed",
            seed = seed_list[i]
          )
          analysis_model <- function(x) {
            lm(x[, 3] ~ x[, 1] * x[, 2])
          }
          true_values <- c(0, 1, 0,-2)
        }
        if (dgp == "cg-mar") {
          df <- MIBench:::dgp_fct(
            n_data = n_data,
            missingness = "mar3",
            dgp_name = "hd",
            seed = seed_list[i]
          )
          analysis_model <- function(x) {
            glm(x[, 1] ~ x[, 3] + x[, 4], family = binomial(link = logit))
          }
          true_values <- c(-1.92, 1.92, 1.92)
        }
        if (dgp == "mo-mcar") {
          df <- MIBench:::dgp_fct(
            n_data = n_data,
            missingness = "mcar1",
            dgp_name = "tbm",
            seed = seed_list[i]
          )
          analysis_model <- function(x) {
            lm(x[, 1] ~ x[, 2] + x[, 3] + x[, 4] + x[,
                                                     11])
          }
          true_values <- c(0, 1, 1, 1, 1)
        }
        if (compare) {
          model <- analysis_model(df$D_full)
          tmp <- cbind(coef(model), confint(model))
          colnames(tmp) <- c("estimate", "2.5 %", "97.5 %")
          rownames(tmp) <- 1:length(coef(model))
          res_list_complete[[i]] <- tmp
        }
        imputations <- MIBench:::quiet(MIalgorithm(df$D, m = m))
        #imputations <- lapply(imputations, function(x) sapply(x, as.numeric))

        if (store_runs) {
          saveRDS(imputations, paste0(path, "imputations_",
                                      i, ".RDS"))

        }
        cli::cli_progress_update()
      }
    }
    if (load_runs) {
      path <- paste0("experiments/", dgp, "/", algorithm_prefix,
                     "/")
      res_list <-
        lapply(start_i:n_iter, function(i)
          readRDS(paste0(path,
                         "analysis_", i, ".RDS")))
      if (dgp == "khjs-mar") {
        true_values <- c(0,-0.11,-0.089)
      }
      if (dgp == "mixed-mcar") {
        true_values <- c(0, 1, 0,-2)
      }
      if (dgp == "cg-mar") {
        true_values <- c(-1.92, 1.92, 1.92)
      }
      if (dgp == "mo-mcar") {
        true_values <- c(0, 1, 1, 1, 1)
      }
    }
    if (start_i == 1) {
      message("Please analyze the results in a separate step.")
    }
    else {
      message("Please analyze the results in a separate step.")
    }
  }
