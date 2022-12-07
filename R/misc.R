#' @importFrom magrittr "%>%"
# Function to generate data according to King et al 2001

kropko_data <- function(n = 1000, missingness = .25, ...) {
  store_seed <- .Random.seed

  .dgp <-
    function(N,
             N_FULL = 3L,
             N_PARTIAL = 1,
             restrict = "triangular",
             ncat = 3,
             type = "continuous",
             pr_miss = .25,
             imp_meth = "ppd",
             strong = 0) {
      if (type != "continuous" &
          type != "binary" &
          type != "ordinal" &
          type != "nominal")
        stop("Type must be continuous, binary, ordinal, or nominal")
      cpc <- (restrict == "none")

      rdf.names <-
        c(paste("x_", 1:(N_FULL + N_PARTIAL), sep = ""), "y_1")

      if (type == "continuous" |
          type == "binary")
        rdf <-
        mi::rdata.frame(
          N = N,
          n_full = N_FULL,
          n_partial = (N_PARTIAL + 1),
          restrictions = restrict,
          types = c(rep("continuous", (
            N_FULL + N_PARTIAL
          )), type),
          pr_miss = c(rep(.1, N_PARTIAL), .25),
          strong = strong,
          estimate_CPCs = cpc
        )
      if (type == "ordinal" |
          type == "nominal")
        rdf <-
        mi::rdata.frame(
          N = N,
          n_full = N_FULL,
          n_partial = (N_PARTIAL + 1),
          restrictions = restrict,
          types = c(rep("continuous", (
            N_FULL + N_PARTIAL
          )), type),
          pr_miss = c(rep(.1, N_PARTIAL), .25),
          n_cat = ncat,
          strong = strong,
          estimate_CPCs = cpc
        )

      if (cpc)
        nmar <- sqrt(mean(rdf$empirical_CPCs ^ 2))
      else
        nmar <- NA
      data <- rdf$obs

      true <- rdf$true
      colnames(data) <- colnames(true) <- rdf.names

      if (type == "continuous") {
        analysis_model <- function(x) {
          arm::bayesglm(y_1 ~ x_1 + x_2 + x_3 + x_4, data = x)
        }
      }

      if (type == "binary") {
        analysis_model <- function(x) {
          arm::bayesglm(y_1 ~ x_1 + x_2 + x_3 + x_4,
                        family = binomial(link = "logit"),
                        data = x)
        }
      }

      if (type == "ordinal") {
        analysis_model <- function(x) {
          arm::bayespolr(y_1 ~ x_1 + x_2 + x_3 + x_4,
                         drop.unused.levels = FALSE,
                         data = x)
        }
      }


      if (type == "nominal") {
        analysis_model <- function(x) {
          nnet::multinom(y_1 ~ x_1 + x_2 + x_3 + x_4,
                         data = x,
                         maxit = 1000)
        }
      }


      return(list(true,
                  data,
                  analysis_model,
                  type))

    }

  tmp <- .dgp(N = n, pr_miss = missingness, ...)


  res <- list(
    D = tmp[[1]],
    D_mis = tmp[[2]],
    analysis_model = tmp[[3]],
    true_values = NA,
    dgp_name = paste0("kropko_data_", tmp[[4]], "_", gsub("\\.", "", paste(missingness))),
    missingness_patterns = c(0.25)
  )
  attr(res, "seed") <- store_seed

  return(res)
}



anes_data <- function(dv = "vote", missingness = .1) {
  #Load complete data
  store_seed <- .Random.seed

  anes <- anes2008_complete_cases

  miss.data <-
    cbind(1, scale(model.matrix(runif(nrow(
      anes
    )) ~ anes[, 2] + anes[, 3] + anes[, 9] + anes[, 11]))[,-1])
  coef <-
    matrix(rnorm(ncol(anes[,-c(2, 3, 9, 11)]) * ncol(miss.data)), ncol(miss.data), ncol(anes[,-c(2, 3, 9, 11)]))
  miss.eta <- miss.data %*% coef
  miss.error <-
    mi::rdata.frame(
      nrow(anes),
      restrictions = "none",
      n_full = ncol(anes[,-c(2, 3, 9, 11)]),
      n_partial = 0
    )$true
  miss.eta <- miss.eta + .3 * miss.error
  miss.pr <-
    apply(miss.eta, 2, plogis) - matrix(runif(nrow(miss.eta) * ncol(miss.eta)), nrow(miss.eta), ncol(miss.eta))
  miss.indic <-
    apply(
      miss.pr,
      2,
      FUN = function(x) {
        x >= quantile(x, (1 - missingness))
      }
    )
  miss.indic <-
    cbind(
      miss.indic[, 1],
      rep(FALSE, nrow(miss.indic)),
      rep(FALSE, nrow(miss.indic)),
      miss.indic[, 2:6],
      rep(FALSE, nrow(miss.indic)),
      miss.indic[, 7],
      rep(FALSE, nrow(miss.indic))
    )
  anes.miss <- anes
  is.na(anes.miss) <- miss.indic

  if (dv == "vote") {
    analysis_model <- function(x) {
      nnet::multinom(
        vote ~ age + female + as.numeric(education) + married + white + as.numeric(income) + religion,
        data = x
      )
    }
    true_values <- coef(analysis_model(anes))
  }

  if (dv == "time") {
    analysis_model <- function(x) {
      lm(
        time ~ age + female + as.numeric(education) + married + white + as.numeric(income) + religion,
        data = x
      )
    }
    true_values <- coef(analysis_model(anes))
  }

  if (dv == "imp_enviro") {
    analysis_model <- function(x) {
      glm(
        imp_enviro ~ age + female + as.numeric(education) + married + white + as.numeric(income) + religion,
        family = binomial(link = "logit"),
        data = x
      )
    }
    true_values <- coef(analysis_model(anes))
  }
  if (dv == "jobs_r") {
    analysis_model <- function(x) {
      MASS::polr(
        jobs_r ~ age + female + as.numeric(education) + married + white + as.numeric(income) + religion,
        data = x
      )
    }
    true_values <- coef(analysis_model(anes))
  }


  res <- list(
    D = anes,
    D_mis = anes.miss,
    analysis_model = analysis_model,
    true_values = true_values,
    dgp_name = paste0("anes_data_", dv, "_", gsub("\\.", "", paste(missingness))),
    missingness_patterns = c(0.1)
  )
  attr(res, "seed") <- store_seed

  return(res)
}

tbm_data <-
  function(n = 500,
           missingness = "mcar1",
           complexity = "additive") {
    ## Function to create datasets:


    data.gen  <- function(DGP, n, w.rand.error = TRUE, mean.vec) {
      X1 <- rgamma(n, 8, 2)
      X2 <- rgamma(n, 10, 1)
      X3.5 <-  MASS::mvrnorm(n, c(2, 3, 6), diag(c(1.5, 0.5, 3.3)))
      X6.8 <-
        t(rmultinom(n, size = 1, prob = c(1 / 3, 1 / 3, 1 / 3)))
      X9.10 <-
        MASS::mvrnorm(n, mu = c(-0.3, 2), Sigma = matrix(c(1.5, 0.685, 0.685, 5.5), 2))
      X11.40 <- MASS::mvrnorm(n, mu = mean.vec, Sigma = diag(30))

      ### Create outcomes
      Y  <- switch(
        DGP
        ,
        "additive" = X1  +  X2 + X3.5[, 1] + X9.10[, 2]
        ,
        "multiplicative" =  X1 * X2 * X3.5[, 1] * X9.10[, 2]
        ,
        "add+multip" = 1.5 + (3.2) * X1 + 0.1 * X2 - (2.8) * X9.10[, 2] + 1.2 *
          X1 * X2 -
          1.2 * X1 * X9.10[, 2] - 0.2 * X2 * X9.10[, 2] + 2.4 * X1 *
          X2 * X9.10[, 2]
        ,
        "complicated" =  ifelse(
          X9.10[, 2] < 2.5,
          c(
            X1 - X1 ^ 2 - X2 ^ 2  - 15 * X1 * X2 * X9.10[, 2] + poly(X9.10[, 2], 3, raw =
                                                                       TRUE) %*% c(10, -5, 0.9)
          ),
          1750 + 350 * X9.10[, 2]
        )
      )
      if (w.rand.error) {
        Y  <- Y + sd(Y) * rnorm(n)
      }

      ### Create dataset
      temp.data  <-  data.frame(
        Y = Y
        ,
        X1 = X1
        ,
        X2 = X2
        ,
        X3 = X3.5[, 1]
        ,
        X4 = X3.5[, 2]
        ,
        X5 = X3.5[, 3]
        ,
        X6 = X6.8[, 1]
        ,
        X7 = X6.8[, 2]
        ,
        X8 = X6.8[, 3]
        ,
        X9 = X9.10[, 1]
        ,
        X10 = X9.10[, 2]
      )
      temp.data <- cbind(temp.data, as.data.frame(X11.40))
      return(as.matrix(temp.data))
    }

    store_seed <- .Random.seed
    #Common mean vector for irrelevant variables
    irr_mean_vec <- sample(2:10, 30, replace = TRUE)
    D <- data.gen(DGP = complexity,
                  n = n,
                  mean.vec = irr_mean_vec)
    ## 100 train datasets

    if (missingness == "complete") {
      D_mis <- D
    } else if (missingness == "mcar1") {
      U_M <- array(runif(prod(dim(D))), dim(D))
      M <- U_M > 0.19
      M[, 1] <- TRUE
      D_mis <- D
      D_mis[!M] <- NA



    }

    res <- list(
      D = D,
      D_mis = D_mis,
      analysis_model = function(x) {
        lm(x[, 1] ~ x[, 2] + x[, 3] + x[, 4] + x[, 11])
      },
      true_values = c(0, 1, 1, 1, 1),
      dgp_name = paste0("tbm_data_", missingness),
      missingness_patterns = c("mcar1")
    )
    attr(res, "seed") <- store_seed


    return(res)

  }




mixed_data <-
  function(n = 1000,
           missingness = "mcar1",
           coefs = c(0, 1, 0, -2)) {
    store_seed <- .Random.seed
    x <- rnorm(n)
    bin <- rbinom(n, 1, 0.5)

    y <-
      coefs[1] + coefs[2] * x + coefs[3] * bin + coefs[4] * x * bin + rnorm(n, 0, 0.2)

    bimod <- c(rnorm(n / 2,-4), rnorm(n / 2, 4))


    D <- cbind(x, bin, y)



    if (missingness == "complete") {
      D_mis <- D
    } else if (missingness == "mcar1") {
      U_M <- array(runif(prod(dim(D))), dim(D))
      M <- U_M > 0.19
      M[, 3] <- T
      D_mis <- D
      D_mis[!M] <- NA



    }
    res <- list(
      D = D,
      D_mis = D_mis,
      analysis_model = function(x) {
        lm(x[, 3] ~ x[, 1] * x[, 2])
      },
      true_values = coefs,
      dgp_name = paste0("mixed_data_", missingness),
      missingness_patterns = c("mcar1")
    )
    attr(res, "seed") <- store_seed


    return(res)

  }


marbach_data <-
  function(n = 1000,
           missingness = "mar1") {
    store_seed <- .Random.seed
    coefs <- c(0, 0, 0, 1)
    x <- runif(n,-5, 5)
    z <- rbinom(n, 1, 0.5)

    y <-
      coefs[1] + coefs[2] * x + coefs[3] * z + coefs[4] * x * z + rnorm(n, 0, 1)


    D <- cbind(y, x, z)



    if (missingness == "complete") {
      D_mis <- D
    } else if (missingness == "mar1") {
      missing_p <- runif(2, 0.1, 0.5)

      p <- ifelse(z == 1, max(missing_p), min(missing_p))
      y_miss <- rbinom(n, 1, p)



      M <- array(TRUE, dim(D))


      M[y_miss == 1, 1] <- FALSE


      D_mis <- D
      D_mis[!M] <- NA

      res <- list(
        D = D,
        D_mis = D_mis,
        analysis_model = function(x) {
          lm(x[, 1] ~ x[, 2] * x[, 3])
        },
        true_values = coefs,
        dgp_name = paste0("marbach_data_", missingness),
        missingness_patterns = c("mar1")
      )
      attr(res, "seed") <- store_seed
      return(res)

    }


  }



hd_data <- function(n = 500, missingness = "mar1") {
  store_seed <- .Random.seed
  mus <- rep(0, 5)

  varcov <- matrix(0.8, nrow = 5, ncol = 5)

  diag(varcov) <- 1

  D_latent <- MASS::mvrnorm(n, mus, varcov)

  D <- (D_latent > 0) * 1


  if (missingness == "complete") {
    D_mis <- D
  } else if (missingness == "mar1") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    U_M[rowSums(D[, 3:5]) > 0, ] <- 1
    U_M[, 3:5] <- 1
    M <- U_M > 0.2
    D_mis <- D
    D_mis[!M] <- NA



  } else if (missingness == "mar2") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    U_M[rowSums(D[, 3:5]) > 0, ] <- 1
    U_M[, 3:5] <- 1
    M <- U_M > 0.5
    D_mis <- D
    D_mis[!M] <- NA



  } else if (missingness == "mar3") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    U_M[rowSums(D[, 3:5]) > 0, ] <- 1
    U_M[, 3:5] <- 1
    M <- U_M > 0.8
    D_mis <- D
    D_mis[!M] <- NA



  }


  res <- list(
    D = D,
    D_mis = D_mis,
    analysis_model = function(x) {
      glm(x[, 1] ~ x[, 3] + x[, 4], family = binomial(link = logit))
    },
    true_values = c(-1.912, 1.912, 1.912),
    dgp_name = paste0("hd_data_", missingness),
    missingness_patterns = c("mar1", "mar2", "mar3")
  )
  attr(res, "seed") <- store_seed


  return(res)
}




amelia_data <- function(n = 500, missingness = "mcar1") {
  # set.seed(seed)
  store_seed <- .Random.seed

  mus <- rep(0, 5)
  sds <- rep(1, 5)

  sds_mat <- diag(sds)

  cor_mat <- matrix(
    c(
      1 ,
      -.12,
      -.1,
      .5,
      .1,-.12,
      1,
      .1,
      -.6,
      .1,-.1,
      .1,
      1,
      -.5,
      .1,
      .5,
      -.6,
      -.5,
      1,
      .1,
      .1,
      .1,
      .1,
      .1,
      1
    ),
    nrow = 5,
    ncol = 5,
    byrow = T
  )

  varcov <- sds_mat %*% cor_mat %*% sds_mat

  D <- MASS::mvrnorm(n, mus, varcov)

  if (missingness == "complete") {
    D_mis <- D
  } else if (missingness == "mcar1") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    M <- U_M > 0.06
    M[, 4] <- T
    D_mis <- D
    D_mis[!M] <- NA



  } else if (missingness == "mcar2") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    M <- U_M > 0.19
    M[, 4] <- T

    D_mis <- D

    D_mis[!M] <- NA


  } else if (missingness == "mar1") {
    U_M <- array(runif(prod(dim(D))), dim(D))

    M <- U_M > 0.06

    M[, 4] <- T

    M[, 2] <- !(D[, 4] < -1 & U_M[, 2] < 0.9)
    M[, 3] <- !(D[, 4] < -1 & U_M[, 3] < 0.9)

    D_mis <- D

    D_mis[!M] <- NA


  } else if (missingness == "mar2") {
    U_M <- array(runif(prod(dim(D))), dim(D))

    M <- U_M > 0.12

    M[, 4] <- T

    M[, 2] <- !(D[, 4] < -0.4 & U_M[, 2] < 0.9)
    M[, 3] <- !(D[, 4] < -0.4 & U_M[, 3] < 0.9)

    D_mis <- D

    D_mis[!M] <- NA


  } else if (missingness == "ni") {
    U_M <- array(runif(prod(dim(D))), dim(D))

    M <- U_M > 0.06

    M[, 1] <- D[, 1] > -0.95

    M[, 2] <- !(D[, 4] < -0.52)
    M[, 3] <- !(D[, 3] > 0.48)

    D_mis <- D

    D_mis[!M] <- NA


  }

  res <- list(
    D = D,
    D_mis = D_mis,
    analysis_model = function(x) {
      lm(x[, 1] ~ x[, 2] + x[, 3])
    },
    true_values = c(0, -0.11, -0.089),
    dgp_name = paste0("amelia_data_", missingness),
    missingness_patterns = c("mcar1", "mcar2", "mar1", "mar2", "ni")
  )
  attr(res, "seed") <- store_seed


  return(res)
}





summarize_mi_analysis <-
  function(res, true_values) {
    nsim <- length(res)


    tmp <- do.call(abind::abind, c(res, along = 3))

    RB <- rowMeans(tmp[, "estimate", ]) - true_values
    PB <-
      100 * abs((rowMeans(tmp[, "estimate", ]) - true_values) / true_values)
    NB <-
      (rowMeans(tmp[, "estimate", ]) - true_values) / apply(tmp[, "estimate", ], 1, sd)
    CR <-
      rowMeans(tmp[, "2.5 %", ] < replicate(nsim, true_values) &
                 replicate(nsim, true_values) < tmp[, "97.5 %", ])
    AW <- rowMeans(tmp[, "97.5 %", ] - tmp[, "2.5 %", ])
    RMSE <-
      sqrt(rowMeans((tmp[, "estimate", ] - replicate(nsim, true_values)) ^ 2))
    data.frame(RB, PB, NB, CR, AW, RMSE)
  }




lambda_fct <- function(x) {
  (x + 1) / (x + 3)
}

pool_mi <- function(mi_obj, analysis_model) {
  fit_list <- lapply(mi_obj, analysis_model)

  m <- length(fit_list)

  k <- length(coef(fit_list[[1]]))
  df_com <- summary(fit_list[[1]])$df[2]

  Q_bar <- rowMeans(sapply(fit_list, coef))

  U_bar <- array(rowMeans(sapply(fit_list, vcov)), dim = c(k, k))

  B <-
    array(rowSums(sapply(fit_list, function(x)
      (coef(x) - Q_bar) %*% t(coef(x) - Q_bar))) / (m - 1), dim = c(k, k))

  T_var <- U_bar + (1 + 1 / m) * B

  return(
    list(
      pooled_coef = Q_bar,
      within_var = U_bar,
      across_var = B,
      total_var = T_var,
      m = m,
      k = k,
      df_com = df_com
    )
  )
}


ci_mi <- function(pool_obj, congenial = FALSE) {
  if (congenial) {
    gamma_m <-
      (1 + 1 / pool_obj$m) * sum(diag(pool_obj$across_var %*% solve(pool_obj$total_var))) / pool_obj$k

    df_m <- (pool_obj$m - 1) * (gamma_m) ^ -2

    df_tilde_m <-
      pool_obj$df_com * (((
        lambda_fct(pool_obj$df_com) * (1 - gamma_m)
      ) ^ -1 + pool_obj$df_com / df_m) ^ -1)

    ci <- array(pool_obj$pooled_coef + c(
      qt(c(0.025), df = df_tilde_m) * sqrt(diag(pool_obj$total_var)),-qt(c(0.025), df = df_tilde_m) * sqrt(diag(pool_obj$total_var))
    ),
    dim = c(pool_obj$k, 2))
  } else {
    gamma_m <-
      (1 + 1 / pool_obj$m) * sum(diag(pool_obj$across_var %*% solve(2 * pool_obj$total_var))) / pool_obj$k

    df_m <- (pool_obj$m - 1) * (gamma_m) ^ -2

    df_tilde_m <-
      pool_obj$df_com * (((
        lambda_fct(pool_obj$df_com) * (1 - gamma_m)
      ) ^ -1 + pool_obj$df_com / df_m) ^ -1)

    ci <- array(pool_obj$pooled_coef + c(
      qt(c(0.025), df = df_tilde_m) * sqrt(diag(2 * pool_obj$total_var)),-qt(c(0.025), df = df_tilde_m) * sqrt(diag(2 *
                                                                                                                      pool_obj$total_var))
    ),
    dim = c(pool_obj$k, 2))


  }


  res <- cbind(pool_obj$pooled_coef, ci)
  colnames(res) <- c("estimate", "2.5 %", "97.5 %")
  rownames(res) <- 1:pool_obj$k
  return(res)
}


analyze_mi <- function(mi_obj, analysis_model, congenial = FALSE) {
  pool_obj <- pool_mi(mi_obj, analysis_model)

  res <- ci_mi(pool_obj, congenial = congenial)

  return(res)
}


dgp_fct <-
  function(n_data = 500,
           dgp_name = "amelia",
           missingness = "mcar1",
           seed = NULL,
           ...) {
    if (!is.null(seed)) {
      set.seed(seed)
    }

    if (!exists(".Random.seed"))
      set.seed(NULL)


    if (dgp_name == "amelia") {
      data_list <- amelia_data(n_data, missingness)
    } else if (dgp_name == "hd") {
      data_list <- hd_data(n_data, missingness)
    } else if (dgp_name == "mixed") {
      data_list <- mixed_data(n_data, missingness)
      colnames(data_list$D) <- NULL
    } else if (dgp_name == "tbm") {
      data_list <- tbm_data(n_data, missingness, ...)
    } else {
      stop("The data generating process you selected is not yet implemented...")
    }
    names(data_list) <- c("D_full", "D")
    class(data_list) <- "mi_experiment"
    attr(data_list, "n_data") <- n_data
    attr(data_list, "dgp_name") <- dgp_name
    attr(data_list, "missingness") <- missingness
    attr(data_list, "seed") <- seed
    return(data_list)
  }


quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
