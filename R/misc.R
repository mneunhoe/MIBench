#' @importFrom magrittr "%>%"
# Function to generate data according to King et al 2001

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
                                                                       TRUE) %*% c(10,-5, 0.9)
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


    #Common mean vector for irrelevant variables
    irr_mean_vec <- sample(2:10, 30, replace = TRUE)
    D <- data.gen(DGP = complexity, n = n, mean.vec = irr_mean_vec)
    ## 100 train datasets

    if (missingness == "complete") {
      return(list(D, D))
    } else if (missingness == "mcar1") {
      U_M <- array(runif(prod(dim(D))), dim(D))
      M <- U_M > 0.19
      M[, 1] <- T
      D_mis <- D
      D_mis[!M] <- NA

      return(list(D, D_mis))

    }

  }





mixed_data <-
  function(n = 1000,
           missingness = "mcar1",
           coefs = c(0, 1, 0,-2)) {
    x <- rnorm(n)
    bin <- rbinom(n, 1, 0.5)

    y <-
      coefs[1] + coefs[2] * x + coefs[3] * bin + coefs[4] * x * bin + rnorm(n, 0, 0.2)

    bimod <- c(rnorm(n / 2, -4), rnorm(n / 2, 4))


    D <- cbind(x, bin, y)



    if (missingness == "complete") {
      return(list(D, D))
    } else if (missingness == "mcar1") {
      U_M <- array(runif(prod(dim(D))), dim(D))
      M <- U_M > 0.19
      M[, 3] <- T
      D_mis <- D
      D_mis[!M] <- NA

      return(list(D, D_mis))

    }


  }


hd_data <- function(n = 500, missingness = "mar1") {
  mus <- rep(0, 5)

  varcov <- matrix(0.8, nrow = 5, ncol = 5)

  diag(varcov) <- 1

  D_latent <- MASS::mvrnorm(n, mus, varcov)

  D <- (D_latent > 0) * 1


  if (missingness == "complete") {
    return(list(D, D))
  } else if (missingness == "mar1") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    U_M[rowSums(D[, 3:5]) > 0,] <- 1
    U_M[, 3:5] <- 1
    M <- U_M > 0.2
    D_mis <- D
    D_mis[!M] <- NA

    return(list(D, D_mis))

  } else if (missingness == "mar2") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    U_M[rowSums(D[, 3:5]) > 0,] <- 1
    U_M[, 3:5] <- 1
    M <- U_M > 0.5
    D_mis <- D
    D_mis[!M] <- NA

    return(list(D, D_mis))

  } else if (missingness == "mar3") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    U_M[rowSums(D[, 3:5]) > 0,] <- 1
    U_M[, 3:5] <- 1
    M <- U_M > 0.8
    D_mis <- D
    D_mis[!M] <- NA

    return(list(D, D_mis))

  }


}



amelia_data <- function(n = 500, missingness = "mcar1") {
  # set.seed(seed)

  mus <- rep(0, 5)
  sds <- rep(1, 5)

  sds_mat <- diag(sds)

  cor_mat <- matrix(
    c(
      1 ,-.12,-.1,
      .5,
      .1,
      -.12,
      1,
      .1,-.6,
      .1,
      -.1,
      .1,
      1,-.5,
      .1,
      .5,-.6,-.5,
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
    return(list(D, D))
  } else if (missingness == "mcar1") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    M <- U_M > 0.06
    M[, 4] <- T
    D_mis <- D
    D_mis[!M] <- NA

    return(list(D, D_mis))

  } else if (missingness == "mcar2") {
    U_M <- array(runif(prod(dim(D))), dim(D))
    M <- U_M > 0.19
    M[, 4] <- T

    D_mis <- D

    D_mis[!M] <- NA

    return(list(D, D_mis))
  } else if (missingness == "mar1") {
    U_M <- array(runif(prod(dim(D))), dim(D))

    M <- U_M > 0.06

    M[, 4] <- T

    M[, 2] <- !(D[, 4] < -1 & U_M[, 2] < 0.9)
    M[, 3] <- !(D[, 4] < -1 & U_M[, 3] < 0.9)

    D_mis <- D

    D_mis[!M] <- NA

    return(list(D, D_mis))
  } else if (missingness == "mar2") {
    U_M <- array(runif(prod(dim(D))), dim(D))

    M <- U_M > 0.12

    M[, 4] <- T

    M[, 2] <- !(D[, 4] < -0.4 & U_M[, 2] < 0.9)
    M[, 3] <- !(D[, 4] < -0.4 & U_M[, 3] < 0.9)

    D_mis <- D

    D_mis[!M] <- NA

    return(list(D, D_mis))
  } else if (missingness == "ni") {
    U_M <- array(runif(prod(dim(D))), dim(D))

    M <- U_M > 0.06

    M[, 1] <- D[, 1] > -0.95

    M[, 2] <- !(D[, 4] < -0.52)
    M[, 3] <- !(D[, 3] > 0.48)

    D_mis <- D

    D_mis[!M] <- NA

    return(list(D, D_mis))
  }


}

summarize_mi_analysis <-
  function(res, true_values) {
    nsim <- length(res)


    tmp <- do.call(abind::abind, c(res, along = 3))

    RB <- rowMeans(tmp[, "estimate",]) - true_values
    PB <-
      100 * abs((rowMeans(tmp[, "estimate",]) - true_values) / true_values)
    NB <-
      (rowMeans(tmp[, "estimate",]) - true_values) / apply(tmp[, "estimate",], 1, sd)
    CR <-
      rowMeans(tmp[, "2.5 %",] < replicate(nsim, true_values) &
                 replicate(nsim, true_values) < tmp[, "97.5 %",])
    AW <- rowMeans(tmp[, "97.5 %",] - tmp[, "2.5 %",])
    RMSE <-
      sqrt(rowMeans((tmp[, "estimate",] - replicate(nsim, true_values)) ^ 2))
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

  if(congenial){
    gamma_m <- (1 + 1 / pool_obj$m) * sum(diag(pool_obj$across_var %*% solve(pool_obj$total_var))) / pool_obj$k

    df_m <- (pool_obj$m - 1) * (gamma_m) ^ -2

    df_tilde_m <-
      pool_obj$df_com * (((lambda_fct(pool_obj$df_com) * (1 - gamma_m)) ^ -1 + pool_obj$df_com / df_m) ^ -1)

    ci <- array(pool_obj$pooled_coef + c(
      qt(c(0.025), df = df_tilde_m) * sqrt(diag(pool_obj$total_var)),
      -qt(c(0.025), df = df_tilde_m) * sqrt(diag(pool_obj$total_var))
    ), dim = c(pool_obj$k, 2))
  } else {

    gamma_m <- (1 + 1 / pool_obj$m) * sum(diag(pool_obj$across_var %*% solve(2*pool_obj$total_var))) / pool_obj$k

    df_m <- (pool_obj$m - 1) * (gamma_m) ^ -2

    df_tilde_m <-
      pool_obj$df_com * (((lambda_fct(pool_obj$df_com) * (1 - gamma_m)) ^ -1 + pool_obj$df_com / df_m) ^ -1)

    ci <- array(pool_obj$pooled_coef + c(
      qt(c(0.025), df = df_tilde_m) * sqrt(diag(2*pool_obj$total_var)),
      -qt(c(0.025), df = df_tilde_m) * sqrt(diag(2*pool_obj$total_var))
    ), dim = c(pool_obj$k, 2))


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
      data_list <- tbm_data(n_data, missingness,...)
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
