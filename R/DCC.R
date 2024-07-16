### External Functions

#' Estimate a Dynamic Control Chart
#'
#' This function estimates a model from baseline data and simulates counterfactual data for future time points.
#'
#' @param y Response data
#' @param x Time data
#' @param intervention_start The index in the data where the intervention began.
#' @param model A DCC model specification
#' @param nsamples The number of counterfactual samples to simulate.
#' @param nburnin The number of samples to discard (i.e., samples until Bayesian convergence).
#' @return An object of class dcc which is a list with the following components: y, x, intervention_start, p, full_samples.
#' @export
dcc <- function(y, x, intervention_start, model=normal_model(), nsamples=1000, nburnin=1000) {
  if(class(model) != "dcc_model_specification")
    stop("DCC requires a model specification")

  if(is.null(y) | is.null(x))
    stop("DCC requires y and x data")

  if(length(y) != length(x))
    stop("DCC requires y and x data are the same length")

  if(length(y) < 2)
    stop("DCC requires at least 2 data points")

  if(intervention_start < 2)
    stop("DCC requires at least 1 data point before the intervention start")

  if(intervention_start > length(y))
    stop("DCC requires at least 1 data point after the intervention start")

  if(!all(x == cummin(x)))
    stop("DCC requires x data increase monotonically")

  if(model$name == "beta" & (any(y < 0) | any(y > 1)))
    stop("DCC Beta Model requires y data is between 0 and 1")

  if(model$name == "poisson" & (any(y < 0) | !all.equal(y, as.integer(y))))
    stop("DCC Poisson Model requires y data is zero or a positive integer")

  if(model$name == "normal")
    full_samples <- dcc_normal(y, x, intervention_start, model, nsamples, nburnin)
  else if(model$name == "normal_growth")
    full_samples <- dcc_normal_growth(y, x, intervention_start, model, nsamples, nburnin)
  else if(model$name == "beta")
    full_samples <- dcc_beta(y, x, intervention_start, model, nsamples, nburnin)
  else if(model$name == "poisson")
    full_samples <- dcc_poisson(y, x, intervention_start, model, nsamples, nburnin)
  else
    stop("DCC requires a model specification")

  return(
    structure(
      list(
        y = y,
        x = x,
        intervention_start = intervention_start,
        p = calculate_p_value(y, intervention_start, full_samples),
        full_samples = full_samples
      ),
      class = "dcc"
    )
  )
}

#' Normal DCC Model Specification
#'
#' This function specifies a stable normal data model for a DCC.
#'
#' @param alpha_prior_mean The prior mean of the normal data model.
#' @param alpha_prior_sd The standard deviation as confidence of the prior mean.
#' @param sigma_prior_mode The prior standard deviation of the normal data model.
#' @param sigma_prior_concentration The confidence of the prior standard deviation.
#' @return A Normal DCC model specification with specified hyperparameters.
#' @export
normal_model <- function(alpha_prior_mean=NULL, alpha_prior_sd=NULL, sigma_prior_mode=NULL, sigma_prior_concentration=NULL) {
  return(
    structure(
      list(
        alpha_prior_mu = alpha_prior_mean,
        alpha_prior_sigma = alpha_prior_sd,
        sigma_prior_alpha = sigma_prior_concentration,
        sigma_prior_beta = sigma_prior_mode*(sigma_prior_concentration + 1),
        name = "normal"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Normal Growth DCC Model Specification
#'
#' This function specifies a normal growth data model for a DCC.
#'
#' @param alpha_prior_mean The prior intercept of the normal growth data model.
#' @param alpha_prior_sd The standard deviation as confidence of the prior intercept.
#' @param beta_prior_mean The prior slope of the normal growth data model.
#' @param beta_prior_sd The standard deviation as confidence of the prior slope.
#' @param sigma_prior_mode The prior standard deviation of the normal growth data model.
#' @param sigma_prior_concentration The confidence of the prior standard deviation.
#' @return A Normal Growth DCC model specification with specified hyperparameters.
#' @export
normal_growth_model <- function(alpha_prior_mean=NULL, alpha_prior_sd=NULL, beta_prior_mean=NULL, beta_prior_sd=NULL, sigma_prior_mode=NULL, sigma_prior_concentration=NULL) {
  return(
    structure(
      list(
        alpha_prior_mu = alpha_prior_mean,
        alpha_prior_sigma = alpha_prior_sd,
        beta_prior_mu = beta_prior_mean,
        beta_prior_sigma = beta_prior_sd,
        sigma_prior_alpha = sigma_prior_concentration,
        sigma_prior_beta = sigma_prior_mode*(sigma_prior_concentration + 1),
        name = "normal_growth"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Beta DCC Model Specification
#'
#' This function specifies a stable beta data model for a DCC.
#'
#' @param mu_prior_mean The prior mean of the beta data model.
#' @param mu_prior_concentration The confidence of the prior mean.
#' @param phi_prior_mean The prior concentration of the beta data model.
#' @param phi_prior_concentration The confidence of the prior concentration.
#' @return A Beta DCC model specification with specified hyperparameters.
#' @export
beta_model <- function(mu_prior_mean=NULL, mu_prior_concentration=NULL, phi_prior_mean=NULL, phi_prior_concentration=NULL) {
  return(
    structure(
      list(
        mu_prior_mu = mu_prior_mean,
        mu_prior_phi = mu_prior_concentration,
        phi_prior_mu = phi_prior_mean,
        phi_prior_phi = phi_prior_concentration,
        name = "beta"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Poisson DCC Model Specification
#'
#' This function specifies a stable poisson data model for a DCC.
#'
#' @param lambda_prior_mean The prior mean of the poisson data model.
#' @param lambda_prior_concentration The confidence of the prior mean.
#' @return A Poisson DCC model specification with specified hyperparameters.
#' @export
poisson_model <- function(lambda_prior_mean=NULL, lambda_prior_concentration=NULL) {
  return(
    structure(
      list(
        lambda_prior_mu = lambda_prior_mean,
        lambda_prior_phi = lambda_prior_concentration,
        name = "poisson"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Plot a Dynamic Control Chart
#'
#' This function plots a Dynamic Control Chart.
#'
#' @param dcc An estimated Dynamic Control Chart object.
#' @export
plot.dcc <- function(dcc) {
  df <- dcc$full_samples
  df <- dplyr::group_by(df, model_iter)
  df <- dplyr::arrange(df, model_iter, x)
  df <- dplyr::mutate(df,
    i = 1:dplyr::n(),
    z = ifelse(i >= dcc$intervention_start, 1, 0),
    y_cum = cumsum(ifelse(z == 1, y, 0)),
    y_model_cum = cumsum(ifelse(z == 1, y_model, 0))
  )
  df <- dplyr::ungroup(df)
  df <- dplyr::mutate(df,
    y_effect = y - y_model,
    y_effect_cum = y_cum - y_model_cum
  )
  df <- dplyr::group_by(df, x, i, z)
  df <- dplyr::summarize(df,
    y_mean = mean(y),
    y_model_mean = mean(y_model),
    y_model_ll = quantile(y_model, probs=0.10),
    y_model_ul = quantile(y_model, probs=0.90),
    y_effect_mean = mean(y_effect, na.rm=T),
    y_effect_ll = quantile(y_effect, probs=0.10, na.rm=T),
    y_effect_ul = quantile(y_effect, probs=0.90, na.rm=T),
    y_effect_cum_mean = mean(y_effect_cum, na.rm=T),
    y_effect_cum_ll = quantile(y_effect_cum, probs=0.10, na.rm=T),
    y_effect_cum_ul = quantile(y_effect_cum, probs=0.90, na.rm=T),
    .groups="drop"
  )

  p1 <- ggplot2::ggplot(df) +
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_model_ll, ymax = y_model_ul), fill="black", color=NA, alpha=0.2) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(df[(dcc$intervention_start-1):dcc$intervention_start,]$x)), color="black", linetype="dashed") +
    ggplot2::geom_line(ggplot2::aes(x = x, y = y_model_mean), color="black", alpha=0.6) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y_mean), color="#1f77b4", size=1.5) +
    ggplot2::labs(x = "t", y = "y") +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(df) +
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_effect_ll, ymax = y_effect_ul), fill="#1f77b4", color=NA, alpha=0.2) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(df[(dcc$intervention_start-1):dcc$intervention_start,]$x)), color="black", linetype="dashed") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), color="black", linetype="dashed") +
    ggplot2::geom_line(ggplot2::aes(x = x, y = y_effect_mean), color="#1f77b4", alpha=0.6) +
    ggplot2::labs(x = "t", y = "y effect") +
    ggplot2::theme_minimal()

  p3 <- ggplot2::ggplot(df) +
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_effect_cum_ll, ymax = y_effect_cum_ul), fill="#1f77b4", color=NA, alpha=0.2) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = df[dcc$intervention_start-1,]$x), color="black", linetype="dashed") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), color="black", linetype="dashed") +
    ggplot2::geom_line(ggplot2::aes(x = x, y = y_effect_cum_mean), color="#1f77b4", alpha=0.6) +
    ggplot2::labs(x = "t", y = "y effect cumsum") +
    ggplot2::theme_minimal()

  plot_design <- "AA
                  BC"
  patchwork::wrap_plots(A=p1, B=p2, C=p3, design=plot_design)
}

### Internal Functions

dcc_normal <- function(y, x, intervention_start, model, nsamples, nburnin) {
  pre_interval <- 1:(intervention_start-1)
  post_interval <- intervention_start:length(y)

  parm_names <- LaplacesDemon::as.parm.names(list(alpha=0, sigma=0))

  ld_data <- list(
    y = y[pre_interval],
    x = x[pre_interval],
    mon.names = "LP",
    parm.names = parm_names,
    pos.alpha = grep("alpha", parm_names),
    pos.sigma = grep("sigma", parm_names)
  )

  ld_model <- function(parm, d) {
    alpha <- parm[d$pos.alpha]

    sigma <- LaplacesDemon::interval(parm[d$pos.sigma], 0.001, Inf)
    parm[d$pos.sigma] <- sigma

    # log-prior
    if(is.numeric(model$alpha_prior_mu) & is.numeric(model$alpha_prior_sigma))
      alpha_prior <- dnorm(alpha, model$alpha_prior_mu, model$alpha_prior_sigma, log=T)
    else {
      alpha_prior <- 0
    }

    if(is.numeric(model$sigma_prior_alpha) & is.numeric(model$sigma_prior_beta))
      sigma_prior <- LaplacesDemon::dinvgamma(sigma, model$sigma_prior_alpha, model$sigma_prior_beta, log=T)
    else {
      sigma_prior <- 0
    }

    # log-likelihood
    yhat <- alpha + 0*d$x

    LL <- sum(dnorm(d$y, yhat, sigma, log=T), na.rm=T)

    # log-posterior
    LP <- LL + alpha_prior + sigma_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(yhat), yhat, sigma), parm=parm))
  }

  ld_initial_values <- c(
    mean(y),
    sd(y)
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior)),
    z = rep(c(rep.int(0, intervention_start-1), rep.int(1, length(y) - intervention_start + 1)), times=nrow(ld_posterior)),
    y_model_intercept = rep(ld_posterior$alpha, each=length(y)),
    y_model_noise = rnorm(nrow(ld_posterior)*length(y), 0, rep(ld_posterior$sigma, each=length(y))),
    model_iter = rep(1:nrow(ld_posterior), each=length(y))
  )
  full_samples$y_model <- full_samples$y_model_intercept + full_samples$y_model_noise

  return(full_samples)
}

dcc_normal_growth <- function(y, x, intervention_start, model, nsamples, nburnin) {
  pre_interval <- 1:(intervention_start-1)
  post_interval <- intervention_start:length(y)

  parm_names <- LaplacesDemon::as.parm.names(list(alpha=0, beta=0, sigma=0))

  ld_data <- list(
    y = y[pre_interval],
    x = x[pre_interval],
    mon.names = "LP",
    parm.names = parm_names,
    pos.alpha = grep("alpha", parm_names),
    pos.beta = grep("beta", parm_names),
    pos.sigma = grep("sigma", parm_names)
  )

  ld_model <- function(parm, d) {
    alpha <- parm[d$pos.alpha]

    beta <- parm[d$pos.beta]

    sigma <- LaplacesDemon::interval(parm[d$pos.sigma], 0.001, Inf)
    parm[d$pos.sigma] <- sigma

    # log-prior
    if(is.numeric(model$alpha_prior_mu) & is.numeric(model$alpha_prior_sigma))
      alpha_prior <- dnorm(alpha, model$alpha_prior_mu, model$alpha_prior_sigma, log=T)
    else {
      alpha_prior <- 0
    }

    if(is.numeric(model$beta_prior_mu) & is.numeric(model$beta_prior_sigma))
      beta_prior <- dnorm(beta, model$beta_prior_mu, model$beta_prior_sigma, log=T)
    else {
      beta_prior <- 0
    }

    if(is.numeric(model$sigma_prior_alpha) & is.numeric(model$sigma_prior_beta))
      sigma_prior <- LaplacesDemon::dinvgamma(sigma, model$sigma_prior_alpha, model$sigma_prior_beta, log=T)
    else {
      sigma_prior <- 0
    }

    # log-likelihood
    yhat <- alpha + beta*d$x

    LL <- sum(dnorm(d$y, yhat, sigma, log=T), na.rm=T)

    # log-posterior
    LP <- LL + alpha_prior + beta_prior + sigma_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(yhat), yhat, sigma), parm=parm))
  }

  ld_initial_values <- c(
    y[1],
    (df$y[length(df$y)]-y[1])/(df$x[length(df$x)]-x[1]),
    sd(y)
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior)),
    z = rep(c(rep.int(0, intervention_start-1), rep.int(1, length(y) - intervention_start + 1)), times=nrow(ld_posterior)),
    y_model_intercept = rep(ld_posterior$alpha, each=length(y)),
    y_model_slope = rep(ld_posterior$beta, each=length(y))*rep(x, times=nrow(ld_posterior)),
    y_model_noise = rnorm(nrow(ld_posterior)*length(y), 0, rep(ld_posterior$sigma, each=length(y))),
    model_iter = rep(1:nrow(ld_posterior), each=length(y))
  )
  full_samples$y_model <- full_samples$y_model_intercept + full_samples$y_model_slope + full_samples$y_model_noise

  return(full_samples)
}

dcc_beta <- function(y, x, intervention_start, model, nsamples, nburnin) {
  y <- dplyr::case_when(y <= 0.0001 ~ 0.0001, y >= 0.9999 ~ 0.9999, TRUE ~ y)

  pre_interval <- 1:(intervention_start-1)
  post_interval <- intervention_start:length(y)

  parm_names <- LaplacesDemon::as.parm.names(list(mu=0, phi=0))

  ld_data <- list(
    y = y[pre_interval],
    x = x[pre_interval],
    mon.names = "LP",
    parm.names = parm_names,
    pos.mu = grep("mu", parm_names),
    pos.phi = grep("phi", parm_names)
  )

  ld_model <- function(parm, d) {
    mu <- LaplacesDemon::interval(parm[d$pos.mu], 0.001, 0.999)
    parm[d$pos.mu] <- mu

    phi <- LaplacesDemon::interval(parm[d$pos.phi], 0.001, Inf)
    parm[d$pos.phi] <- phi

    # log-prior
    if(is.numeric(model$mu_prior_mu) & is.numeric(model$mu_prior_phi))
      mu_prior <- dbeta(mu, model$mu_prior_mu*model$mu_prior_phi, (1 - model$mu_prior_mu)*model$mu_prior_phi, log=T)
    else {
      mu_prior <- 0
    }

    if(is.numeric(model$phi_prior_mu) & is.numeric(model$phi_prior_phi))
      phi_prior <- dgamma(phi, shape=model$phi_prior_mu*model$phi_prior_phi, scale=1/model$phi_prior_phi, log=T)
    else {
      phi_prior <- 0
    }

    # log-likelihood
    alpha <- mu*phi
    beta <- (1 - mu)*phi
    LL <- sum(dbeta(d$y, alpha, beta, log=T), na.rm=T)

    # log-posterior
    LP <- LL + mu_prior + phi_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rbeta(length(d$x), alpha, beta), parm=parm))
  }

  ld_initial_values <- c(
    mean(y),
    sd(y)
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior)),
    z = rep(c(rep.int(0, intervention_start-1), rep.int(1, length(y) - intervention_start + 1)), times=nrow(ld_posterior)),
    y_model_mean = rep(ld_posterior$mu, each=length(y)),
    y_model_concentration = rep(ld_posterior$phi, each=length(y)),
    y_model = rbeta(nrow(ld_posterior)*length(y), rep(ld_posterior$mu, each=length(y))*rep(ld_posterior$phi, each=length(y)), (1 - rep(ld_posterior$mu, each=length(y)))*rep(ld_posterior$phi, each=length(y))),
    model_iter = rep(1:nrow(ld_posterior), each=length(y))
  )

  return(full_samples)
}

dcc_poisson <- function(y, x, intervention_start, model, nsamples, nburnin) {
  pre_interval <- 1:(intervention_start-1)
  post_interval <- intervention_start:length(y)

  parm_names <- LaplacesDemon::as.parm.names(list(lambda=0))

  ld_data <- list(
    y = y[pre_interval],
    x = x[pre_interval],
    mon.names = "LP",
    parm.names = parm_names,
    pos.lambda = grep("lambda", parm_names)
  )

  ld_model <- function(parm, d) {
    lambda <- LaplacesDemon::interval(parm[d$pos.lambda], 0, Inf)
    parm[d$pos.lambda] <- lambda

    # log-prior
    if(is.numeric(model$lambda_prior_mu) & is.numeric(model$lambda_prior_phi))
      lambda_prior <- dgamma(lambda, shape=model$lambda_prior_mu*model$lambda_prior_phi, scale=1/model$lambda_prior_phi, log=T)
    else {
      lambda_prior <- 0
    }

    # log-likelihood
    LL <- sum(dpois(d$y, lambda, log=T), na.rm=T)

    # log-posterior
    LP <- LL + lambda_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rpois(length(d$x), lambda), parm=parm))
  }

  ld_initial_values <- c(
    mean(y)
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  colnames(ld_posterior) <- c("lambda")
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior)),
    z = rep(c(rep.int(0, intervention_start-1), rep.int(1, length(y) - intervention_start + 1)), times=nrow(ld_posterior)),
    y_model_mean = rep(ld_posterior$lambda, each=length(y)),
    y_model = rpois(nrow(ld_posterior)*length(y), rep(ld_posterior$lambda, each=length(y))),
    model_iter = rep(1:nrow(ld_posterior), each=length(y))
  )

  return(full_samples)
}


calculate_p_value <- function(y, intervention_start, full_samples) {
  pre_interval = 1:(intervention_start-1)
  post_interval = intervention_start:length(y)

  y_sum <- sum(y[post_interval])

  if(length(post_interval) == 1)
    y_model_sums <- matrix(full_samples$y_model, length(y))[post_interval,]
  else
    y_model_sums <- colSums(matrix(full_samples$y_model, length(y))[post_interval,])

  p_position <- min(sum(c(y_model_sums, y_sum) >= y_sum),
                    sum(c(y_model_sums, y_sum) <= y_sum))
  p <- p_position / (length(y_model_sums) + 1)

  return(p)
}
