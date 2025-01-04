### External Functions

#' Calculate Alpha
#'
#' This function calculates the percentage of values that lie outside of n standard deviations in a normal distribution. This percentage, often called alpha, is used to set the false positive rate in statistical functions.
#'
#' * 3 standard deviations minimizes false positives (1 in 370) at the cost of more missed positives (confirmatory)
#' * 2 standard deviations balances false positives (1 in 22) with missed positives (balanced)
#' * 1.5 standard deviations accepts more false positives (1 in 7) to avoid missed positives (exploratory)
#'
#' @param N Number of standard deviations
#' @return Alpha
#' @export
alpha <- function(N) {
  return( pnorm(-N, 0, 1)*2 )
}

#' Estimate a Dynamic Control Chart
#'
#' This function estimates a Dynamic Control Chart (DCC) using Bayesian statistical methods. The DCC is a Phase II tool for substantiating whether a change resulted in an improved outcome.
#'
#' @param y Outcome data
#' @param x Time series data
#' @param intervention_start Index when change began
#' @param ignored Indices of data points to ignore when estimating the model
#' @param covariates A matrix of covariates for the Outcome data (optional)
#' @param model A DCC model specification
#' @param nsamples Number of posterior samples to simulate
#' @param nburnin Number of posterior samples to discard (i.e., samples until Bayesian convergence)
#' @return An object of class dcc which is a list with the following components: y, x, intervention_start, ignored, p, full_samples
#' @export
dcc <- function(y, x, intervention_start=NULL, ignored=NULL, covariates=NULL, model=normal_model(), nsamples=10000, nburnin=1000) {
  if(class(model) != "dcc_model_specification")
    stop("DCC requires a model specification")

  if(is.null(y) | is.null(x))
    stop("DCC requires y and x data")

  if(length(y) != length(x))
    stop("DCC requires y and x data are the same length")

  if(length(y) < 2)
    stop("DCC requires at least 2 data points")

  if(any(is.na(x)))
    stop("DCC requires no missing x data")

  if(!is.null(intervention_start)) {
    if(any(is.na(y[setdiff(1:(intervention_start-1), ignored)])))
      stop("DCC requires no missing y data before the intervention start")
  } else {
    if(any(is.na(y)))
      stop("DCC requires no missing y data")
  }

  if(!is.null(intervention_start))
    if(intervention_start < 2)
      stop("DCC requires at least 1 data point before the intervention start")

  if(!is.null(intervention_start))
    if(intervention_start > length(x))
      stop("DCC requires at least 1 data point after the intervention start")

  if(!is.null(intervention_start))
    if(!all(x == cummax(x)))
      stop("DCC requires x data increase monotonically")

  if(model$name == "beta" & (any(na.omit(y) < 0) | any(na.omit(y) > 1)))
    stop("DCC Beta Model requires y data is between 0 and 1")

  if(!is.null(covariates))
    if(!is.matrix(covariates))
      stop("DCC requires covariates are in a matrix")

  if(!is.null(covariates))
    if(nrow(covariates) != length(x))
      stop("DCC requires the number of covariate rows match the number of data points")

  if(model$name == "normal")
    full_samples <- dcc_normal(y, x, intervention_start, ignored, model, nsamples, nburnin)
  else if(model$name == "normal_growth")
    full_samples <- dcc_normal_growth(y, x, intervention_start, ignored, model, nsamples, nburnin)
  else if(model$name == "beta")
    full_samples <- dcc_beta(y, x, intervention_start, ignored, model, nsamples, nburnin)
  else
    stop("DCC requires a valid model specification")

  return(
    structure(
      list(
        y = y,
        x = x,
        intervention_start = intervention_start,
        ignored = ignored,
        p = calculate_p_value(y, intervention_start, ignored, full_samples),
        full_samples = full_samples
      ),
      class = "dcc"
    )
  )
}

#' Estimate an Exchangeability Chart
#'
#' This function estimates an Exchangeability Chart (E-Chart) using randomization/permutation methods. The E-Chart is a Phase I tool for examining the exchangeability of outcomes across subgroups.
#'
#' @param y Outcome data
#' @param x Subgroup data
#' @param ignored Subgroups to ignore when estimating the model
#' @param stat_fn Statistical function for aggregating the outcome data by subgroup
#' @param nsamples Number of randomization/permutation samples to simulate
#' @return An object of class echart which is a list with the following components: y, x, stat_fn, ignored, full_samples
#' @export
echart <- function(y, x, stat_fn, ignored=NULL, nsamples=10000) {
  if(is.null(y) | is.null(x))
    stop("E-Chart requires y and x data")

  if(length(y) != length(x))
    stop("E-Chart requires y and x data are the same length")

    if(any(is.na(x)))
    stop("E-Chart requires no missing x data")

  if(any(is.na(y)))
    stop("E-Chart requires no missing y data")

  if(any(!is.numeric(y)))
    stop("E-Chart requires numeric y data")

  if(length(unique(x)) < 2)
    stop("E-Chart requires at least 2 subgroups in x data")

  full_samples <- data.frame(
    x=rep(x, times=nsamples),
    y=rep(y, times=nsamples),
    y_perm=sample(y[which(!(x %in% ignored))], length(x)*nsamples, replace=T),
    permutation_iter=rep(1:nsamples, each=length(x))
  )

  return(
    structure(
      list(
        y = y,
        x = x,
        stat_fn = stat_fn,
        ignored = ignored,
        full_samples = full_samples
      ),
      class = "echart"
    )
  )
}

#' Estimate an Exchangeability Chart from Subgroup Percentages
#'
#' This function estimates an Exchangeability Chart (E-Chart) from subgroup percentages and sizes. This is a helper function which generates binary outcome data based on subgroup percentages and sizes and then passes the data to the echart function.
#'
#' @param y Percentage data
#' @param x Subgroup data
#' @param n Size data
#' @param ignored Subgroups to ignore when estimating the model
#' @param nsamples Number of randomization/permutation samples to simulate
#' @return An object of class echart which is a list with the following components: y, x, stat_fn, ignored, full_samples
#' @export
echart_from_percentages <- function(y, x, n, stat_fn=base::mean, ignored=NULL, nsamples=10000) {
  df <- data.frame(
    y = c( rep( 1, times = sum(round(y*n)) ), rep( 0, times = sum(round((1-y)*n)) ) ),
    x = c( rep( x, times = round(y*n) ), rep( x, times = round((1-y)*n) ) )
  )

  return( echart(df$y, df$x, stat_fn, ignored, nsamples) )
}

#' Depreciated function. Use echart function instead.
#'
#' @export
ucc <- function(y, x, stat_fn, nsamples=10000) {
  stop("Depreciated function. Use echart function instead.")
}

#' Normal DCC Model Specification
#'
#' This function specifies a stable normal data model for a DCC.
#'
#' @param alpha_prior_mean The prior mean of the normal data model
#' @param alpha_prior_sd The confidence of the prior mean
#' @param sigma_prior_mode The prior standard deviation of the normal data model
#' @param sigma_prior_concentration The confidence of the prior standard deviation (sample size)
#' @return A Normal DCC model specification with specified hyperparameters
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
#' @param alpha_prior_mean The prior intercept of the normal growth data model
#' @param alpha_prior_sd The confidence of the prior intercept
#' @param beta_prior_mean The prior slope of the normal growth data model
#' @param beta_prior_sd The confidence of the prior slope
#' @param sigma_prior_mode The prior standard deviation of the normal growth data model
#' @param sigma_prior_concentration The confidence of the prior standard deviation (sample size)
#' @return A Normal Growth DCC model specification with specified hyperparameters
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
#' @param mu_prior_mean The prior mean of the beta data model
#' @param mu_prior_concentration The confidence of the prior mean (sample size)
#' @param phi_prior_mean The prior concentration of the beta data model
#' @param phi_prior_concentration The confidence of the prior concentration (sample size)
#' @param variance_prior_mean The prior variance of the beta data model (alternative to phi_prior_mean)
#' @return A Beta DCC model specification with specified hyperparameters
#' @export
beta_model <- function(mu_prior_mean=NULL, mu_prior_concentration=NULL, phi_prior_mean=NULL, phi_prior_concentration=NULL, variance_prior_mean=NULL) {
  if(is.null(phi_prior_mean) & !is.null(variance_prior_mean)) {
    phi_prior_mean_val <- ((mu_prior_mean*(1-mu_prior_mean))/variance_prior_mean)-1

    if(phi_prior_mean_val < 0)
      stop("Beta DCC Model requires phi greater than zero.")
  } else {
    phi_prior_mean_val <- phi_prior_mean
  }

  return(
    structure(
      list(
        mu_prior_mu = mu_prior_mean,
        mu_prior_phi = mu_prior_concentration,
        phi_prior_mu = phi_prior_mean_val,
        phi_prior_phi = phi_prior_concentration,
        name = "beta"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Plot a Dynamic Control Chart
#'
#' This function plots a Dynamic Control Chart (DCC).
#'
#' @param dcc An estimated DCC object
#' @param alpha False positive rate (Default: 2 standard deviations)
#' @export
plot.dcc <- function(dcc, alpha=0.04550026) {
  if(!is.null(dcc$intervention_start)) {
    df <- dcc$full_samples
    df <- dplyr::group_by(df, model_iter)
    df <- dplyr::arrange(df, model_iter, x)
    df <- dplyr::mutate(df,
      i = 1:dplyr::n(),
      ignored = as.character(as.integer(i %in% dcc$ignored)),
      y_cum = cumsum(ifelse(i >= dcc$intervention_start, y, 0)),
      y_model_cum = cumsum(ifelse(i >= dcc$intervention_start, y_model, 0))
    )
    df <- dplyr::ungroup(df)
    df <- dplyr::mutate(df,
      y_effect = y - y_model,
      y_effect_cum = y_cum - y_model_cum
    )
    df <- dplyr::group_by(df, x, i, ignored)
    df <- dplyr::summarize(df,
      y_mean = mean(y),
      y_model_mean = mean(y_model),
      y_model_ll = quantile(y_model, probs=alpha/2),
      y_model_ul = quantile(y_model, probs=1-alpha/2),
      y_effect_mean = mean(y_effect, na.rm=T),
      y_effect_ll = quantile(y_effect, probs=alpha/2, na.rm=T),
      y_effect_ul = quantile(y_effect, probs=1-alpha/2, na.rm=T),
      y_effect_cum_mean = mean(y_effect_cum, na.rm=T),
      y_effect_cum_ll = quantile(y_effect_cum, probs=alpha/2, na.rm=T),
      y_effect_cum_ul = quantile(y_effect_cum, probs=1-alpha/2, na.rm=T),
      .groups="drop"
    )

    p1 <- ggplot2::ggplot(df) +
      ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_model_ll, ymax = y_model_ul), fill="black", color=NA, alpha=0.2) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = mean(df[(dcc$intervention_start-1):dcc$intervention_start,]$x)), color="black", linetype="dashed") +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y_model_mean), color="black", alpha=0.6) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y_mean, fill = ignored), shape=21, color="#1f77b4", size=1.5) +
      ggplot2::scale_fill_manual(values=c("#1f77b4", "#ffffff")) +
      ggplot2::labs(x = "x", y = "y") +
      ggplot2::theme_minimal() +
      ggplot2::guides(fill="none")

    p2 <- ggplot2::ggplot(df) +
      ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_effect_ll, ymax = y_effect_ul), fill="#1f77b4", color=NA, alpha=0.2) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = mean(df[(dcc$intervention_start-1):dcc$intervention_start,]$x)), color="black", linetype="dashed") +
      ggplot2::geom_hline(ggplot2::aes(yintercept = 0), color="black", linetype="dashed") +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y_effect_mean), color="#1f77b4", alpha=0.6) +
      ggplot2::labs(x = "x", y = "y effect") +
      ggplot2::theme_minimal()

    p3 <- ggplot2::ggplot(df) +
      ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_effect_cum_ll, ymax = y_effect_cum_ul), fill="#1f77b4", color=NA, alpha=0.2) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = df[dcc$intervention_start-1,]$x), color="black", linetype="dashed") +
      ggplot2::geom_hline(ggplot2::aes(yintercept = 0), color="black", linetype="dashed") +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y_effect_cum_mean), color="#1f77b4", alpha=0.6) +
      ggplot2::labs(x = "x", y = "y effect cumsum") +
      ggplot2::theme_minimal()

    plot_design <- "AA
                    BC"
    patchwork::wrap_plots(A=p1, B=p2, C=p3, design=plot_design)
  } else {
    df1 <- dcc$full_samples
    df1 <- dplyr::group_by(df1, model_iter)
    df1 <- dplyr::arrange(df1, model_iter, x)
    df1 <- dplyr::mutate(df1,
      i = 1:dplyr::n(),
      ignored = as.character(as.integer(i %in% dcc$ignored))
    )
    df1 <- dplyr::ungroup(df1)

    df1 <- dplyr::group_by(df1, x, i, ignored)
    df1 <- dplyr::summarize(df1,
      y_mean = mean(y),
      y_model_mean = mean(y_model),
      y_model_ll = quantile(y_model, probs=alpha/2),
      y_model_ul = quantile(y_model, probs=1-alpha/2),
      .groups="drop"
    )

    p1 <- ggplot2::ggplot(df1) +
      ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_model_ll, ymax = y_model_ul), fill="black", color=NA, alpha=0.2) +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y_model_mean), color="black", alpha=0.6) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y_mean, fill = ignored), shape=21, color="#1f77b4", size=1.5) +
      ggplot2::scale_fill_manual(values=c("#1f77b4", "#ffffff")) +
      ggplot2::labs(x = "x", y = "y") +
      ggplot2::theme_minimal() +
      ggplot2::guides(fill="none")

    df2 <- dcc$full_samples

    df2 <- dplyr::group_by(df2, model_iter)
    df2 <- dplyr::arrange(df2, model_iter, x)
    df2 <- dplyr::mutate(df2, i = 1:dplyr::n())
    df2 <- dplyr::filter(df2, !(i %in% dcc$ignored))
    df2 <- dplyr::ungroup(df2)

    p2 <- ggplot2::ggplot() +
      ggplot2::geom_density(ggplot2::aes(x = df2$y_model), linetype="dashed") +
      ggplot2::geom_density(ggplot2::aes(x = unique(df2$y)), fill="#1f77b4", color=NA, alpha=0.4) +
      ggplot2::labs(x = "y") +
      ggplot2::theme_minimal()

    plot_design <- "A
                    B"
    patchwork::wrap_plots(A=p1, B=p2, design=plot_design)
  }
}

#' Plot an Exchangeability Chart
#'
#' This function plots an Exchangeability Chart (E-Chart). The plot can be faceted across strata by specifying stratification data.
#'
#' @param echart An estimated E-Chart object
#' @param stratify_by Stratification data
#' @param alpha False positive rate (Default: 3 standard deviations)
#' @param funnel Sort subgroups by subgroup size (N)
#' @export
plot.echart <- function(echart, stratify_by=NULL, funnel=F, alpha=0.002699796) {
  if(!is.null(stratify_by)) {
    nsamples <- max(echart$full_samples$permutation_iter)

    df <- echart$full_samples
    df$stratify_by <- rep(stratify_by, times=nsamples)

    df <- dplyr::group_by(df, permutation_iter, stratify_by, x)
    df <- dplyr::summarize(df,
      n = dplyr::n(),
      stat = do.call(echart$stat_fn, list(y)),
      stat_perm = do.call(echart$stat_fn, list(y_perm)),
      .groups="drop"
    )

    df <- dplyr::group_by(df, stratify_by, x)
    df <- dplyr::summarize(df,
      n = dplyr::first(n),
      stat = dplyr::first(stat),
      stat_ll = as.vector(quantile(stat_perm, p=alpha/2)),
      stat_ul = as.vector(quantile(stat_perm, p=1-alpha/2)),
      .groups="drop"
    )
    df <- dplyr::mutate(df,
      assignable_cause = stat > stat_ul | stat < stat_ll,
      ignored = as.character(as.integer(x %in% echart$ignored))
    )

    if(funnel) {
      df <- dplyr::mutate(df, x_val = n)
    } else {
      df <- dplyr::group_by(df, stratify_by)
      df <- dplyr::mutate(df, x_val = 1:dplyr::n())
      df <- dplyr::ungroup(df)
    }

    ggplot2::ggplot(df) +
      (if(funnel) {
        ggplot2::geom_ribbon(ggplot2::aes(x = x_val, ymin = stat_ll, ymax = stat_ul), data=function(df) { return( rbind( df, dplyr::mutate(dplyr::summarize(dplyr::group_by(dplyr::arrange(df, x_val), stratify_by), dplyr::across(dplyr::everything(), dplyr::first), .groups = "drop"), x_val = x_val - 0.5), dplyr::mutate(dplyr::summarize(dplyr::group_by(dplyr::arrange(df, x_val), stratify_by), dplyr::across(dplyr::everything(), dplyr::last), .groups = "drop"), x_val = x_val + 0.5) ) ) }, alpha=0.2)
      } else {
        pammtools::geom_stepribbon(ggplot2::aes(x = x_val, ymin = stat_ll, ymax = stat_ul), data=function(df) { return( dplyr::mutate(rbind( df, dplyr::mutate(dplyr::summarize(dplyr::group_by(dplyr::arrange(df, x_val), stratify_by), dplyr::across(dplyr::everything(), dplyr::last), .groups = "drop"), x_val = x_val + 1) ), x_val = x_val - 0.5) ) }, alpha=0.2)
      }) +
      ggplot2::geom_point(ggplot2::aes(x = x_val, y = stat, color=assignable_cause, shape=ignored)) +
      ggplot2::geom_text(ggplot2::aes(x = x_val, y = stat, label=x, color=assignable_cause), size=2, nudge_x = (diff(range(df$x_val))+1)*0.015*length(unique(df$stratify_by)), hjust="left") +
      ggplot2::scale_color_manual(values=c("black", "red")) +
      ggplot2::scale_shape_manual(values=c(19, 1)) +
      ggplot2:::scale_x_continuous(breaks=seq(min(df$x_val), max(df$x_val), by=ceiling(diff(range(df$x_val))/8))) +
      ggplot2::labs(y = "stat", x = (if(funnel) { "N" } else { "index" })) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.ticks.x = ggplot2::element_line(), panel.grid.minor.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank()) +
      ggplot2::guides(shape="none", color = "none") +
      ggplot2::facet_grid(cols=ggplot2::vars(stratify_by))
  } else {
    df <- echart$full_samples
    df <- dplyr::group_by(df, permutation_iter, x)
    df <- dplyr::summarize(df,
      n = dplyr::n(),
      stat = do.call(echart$stat_fn, list(y)),
      stat_perm = do.call(echart$stat_fn, list(y_perm)),
      .groups="drop"
    )

    df <- dplyr::group_by(df, x)
    df <- dplyr::summarize(df,
      n = dplyr::first(n),
      stat = dplyr::first(stat),
      stat_ll = as.vector(quantile(stat_perm, p=alpha/2)),
      stat_ul = as.vector(quantile(stat_perm, p=1-alpha/2)),
      .groups="drop"
    )
    df <- dplyr::mutate(df,
      assignable_cause = stat > stat_ul | stat < stat_ll,
      ignored = as.character(as.integer(x %in% echart$ignored))
    )

    if(funnel) {
      df <- dplyr::mutate(df, x_val = n)
    } else {
      df <- dplyr::mutate(df, x_val = 1:dplyr::n())
    }

    ggplot2::ggplot(df) +
      (if(funnel) {
        ggplot2::geom_ribbon(ggplot2::aes(x = x_val, ymin = stat_ll, ymax = stat_ul), data=function(df) { return( rbind( df, dplyr::mutate(dplyr::first(dplyr::arrange(df, x_val)), x_val = x_val - 0.5), dplyr::mutate(dplyr::last(dplyr::arrange(df, x_val)), x_val = x_val + 0.5) ) ) }, alpha=0.2)
      } else {
        pammtools::geom_stepribbon(ggplot2::aes(x = x_val, ymin = stat_ll, ymax = stat_ul), data=function(df) { return( dplyr::mutate(rbind( df, dplyr::mutate(dplyr::last(dplyr::arrange(df, x_val)), x_val = x_val + 1) ), x_val = x_val - 0.5) ) }, alpha=0.2)
      }) +
      ggplot2::geom_point(ggplot2::aes(x = x_val, y = stat, color=assignable_cause, shape=ignored)) +
      ggplot2::geom_text(ggplot2::aes(x = x_val, y = stat, label=x, color=assignable_cause), size=2, nudge_x = (diff(range(df$x_val))+1)*0.015, hjust="left") +
      ggplot2::scale_color_manual(values=c("black", "red")) +
      ggplot2::scale_shape_manual(values=c(19, 1)) +
      ggplot2:::scale_x_continuous(breaks=seq(min(df$x_val), max(df$x_val), by=ceiling(diff(range(df$x_val))/8))) +
      ggplot2::labs(y = "stat", x = (if(funnel) { "N" } else { "index" })) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.ticks.x = ggplot2::element_line(), panel.grid.minor.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank()) +
      ggplot2::guides(shape="none", color = "none")
  }
}

### Internal Functions

dcc_normal <- function(y, x, intervention_start, ignored, model, nsamples, nburnin) {
  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- 1:length(y)
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

  parm_names <- LaplacesDemon::as.parm.names(list(alpha=0, sigma=0))

  ld_data <- list(
    y = y[pre_interval],
    x = x[pre_interval],
    mon.names = "LP",
    parm.names = parm_names,
    pos.alpha = grep("^alpha", parm_names),
    pos.sigma = grep("^sigma", parm_names)
  )

  ld_model <- function(parm, d) {
    alpha <- parm[d$pos.alpha]

    sigma <- LaplacesDemon::interval(parm[d$pos.sigma], 0.0001, Inf)
    parm[d$pos.sigma] <- sigma

    # log-prior
    if(is.numeric(model$alpha_prior_mu) & is.numeric(model$alpha_prior_sigma))
      alpha_prior <- dnorm(alpha, model$alpha_prior_mu, model$alpha_prior_sigma, log=T)
    else {
      alpha_prior <- dnorm(alpha, mean(y, na.rm=T), sd(y, na.rm=T), log=T)
    }

    if(is.numeric(model$sigma_prior_alpha) & is.numeric(model$sigma_prior_beta))
      sigma_prior <- LaplacesDemon::dinvgamma(sigma, model$sigma_prior_alpha, model$sigma_prior_beta, log=T)
    else {
      sigma_prior <- LaplacesDemon::dhalfcauchy(sigma, sd(y, na.rm=T)*1.1)
    }

    # log-likelihood
    yhat <- alpha + 0*d$x

    LL <- sum(dnorm(d$y, yhat, sigma, log=T))

    # log-posterior
    LP <- LL + alpha_prior + sigma_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(yhat), yhat, sigma), parm=parm))
  }

  ld_initial_values <- c(
    mean(y[pre_interval]),
    sd(y[pre_interval])
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
    y_model_intercept = rep(ld_posterior$alpha, each=length(y)),
    y_model_noise = rnorm(nrow(ld_posterior)*length(y), 0, rep(ld_posterior$sigma, each=length(y))),
    model_iter = rep(1:nrow(ld_posterior), each=length(y))
  )
  full_samples$y_model <- full_samples$y_model_intercept + full_samples$y_model_noise

  return(full_samples)
}

dcc_normal_growth <- function(y, x, intervention_start, ignored, model, nsamples, nburnin) {
  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- 1:length(y)
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

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

    sigma <- LaplacesDemon::interval(parm[d$pos.sigma], 0.0001, Inf)
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

    LL <- sum(dnorm(d$y, yhat, sigma, log=T))

    # log-posterior
    LP <- LL + alpha_prior + beta_prior + sigma_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(yhat), yhat, sigma), parm=parm))
  }

  ld_initial_values <- c(
    y[pre_interval][1],
    (y[pre_interval][length(y[pre_interval])]-y[pre_interval][1])/(x[pre_interval][length(x[pre_interval])]-x[pre_interval][1]),
    sd(y[pre_interval])
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
    y_model_intercept = rep(ld_posterior$alpha, each=length(y)),
    y_model_slope = rep(ld_posterior$beta, each=length(y))*rep(x, times=nrow(ld_posterior)),
    y_model_noise = rnorm(nrow(ld_posterior)*length(y), 0, rep(ld_posterior$sigma, each=length(y))),
    model_iter = rep(1:nrow(ld_posterior), each=length(y))
  )
  full_samples$y_model <- full_samples$y_model_intercept + full_samples$y_model_slope + full_samples$y_model_noise

  return(full_samples)
}

dcc_beta <- function(y, x, intervention_start, ignored, model, nsamples, nburnin) {
  y <- dplyr::case_when(y < 0.0001 ~ 0.0001, y > 0.9999 ~ 0.9999, TRUE ~ y)

  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- 1:length(y)
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

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
    mu <- LaplacesDemon::interval(parm[d$pos.mu], 0.0001, 0.9999)
    parm[d$pos.mu] <- mu

    phi <- LaplacesDemon::interval(parm[d$pos.phi], 0.0001, Inf)
    parm[d$pos.phi] <- phi

    # log-prior
    if(is.numeric(model$mu_prior_mu) & is.numeric(model$mu_prior_phi)) {
      mu_prior <- dbeta(mu, model$mu_prior_mu*model$mu_prior_phi, (1 - model$mu_prior_mu)*model$mu_prior_phi, log=T)
    } else {
      mu_prior <- 0
    }

    if(is.numeric(model$phi_prior_mu) & is.numeric(model$phi_prior_phi)){
      phi_prior <- dgamma(phi, shape=model$phi_prior_mu*model$phi_prior_phi, scale=1/model$phi_prior_phi, log=T)
    } else {
      phi_prior <- 0
    }

    # log-likelihood
    alpha <- mu*phi
    beta <- (1 - mu)*phi
    LL <- sum(dbeta(d$y, alpha, beta, log=T))

    # log-posterior
    LP <- LL + mu_prior + phi_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rbeta(length(d$x), alpha, beta), parm=parm))
  }

  ld_initial_values <- c(
    mean(y[pre_interval]),
    sd(y[pre_interval])
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
    y_model_mean = rep(ld_posterior$mu, each=length(y)),
    y_model_concentration = rep(ld_posterior$phi, each=length(y)),
    y_model = rbeta(nrow(ld_posterior)*length(y), rep(ld_posterior$mu, each=length(y))*rep(ld_posterior$phi, each=length(y)), (1 - rep(ld_posterior$mu, each=length(y)))*rep(ld_posterior$phi, each=length(y))),
    model_iter = rep(1:nrow(ld_posterior), each=length(y))
  )

  return(full_samples)
}

dcc_beta_control <- function(y, x, intervention_start, ignored, covariates, model, nsamples, nburnin) {
  y <- dplyr::case_when(y < 0.0001 ~ 0.0001, y > 0.9999 ~ 0.9999, TRUE ~ y)
  nmatches <- 4

  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- 1:length(y)
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

  parm_names <- LaplacesDemon::as.parm.names(list(alpha=0, beta=rep(0, M), gamma=rep(0, M), phi=0))

  ld_data <- list(
    y = y[pre_interval],
    x = x[pre_interval],
    s = covariates[pre_interval,],
    mon.names = "LP",
    parm.names = parm_names,
    pos.alpha = grep("alpha", parm_names),
    pos.beta = grep("beta", parm_names),
    pos.gamma = grep("gamma", parm_names),
    pos.phi = grep("phi", parm_names)
  )

  ld_model <- function(parm, d) {
    alpha <- parm[d$pos.alpha]

    beta <- LaplacesDemon::interval(parm[d$pos.beta], 0.0001, Inf)
    parm[d$pos.beta] <- beta

    gamma <- LaplacesDemon::interval(parm[d$pos.gamma], 0.0001, 0.9999)
    parm[d$pos.gamma] <- gamma

    beta_sigma <- rep(1.96, M)
    beta_sigma[round(gamma) == 0] <- 0.1

    phi <- interval(parm[d$pos.phi], 0.0001, Inf)
    parm[d$pos.phi] <- phi

    # log-prior
    alpha_prior <- dnorm(alpha, qnorm(mean(y, na.rm=T)), 1.96, log=T)
    beta_prior <- sum(dnorm(beta, 0, beta_sigma, log=T))
    gamma_prior <- sum(dbinom(round(gamma), 1, nmatches/M, log=T))
    phi_prior <- dgamma(phi, shape=5/20, scale=20, log=T)

    # log-likelihood
    betas <- as.vector(tcrossprod(d$s, t(round(gamma)*beta)))
    mu <- LaplacesDemon::interval(pnorm(alpha + betas), 0.0001, 0.9999)

    a <- mu*phi
    b <- (1 - mu)*phi

    LL <- sum(dbeta(d$y, a, b, log=T), na.rm=T)

    # log-posterior
    LP <- LL + alpha_prior + beta_prior + gamma_prior + phi_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rbeta(length(d$x), a, b), parm=parm))
  }

  ld_model <- function(parm, d) {
    mu <- LaplacesDemon::interval(parm[d$pos.mu], 0.0001, 0.9999)
    parm[d$pos.mu] <- mu

    phi <- LaplacesDemon::interval(parm[d$pos.phi], 0.0001, Inf)
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
    LL <- sum(dbeta(d$y, alpha, beta, log=T))

    # log-posterior
    LP <- LL + mu_prior + phi_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rbeta(length(d$x), alpha, beta), parm=parm))
  }

  ld_initial_values <- c(
    mean(y[pre_interval]),
    sd(y[pre_interval])
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
    y_model_mean = rep(ld_posterior$mu, each=length(y)),
    y_model_concentration = rep(ld_posterior$phi, each=length(y)),
    y_model = rbeta(nrow(ld_posterior)*length(y), rep(ld_posterior$mu, each=length(y))*rep(ld_posterior$phi, each=length(y)), (1 - rep(ld_posterior$mu, each=length(y)))*rep(ld_posterior$phi, each=length(y))),
    model_iter = rep(1:nrow(ld_posterior), each=length(y))
  )

  return(full_samples)
}

calculate_p_value <- function(y, intervention_start, ignored, full_samples) {
  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- 1:length(y)
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

  post_interval <- setdiff(post_interval, which(is.na(y)))

  if(length(post_interval) == 0)
    return(NA)

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
