### External Functions

#' Calculate Alpha Coefficient
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
alpha_coef <- function(N) {
  return( pnorm(-N, 0, 1)*2 )
}

#' Estimate a Dynamic Control Chart
#'
#' This function estimates a Dynamic Control Chart (DCC) using Bayesian statistical methods. The DCC is a tool for substantiating whether a change resulted in an improved outcome. The DCC can also be used for forecasting an outcome during goal setting.
#'
#' @param y Outcome data
#' @param x Time series data
#' @param intervention_start Index when change began
#' @param ignored Indices of data points to ignore when estimating the model
#' @param covariates A matrix of contemporaneous covariates for the outcome data (optional)
#' @param model A DCC model specification
#' @param nsamples Number of posterior samples to simulate
#' @param nburnin Number of posterior samples to discard (i.e., samples until Bayesian convergence)
#' @return An object of class dcc which is a list with the following components: y, x, intervention_start, ignored, covariates, p, full_samples
#' @export
dcc <- function(y, x, intervention_start=NULL, ignored=NULL, covariates=NULL, model=normal_model(), nsamples=10000, nburnin=1000, param_diagnostic=T) {
  if(class(model) != "dcc_model_specification")
    stop("DCC requires a model specification")

  if(is.null(y) | is.null(x))
    stop("DCC requires y and x data")

  if(!is.numeric(x) | !is.numeric(y))
    stop("DCC requires numeric y and x data")

  if(length(y) != length(x))
    stop("DCC requires y and x data are the same length")

  if(length(y) < 2)
    stop("DCC requires at least 2 data points")

  if(any(is.na(x)))
    stop("DCC requires no missing x data")

  if(!all(x == cummax(x)))
    stop("DCC requires x data increase monotonically")

  if(!is.null(intervention_start)) {
    if(intervention_start < 2)
      stop("DCC requires at least 1 data point before the intervention start")

    if(intervention_start > length(x))
      stop("DCC requires at least 1 data point after the intervention start")

    if(any(is.na(y[1:(intervention_start-1)])))
      stop("DCC requires no missing y data before the intervention start")
  } else {
    if(any(is.na(y)))
      stop("DCC requires no missing y data")
  }

  if(model$name %in% c("beta", "beta_forecast", "beta_growth", "beta_control") & (any(na.omit(y) < 0) | any(na.omit(y) > 1)))
    stop("DCC Beta Models require y data is between 0 and 1")

  if(!is.null(covariates)) {
    if(!is.matrix(covariates))
      stop("DCC requires covariates are in a matrix")

    if(nrow(covariates) != length(x))
      stop("DCC requires the number of covariate rows match the number of data points")

    if(!is.numeric(covariates))
      stop("DCC requires numeric covariate data")

    if(any(is.na(covariates)))
      stop("DCC requires no missing covariate data")
  }

  if(model$name == "normal")
    full_samples <- dcc_normal(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic)
  else if(model$name == "normal_forecast")
    full_samples <- dcc_normal_forecast(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic)
  else if(model$name == "normal_growth")
    full_samples <- dcc_normal_growth(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic)
  else if(model$name == "normal_control")
    full_samples <- dcc_normal_control(y, x, intervention_start, ignored, covariates, model, nsamples, nburnin, param_diagnostic)
  else if(model$name == "beta")
    full_samples <- dcc_beta(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic)
  else if(model$name == "beta_forecast")
    full_samples <- dcc_beta_forecast(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic)
  else if(model$name == "beta_growth")
    full_samples <- dcc_beta_growth(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic)
  else if(model$name == "beta_control")
    full_samples <- dcc_beta_control(y, x, intervention_start, ignored, covariates, model, nsamples, nburnin, param_diagnostic)
  else
    stop("DCC requires a valid model specification")

  return(
    structure(
      list(
        y = y,
        x = x,
        intervention_start = intervention_start,
        ignored = ignored,
        covariates = covariates,
        p = calculate_p_value(y, intervention_start, ignored, full_samples),
        full_samples = full_samples
      ),
      class = "dcc"
    )
  )
}

#' Estimate an Exchangeability Chart
#'
#' This function estimates an Exchangeability Chart (E-Chart) using randomization/permutation methods. The E-Chart is a tool for examining the exchangeability of outcomes across subgroups.
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

  if(length(unique(x)) < 2)
    stop("E-Chart requires at least 2 subgroups in x data")

  if(any(is.na(y)))
    stop("E-Chart requires no missing y data")

  if(any(!is.numeric(y)))
    stop("E-Chart requires numeric y data")

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
#' This function specifies a stationary normal data model for a DCC.
#'
#' @param alpha_prior The prior mean of the normal data model
#' @param alpha_prior_confidence The confidence of the prior mean (scale: standard deviation)
#' @param sigma_prior The prior standard deviation of the normal data model
#' @param sigma_prior_confidence The confidence of the prior standard deviation (rate: sample size)
#' @return A Normal DCC model specification with specified hyperparameters
#' @export
normal_model <- function(alpha_prior=NULL, alpha_prior_confidence=NULL, sigma_prior=NULL, sigma_prior_confidence=NULL) {
  return(
    structure(
      list(
        alpha_prior_mu = alpha_prior,
        alpha_prior_sigma = alpha_prior_confidence,
        sigma_prior_alpha = sigma_prior_confidence,
        sigma_prior_beta = sigma_prior*(sigma_prior_confidence + 1),
        name = "normal"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Normal Forecast DCC Model Specification
#'
#' This function specifies a non-stationary normal data model (i.e., local level model) for a DCC.
#'
#' @param omega_prior The prior standard deviation of the innovations (i.e., the change in the non-stationary mean from point to point)
#' @param omega_prior_confidence The confidence of the prior standard deviation of the innovations (rate: sample size)
#' @param sigma_prior The prior standard deviation of the normal data model
#' @param sigma_prior_confidence The confidence of the prior standard deviation (rate: sample size)
#' @return A Normal Forecast DCC model specification with specified hyperparameters
#' @export
normal_forecast_model <- function(omega_prior=NULL, omega_prior_confidence=NULL, sigma_prior=NULL, sigma_prior_confidence=NULL) {
  return(
    structure(
      list(
        omega_prior_alpha = omega_prior_confidence,
        omega_prior_beta = omega_prior*(omega_prior_confidence + 1),
        sigma_prior_alpha = sigma_prior_confidence,
        sigma_prior_beta = sigma_prior*(sigma_prior_confidence + 1),
        name = "normal_forecast"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Normal Growth DCC Model Specification
#'
#' This function specifies a normal growth data model for a DCC.
#'
#' @param alpha_prior The prior intercept of the normal growth data model
#' @param alpha_prior_confidence The confidence of the prior intercept (scale: standard deviation)
#' @param beta_prior The prior slope of the normal growth data model
#' @param beta_prior_confidence The confidence of the prior slope (scale: standard deviation)
#' @param sigma_prior The prior standard deviation of the normal growth data model
#' @param sigma_prior_confidence The confidence of the prior standard deviation (rate: sample size)
#' @return A Normal Growth DCC model specification with specified hyperparameters
#' @export
normal_growth_model <- function(alpha_prior=NULL, alpha_prior_confidence=NULL, beta_prior=NULL, beta_prior_confidence=NULL, sigma_prior=NULL, sigma_prior_confidence=NULL) {
  return(
    structure(
      list(
        alpha_prior_mu = alpha_prior,
        alpha_prior_sigma = alpha_prior_confidence,
        beta_prior_mu = beta_prior,
        beta_prior_sigma = beta_prior_confidence,
        sigma_prior_alpha = sigma_prior_confidence,
        sigma_prior_beta = sigma_prior*(sigma_prior_confidence + 1),
        name = "normal_growth"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Normal Control DCC Model Specification
#'
#' This function specifies a normal contemporaneous regression data model with stochastic search variable selection (SSVS) for a DCC.
#'
#' @param inclusion_probabilities The prior inclusion probabilities of the contemporaneous covariates (optional)
#' @param sigma_prior The prior standard deviation of the normal data model
#' @param sigma_prior_confidence The confidence of the prior standard deviation (rate: sample size)
#' @param nmatches Number of expected contemporaneous covariate matches in the stochastic search (alternative to specifying inclusion probabilities)
#' @return A Normal Control DCC model specification with specified hyperparameters
#' @export
normal_control_model <- function(inclusion_probabilities=NULL, sigma_prior=NULL, sigma_prior_confidence=NULL, nmatches=1) {
  return(
    structure(
      list(
        inclusion_probabilities = inclusion_probabilities,
        sigma_prior_alpha = sigma_prior_confidence,
        sigma_prior_beta = sigma_prior*(sigma_prior_confidence + 1),
        nmatches = nmatches,
        name = "normal_control"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Beta DCC Model Specification
#'
#' This function specifies a stationary beta data model for a DCC.
#'
#' @param mu_prior The prior mean of the beta data model
#' @param mu_prior_confidence The confidence of the prior mean (rate: sample size)
#' @param phi_prior The prior concentration (precision) of the beta data model
#' @param phi_prior_confidence The confidence of the prior concentration (rate)
#' @param variance_prior The prior variance of the beta data model (alternative to specifying phi prior)
#' @return A Beta DCC model specification with specified hyperparameters
#' @export
beta_model <- function(mu_prior=NULL, mu_prior_confidence=NULL, phi_prior=NULL, phi_prior_confidence=NULL, variance_prior=NULL) {
  phi_prior_val <- phi_prior

  if(is.null(phi_prior_val) & !is.null(variance_prior)) {
    phi_prior_val <- ((mu_prior*(1-mu_prior))/variance_prior)-1
  }

  phi_prior_val <- LaplacesDemon::interval(phi_prior_val, 0.0001, Inf)

  return(
    structure(
      list(
        mu_prior_mu = mu_prior,
        mu_prior_phi = mu_prior_confidence,
        phi_prior_mu = phi_prior_val,
        phi_prior_phi = phi_prior_confidence,
        name = "beta"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Beta Forecast DCC Model Specification
#'
#' This function specifies a non-stationary beta data model (i.e., local level model) for a DCC.
#'
#' @param phi_prior The prior concentration (precision) of the beta data model
#' @param phi_prior_confidence The confidence of the prior concentration (rate)
#' @return A Beta Forecast DCC model specification with specified hyperparameters
#' @export
beta_forecast_model <- function(phi_prior=NULL, phi_prior_confidence=NULL) {
  phi_prior_val <- LaplacesDemon::interval(phi_prior, 0.0001, Inf)

  return(
    structure(
      list(
        phi_prior_mu = phi_prior_val,
        phi_prior_phi = phi_prior_confidence,
        name = "beta_forecast"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Beta Growth DCC Model Specification
#'
#' This function specifies a beta growth data model for a DCC.
#'
#' @param phi_prior The prior concentration (precision) of the beta data model
#' @param phi_prior_confidence The confidence of the prior concentration (rate)
#' @return A Beta Growth DCC model specification with specified hyperparameters
#' @export
beta_growth_model <- function(phi_prior=NULL, phi_prior_confidence=NULL) {
  phi_prior_val <- LaplacesDemon::interval(phi_prior, 0.0001, Inf)

  return(
    structure(
      list(
        phi_prior_mu = phi_prior_val,
        phi_prior_phi = phi_prior_confidence,
        name = "beta_growth"
      ),
      class = "dcc_model_specification"
    )
  )
}

#' Beta Control DCC Model Specification
#'
#' This function specifies a beta contemporaneous regression data model with stochastic search variable selection (SSVS) for a DCC.
#'
#' @param inclusion_probabilities The prior inclusion probabilities of the contemporaneous covariates (optional)
#' @param phi_prior The prior concentration (precision) of the beta data model
#' @param phi_prior_confidence The confidence of the prior concentration (rate)
#' @param nmatches Number of expected contemporaneous covariate matches in the stochastic search (alternative to specifying inclusion probabilities)
#' @return A Beta Control DCC model specification with specified hyperparameters
#' @export
beta_control_model <- function(inclusion_probabilities=NULL, phi_prior=NULL, phi_prior_confidence=NULL, nmatches=1) {
  return(
    structure(
      list(
        inclusion_probabilities = inclusion_probabilities,
        phi_prior_mu = phi_prior,
        phi_prior_phi = phi_prior_confidence,
        nmatches = nmatches,
        name = "beta_control"
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
#' @param export_data Return a data frame instead of a plot object
#' @export
plot.dcc <- function(dcc, alpha=0.04550026, export_data=F) {
  if(!is.null(dcc$intervention_start)) {
    cache_key <- rlang::hash(list(dcc$full_samples, dcc$ignored, dcc$intervention_start, alpha))

    df <- cache_get(cache_key)
    if(is.null(df)) {
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

      cache_set(cache_key, df)
    }

    if(export_data)
      return(df)

    p1 <- ggplot2::ggplot(df) +
      ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_model_ll, ymax = y_model_ul), fill="black", color=NA, alpha=0.2) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = mean(df[(dcc$intervention_start-1):dcc$intervention_start,]$x)), color="black", linetype="dashed") +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y_model_mean), color="black", alpha=0.6) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y_mean, fill = ignored), data=function(df) { return(dplyr::filter(df, !is.na(y_mean))) }, shape=21, color="#1f77b4", size=1.5) +
      ggplot2::scale_fill_manual(values=c("#1f77b4", "#ffffff")) +
      ggplot2::labs(x = "x", y = "y") +
      ggplot2::theme_minimal() +
      ggplot2::guides(fill="none")

    if(!all(is.na(dcc$y[dcc$intervention_start:length(dcc$y)]))) {
      p2 <- ggplot2::ggplot(df) +
        ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_effect_ll, ymax = y_effect_ul), fill="#1f77b4", color=NA, alpha=0.2) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = mean(df[(dcc$intervention_start-1):dcc$intervention_start,]$x)), color="black", linetype="dashed") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), color="black", linetype="dashed") +
        ggplot2::geom_line(ggplot2::aes(x = x, y = y_effect_mean), data=function(df) { return(dplyr::filter(df, !is.na(y_effect_mean))) }, color="#1f77b4", alpha=0.6) +
        ggplot2::labs(x = "x", y = "y effect") +
        ggplot2::theme_minimal()

      p3 <- ggplot2::ggplot(df) +
        ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_effect_cum_ll, ymax = y_effect_cum_ul), fill="#1f77b4", color=NA, alpha=0.2) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = df[dcc$intervention_start-1,]$x), color="black", linetype="dashed") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), color="black", linetype="dashed") +
        ggplot2::geom_line(ggplot2::aes(x = x, y = y_effect_cum_mean), data=function(df) { return(dplyr::filter(df, !is.na(y_effect_cum_mean))) }, color="#1f77b4", alpha=0.6) +
        ggplot2::labs(x = "x", y = "y effect cumsum") +
        ggplot2::theme_minimal()
    } else {
      p2 <- ggplot2::ggplot(df) +
        ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_effect_ll, ymax = y_effect_ul), fill="#1f77b4", color=NA, alpha=0.2) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = mean(df[(dcc$intervention_start-1):dcc$intervention_start,]$x)), color="black", linetype="dashed") +
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), color="black", linetype="dashed") +
        ggplot2::geom_line(ggplot2::aes(x = x, y = y_effect_mean), data=function(df) { return(dplyr::filter(df, !is.na(y_effect_mean))) }, color="#1f77b4", alpha=0.6) +
        ggplot2::labs(x = "x", y = "y effect") +
        ggplot2::coord_cartesian(xlim=c( df[1,]$x, df[dcc$intervention_start,]$x )) +
        ggplot2::theme_minimal()

      p3 <- ggplot2::ggplot(df) +
        ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_model_ll, ymax = y_model_ul), fill="black", color=NA, alpha=0.2) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = mean(df[(dcc$intervention_start-1):dcc$intervention_start,]$x)), color="black", linetype="dashed") +
        ggplot2::geom_line(ggplot2::aes(x = x, y = y_model_mean), color="black", alpha=0.6) +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y_mean, fill = ignored), data=function(df) { return(dplyr::filter(df, !is.na(y_mean))) }, shape=21, color="#1f77b4", size=1.5) +
        ggplot2::geom_text(ggplot2::aes(x = x, y = y_mean, label = round(y_mean, 2)), data=df[dcc$intervention_start-1,], size=2.5, nudge_x = (diff(range(df[(dcc$intervention_start-1):length(df$x),]$x))+1)*0.03, hjust="left", vjust="top", color="#1f77b4") +
        ggplot2::scale_fill_manual(values=c("#1f77b4", "#ffffff")) +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y_model_ll), data=df[dcc$intervention_start,], size=1.5, color="#777") +
        ggplot2::geom_text(ggplot2::aes(x = x, y = y_model_ll, label = round(y_model_ll, 2)), data=df[dcc$intervention_start,], size=2.5, nudge_x = (diff(range(df[(dcc$intervention_start-1):length(df$x),]$x))+1)*0.03, hjust="left", vjust="top") +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y_model_ul), data=df[dcc$intervention_start,], size=1.5, color="#777") +
        ggplot2::geom_text(ggplot2::aes(x = x, y = y_model_ul, label = round(y_model_ul, 2)), data=df[dcc$intervention_start,], size=2.5, nudge_x = (diff(range(df[(dcc$intervention_start-1):length(df$x),]$x))+1)*0.03, hjust="left", vjust="bottom") +
        ggplot2::labs(x = "x", y = "y") +
        ggplot2::coord_cartesian(xlim=c( df[(dcc$intervention_start-1),]$x, df[length(df$x),]$x+(diff(range(df[(dcc$intervention_start-1):length(df$x),]$x))+1)*0.25 )) +
        ggplot2::theme_minimal() +
        ggplot2::guides(fill="none")
    }

    plot_design <- "AA
                    BC"
    patchwork::wrap_plots(A=p1, B=p2, C=p3, design=plot_design)
  } else {
    cache_key <- rlang::hash(list("df1", dcc$full_samples, dcc$ignored, dcc$intervention_start, alpha))

    df1 <- cache_get(cache_key)
    if(is.null(df1)) {
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

      cache_set(cache_key, df1)
    }

    if(export_data)
      return(df1)

    p1 <- ggplot2::ggplot(df1) +
      ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_model_ll, ymax = y_model_ul), fill="black", color=NA, alpha=0.2) +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y_model_mean), color="black", alpha=0.6) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y_mean, fill = ignored), shape=21, color="#1f77b4", size=1.5) +
      ggplot2::scale_fill_manual(values=c("#1f77b4", "#ffffff")) +
      ggplot2::labs(x = "x", y = "y") +
      ggplot2::theme_minimal() +
      ggplot2::guides(fill="none")

    cache_key <- rlang::hash(list("df2", dcc$full_samples, dcc$ignored, dcc$intervention_start, alpha))

    df2 <- cache_get(cache_key)
    if(is.null(df2)) {
      df2 <- dcc$full_samples

      df2 <- dplyr::group_by(df2, model_iter)
      df2 <- dplyr::arrange(df2, model_iter, x)
      df2 <- dplyr::mutate(df2, i = 1:dplyr::n())
      df2 <- dplyr::filter(df2, !(i %in% dcc$ignored))
      df2 <- dplyr::ungroup(df2)

      cache_set(cache_key, df2)
    }

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

#' Plot the Parameters of a Dynamic Control Chart
#'
#' This function plots the posterior distributions for the parameters for a Dynamic Control Chart (DCC).
#'
#' @param dcc An estimated DCC object
#' @param traces Plot the traces instead of the posterior distributions
#' @export
param_plot <- function(dcc, bins=30, traces=F) {
  if(dim(dplyr::select(dcc$full_samples, dplyr::starts_with("param")))[2] == 0)
    stop("No parameter data found. Set param_diagnostic = TRUE in dcc function. ")

  df <- dcc$full_samples
  df <- dplyr::select(df, model_iter, dplyr::starts_with("param"))
  df <- dplyr::group_by(df, model_iter)
  df <- dplyr::summarize(df, dplyr::across(dplyr::everything(), dplyr::first))
  df <- tidyr::pivot_longer(df, -model_iter, names_prefix = "param_", names_to = "param", values_to = "param value")

  if(!traces) {
    ggplot2::ggplot(df) +
      ggplot2::geom_histogram(ggplot2::aes(x = `param value`), bins=bins) +
      ggplot2::theme_minimal() +
      ggplot2::facet_wrap(ggplot2::vars(param), scales="free")
  } else {
    ggplot2::ggplot(df) +
      ggplot2::geom_line(ggplot2::aes(x = model_iter, y = `param value`)) +
      ggplot2::theme_minimal() +
      ggplot2::facet_wrap(ggplot2::vars(param), scales="free_y")
  }
}

#' Plot an Exchangeability Chart
#'
#' This function plots an Exchangeability Chart (E-Chart). The plot can be faceted across strata by specifying stratification data.
#'
#' @param echart An estimated E-Chart object
#' @param stratify_by Stratification data
#' @param funnel Sort subgroups by subgroup size (N)
#' @param alpha False positive rate (Default: 3 standard deviations)
#' @param export_data Return a data frame instead of a plot object
#' @export
plot.echart <- function(echart, stratify_by=NULL, funnel=F, alpha=0.002699796, export_data=F) {
  if(!is.null(stratify_by)) {
    cache_key <- rlang::hash(list(echart$full_samples, echart$stat_fn, echart$ignored, stratify_by, alpha))

    df <- cache_get(cache_key)
    if(is.null(df)) {
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

      cache_set(cache_key, df)
    }

    if(export_data)
      return(df)

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
    cache_key <- rlang::hash(list(echart$full_samples, echart$stat_fn, echart$ignored, stratify_by, alpha))

    df <- cache_get(cache_key)
    if(is.null(df)) {
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

      cache_set(cache_key, df)
    }

    if(export_data)
      return(df)

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

dcc_normal <- function(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic) {
  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- c()
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
      alpha_prior <- 0
    }

    if(is.numeric(model$sigma_prior_alpha) & is.numeric(model$sigma_prior_beta))
      sigma_prior <- LaplacesDemon::dinvgamma(sigma, model$sigma_prior_alpha, model$sigma_prior_beta, log=T)
    else {
      sigma_prior <- 0
    }

    # log-likelihood
    yhat <- alpha

    LL <- sum(dnorm(d$y, yhat, sigma, log=T))

    # log-posterior
    LP <- LL + alpha_prior + sigma_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(yhat), yhat, sigma), parm=parm))
  }

  ld_initial_values <- c(
    mean(y[pre_interval]), # alpha: mean y
    sd(y[pre_interval]) # sigma: sd y
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior))
  )

  if(param_diagnostic) {
    full_samples$param_alpha <- rep(ld_posterior$alpha, each=length(y))
    full_samples$param_sigma <- rep(ld_posterior$sigma, each=length(y))
  }

  full_samples$y_model_intercept <- rep(ld_posterior$alpha, each=length(y))
  full_samples$y_model_noise <- rnorm(nrow(ld_posterior)*length(y), 0, rep(ld_posterior$sigma, each=length(y)))
  full_samples$y_model <- full_samples$y_model_intercept + full_samples$y_model_noise

  full_samples$model_iter <- rep(1:nrow(ld_posterior), each=length(y))

  return(full_samples)
}

dcc_normal_forecast <- function(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic) {
  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- c()
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

  N <- length(x)
  parm_names <- LaplacesDemon::as.parm.names(list(nu=rep(0, N), omega=0, sigma=0))

  y_na_post <- y
  y_na_post[post_interval] <- NA

  ld_data <- list(
    y = y_na_post,
    x = x,
    mon.names = "LP",
    parm.names = parm_names,
    pos.nu = grep("^nu", parm_names),
    pos.omega = grep("^omega", parm_names),
    pos.sigma = grep("^sigma", parm_names)
  )

  ld_model <- function(parm, d) {
    nu <- parm[d$pos.nu]

    omega <- LaplacesDemon::interval(parm[d$pos.omega], 0.0001, Inf)
    parm[d$pos.omega] <- omega

    sigma <- LaplacesDemon::interval(parm[d$pos.sigma], 0.0001, Inf)
    parm[d$pos.sigma] <- sigma

    # log-prior
    nu_prior <- sum(0, dnorm(nu[-1], nu[-N], omega, log=T))

    if(is.numeric(model$omega_prior_alpha) & is.numeric(model$omega_prior_beta))
      omega_prior <- LaplacesDemon::dinvgamma(omega, model$omega_prior_alpha, model$omega_prior_beta, log=T)
    else {
      omega_prior <- 0
    }

    if(is.numeric(model$sigma_prior_alpha) & is.numeric(model$sigma_prior_beta))
      sigma_prior <- LaplacesDemon::dinvgamma(sigma, model$sigma_prior_alpha, model$sigma_prior_beta, log=T)
    else if(length(y[pre_interval]) > 1) {
      sigma_prior <- LaplacesDemon::dinvgamma(sigma, 0.1, sd(diff(y[pre_interval]))*(0.1 + 1), log=T) # rate: 0.1
    } else {
      sigma_prior <- 0
    }

    # log-likelihood
    yhat <- nu

    LL <- sum(dnorm(d$y[pre_interval], yhat[pre_interval], sigma, log=T))

    # log-posterior
    LP <- LL + nu_prior + omega_prior + sigma_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(yhat[pre_interval]), yhat[pre_interval], sigma), parm=parm))
  }

  ld_initial_values <- c(
    rep(mean(y[pre_interval]), N), # nu : mean y
    sd(y[pre_interval]), # omega: sd y
    sd(y[pre_interval]) # sigma: sd y
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior))
  )

  if(param_diagnostic) {
    full_samples <- cbind(full_samples, magrittr::set_colnames(matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("nu"))), each=length(y)), ncol=N, byrow=F), paste0("param_nu_", 1:N)))
    full_samples$param_omega <- rep(ld_posterior$omega, each=length(y))
    full_samples$param_sigma <- rep(ld_posterior$sigma, each=length(y))
  }

  full_samples <- cbind(full_samples, data.frame(y_model_level=c(t(as.matrix(dplyr::select(ld_posterior, starts_with("nu")))))))
  full_samples$y_model_noise <- rnorm(nrow(ld_posterior)*length(y), 0, rep(ld_posterior$sigma, each=length(y)))
  full_samples$y_model <- full_samples$y_model_level + full_samples$y_model_noise

  full_samples$model_iter <- rep(1:nrow(ld_posterior), each=length(y))

  return(full_samples)
}

dcc_normal_growth <- function(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic) {
  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- c()
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
    y[pre_interval][1], # alpha: first y
    0, # beta: 0
    sd(y[pre_interval]) # sigma: sd y
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior))
  )

  if(param_diagnostic) {
    full_samples$param_alpha <- rep(ld_posterior$alpha, each=length(y))
    full_samples$param_beta <- rep(ld_posterior$beta, each=length(y))
    full_samples$param_sigma <- rep(ld_posterior$sigma, each=length(y))
  }

  full_samples$y_model_intercept <- rep(ld_posterior$alpha, each=length(y))
  full_samples$y_model_slope <- rep(ld_posterior$beta, each=length(y))*rep(x, times=nrow(ld_posterior))
  full_samples$y_model_noise <- rnorm(nrow(ld_posterior)*length(y), 0, rep(ld_posterior$sigma, each=length(y)))
  full_samples$y_model <- full_samples$y_model_intercept + full_samples$y_model_slope + full_samples$y_model_noise

  full_samples$model_iter <- rep(1:nrow(ld_posterior), each=length(y))

  return(full_samples)
}

dcc_normal_control <- function(y, x, intervention_start, ignored, covariates, model, nsamples, nburnin, param_diagnostic) {
  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- c()
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

  s <- scale(covariates)
  M <- dim(s)[[2]] # number of columns
  parm_names <- LaplacesDemon::as.parm.names(list(alpha=0, beta=rep(0, M), gamma=rep(0, M), rho=0, sigma=0))

  ld_data <- list(
    y = y[pre_interval],
    x = x[pre_interval],
    s = s[pre_interval,],
    mon.names = "LP",
    parm.names = parm_names,
    pos.alpha = grep("alpha", parm_names),
    pos.beta = grep("beta", parm_names),
    pos.gamma = grep("gamma", parm_names),
    pos.rho = grep("rho", parm_names),
    pos.sigma = grep("sigma", parm_names)
  )

  ld_model <- function(parm, d) {
    alpha <- parm[d$pos.alpha]

    beta <- LaplacesDemon::interval(parm[d$pos.beta], 0, Inf) # no inverse covariates
    parm[d$pos.beta] <- beta

    gamma <- LaplacesDemon::interval(parm[d$pos.gamma], 0.0001, 0.9999)
    parm[d$pos.gamma] <- gamma

    rho <- LaplacesDemon::interval(parm[d$pos.rho], 0.0001, Inf)
    parm[d$pos.rho] <- rho
    rhos <- rep(rho, M)
    rhos[round(gamma) == 0] <- 0.1 # precision = 10

    sigma <- LaplacesDemon::interval(parm[d$pos.sigma], 0.0001, Inf)
    parm[d$pos.sigma] <- sigma

    # log-prior
    alpha_prior <- 0
    beta_prior <- sum(dnorm(beta, 0, rhos, log=T))

    if(!is.null(model$inclusion_probabilities)) {
      gamma_prior <- sum(dbinom(round(gamma), 1, model$inclusion_probabilities, log=T))
    } else {
      gamma_prior <- sum(dbinom(round(gamma), 1, model$nmatches/M, log=T))
    }

    rho_prior <- 0

    if(is.numeric(model$sigma_prior_alpha) & is.numeric(model$sigma_prior_beta))
      sigma_prior <- LaplacesDemon::dinvgamma(sigma, model$sigma_prior_alpha, model$sigma_prior_beta, log=T)
    else {
      sigma_prior <- 0
    }

    # log-likelihood
    betas <- as.vector(tcrossprod(d$s, t(round(gamma)*beta)))
    yhat <- alpha + betas

    LL <- sum(dnorm(d$y, yhat, sigma, log=T), na.rm=T)

    # log-posterior
    LP <- LL + alpha_prior + beta_prior + gamma_prior + rho_prior + sigma_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(length(yhat), yhat, sigma), parm=parm))
  }

  ld_initial_values <- c(
    mean(y[pre_interval]), # alpha: mean y
    rep(0, M), # betas: 0
    rep(1, M), # gammas: 1
    sd(y[pre_interval]), # rho: sd y
    sd(y[pre_interval]) # sigma:  sd y
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior))
  )
  full_samples <- cbind(full_samples, magrittr::set_colnames(matrix(rep(t(s), nrow(ld_posterior)), ncol=M, byrow=T), paste0("s_", 1:M)))

  if(param_diagnostic) {
    full_samples$param_alpha <- rep(ld_posterior$alpha, each=length(y))
    full_samples <- cbind(full_samples, magrittr::set_colnames(matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("beta"))), each=length(y)), ncol=M, byrow=F), paste0("param_beta_", 1:M)))
    full_samples <- cbind(full_samples, magrittr::set_colnames(matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("gamma"))), each=length(y)), ncol=M, byrow=F), paste0("param_gamma_", 1:M)))
    full_samples$param_rho <- rep(ld_posterior$rho, each=length(y))
    full_samples$param_sigma <- rep(ld_posterior$sigma, each=length(y))
  }

  full_samples$y_model_intercept <- rep(ld_posterior$alpha, each=length(y))
  full_samples$y_model_covariates <- rowSums( matrix(rep(t(s), nrow(ld_posterior)), ncol=M, byrow=T) * ( round(matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("gamma"))), each=length(y)), ncol=M, byrow=F)) * matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("beta"))), each=length(y)), ncol=M, byrow=F) ) )
  full_samples$y_model_noise <- rnorm(nrow(ld_posterior)*length(y), 0, rep(ld_posterior$sigma, each=length(y)))
  full_samples$y_model <- full_samples$y_model_intercept + full_samples$y_model_covariates + full_samples$y_model_noise

  full_samples$model_iter <- rep(1:nrow(ld_posterior), each=length(y))

  return(full_samples)
}

dcc_beta <- function(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic) {
  y <- dplyr::case_when(y < 0.0001 ~ 0.0001, y > 0.9999 ~ 0.9999, TRUE ~ y)

  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- c()
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
    a <- mu*phi
    b <- (1 - mu)*phi
    LL <- sum(dbeta(d$y, a, b, log=T))

    # log-posterior
    LP <- LL + mu_prior + phi_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rbeta(length(d$x), a, b), parm=parm))
  }

  ld_initial_values <- c(
    mean(y[pre_interval]), # mu: mean y
    ((mean(y[pre_interval])*(1-mean(y[pre_interval])))/var(y[pre_interval]))-1 # phi: sample size from variance y
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior))
  )

  if(param_diagnostic) {
    full_samples$param_mu <- rep(ld_posterior$mu, each=length(y))
    full_samples$param_phi <- rep(ld_posterior$phi, each=length(y))
  }

  full_samples$y_model_mean <- rep(ld_posterior$mu, each=length(y))
  full_samples$y_model_concentration <- rep(ld_posterior$phi, each=length(y))
  full_samples$y_model <- rbeta(
    nrow(ld_posterior)*length(y),
    full_samples$y_model_mean*full_samples$y_model_concentration,
    (1 - full_samples$y_model_mean)*full_samples$y_model_concentration
  )

  full_samples$model_iter <- rep(1:nrow(ld_posterior), each=length(y))

  return(full_samples)
}

dcc_beta_forecast <- function(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic) {
  y <- dplyr::case_when(y < 0.0001 ~ 0.0001, y > 0.9999 ~ 0.9999, TRUE ~ y)

  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- c()
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

  N <- length(x)
  parm_names <- LaplacesDemon::as.parm.names(list(nu=rep(0, N), omega=0, phi=0))

  y_na_post <- y
  y_na_post[post_interval] <- NA

  ld_data <- list(
    y = y_na_post,
    x = x,
    mon.names = "LP",
    parm.names = parm_names,
    pos.nu = grep("^nu", parm_names),
    pos.omega = grep("^omega", parm_names),
    pos.phi = grep("^phi", parm_names)
  )

  ld_model <- function(parm, d) {
    nu <- parm[d$pos.nu]

    omega <- LaplacesDemon::interval(parm[d$pos.omega], 0.0001, Inf)
    parm[d$pos.omega] <- omega

    phi <- LaplacesDemon::interval(parm[d$pos.phi], 0.0001, Inf)
    parm[d$pos.phi] <- phi

    # log-prior
    nu_prior <- sum(0, dnorm(nu[-1], nu[-N], omega, log=T))

    omega_prior <- 0

    if(is.numeric(model$phi_prior_mu) & is.numeric(model$phi_prior_phi)){
      phi_prior <- dgamma(phi, shape=model$phi_prior_mu*model$phi_prior_phi, scale=1/model$phi_prior_phi, log=T)
    } else if(length(y[pre_interval]) > 1) {
      phi_prior_val <- ((mean(y[pre_interval])*(1-mean(y[pre_interval])))/var(y[pre_interval]))-1
      phi_prior <- dgamma(phi, shape=phi_prior_val*0.05, scale=1/0.05, log=T) # rate: 0.05
    } else {
      phi_prior <- 0
    }

    # log-likelihood
    mu <- LaplacesDemon::interval(pnorm(nu), 0.0001, 0.9999)

    a <- mu*phi
    b <- (1 - mu)*phi

    LL <- sum(dbeta(d$y[pre_interval], a[pre_interval], b[pre_interval], log=T), na.rm=T)

    # log-posterior
    LP <- LL + nu_prior + omega_prior + phi_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rbeta(length(a[pre_interval]), a[pre_interval], b[pre_interval]), parm=parm))
  }

  ld_initial_values <- c(
    rep(qnorm(mean(y[pre_interval])), N), # nu: qnorm mean y
    sd(qnorm(y[pre_interval])), # omega: sd qnorm y
    ((mean(y[pre_interval])*(1-mean(y[pre_interval])))/var(y[pre_interval]))-1 # phi: sample size from variance y
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior))
  )

  if(param_diagnostic) {
    full_samples <- cbind(full_samples, magrittr::set_colnames(matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("nu"))), each=length(y)), ncol=N, byrow=F), paste0("param_nu_", 1:N)))
    full_samples$param_omega <- rep(ld_posterior$omega, each=length(y))
    full_samples$param_phi <- rep(ld_posterior$phi, each=length(y))
  }

  full_samples <- cbind(full_samples, data.frame(y_model_level=c(t(as.matrix(dplyr::select(ld_posterior, starts_with("nu")))))))
  full_samples$y_model_mean <- LaplacesDemon::interval(pnorm(full_samples$y_model_level), 0.0001, 0.9999)
  full_samples$y_model_concentration <- rep(ld_posterior$phi, each=length(y))
  full_samples$y_model <- rbeta(
    nrow(ld_posterior)*length(y),
    full_samples$y_model_mean*full_samples$y_model_concentration,
    (1 - full_samples$y_model_mean)*full_samples$y_model_concentration
  )

  full_samples$model_iter <- rep(1:nrow(ld_posterior), each=length(y))

  return(full_samples)
}

dcc_beta_growth <- function(y, x, intervention_start, ignored, model, nsamples, nburnin, param_diagnostic) {
  y <- dplyr::case_when(y < 0.0001 ~ 0.0001, y > 0.9999 ~ 0.9999, TRUE ~ y)

  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- c()
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

  parm_names <- LaplacesDemon::as.parm.names(list(alpha=0, beta=0, phi=0))

  ld_data <- list(
    y = y[pre_interval],
    x = x[pre_interval],
    mon.names = "LP",
    parm.names = parm_names,
    pos.alpha = grep("alpha", parm_names),
    pos.beta = grep("beta", parm_names),
    pos.phi = grep("phi", parm_names)
  )

  ld_model <- function(parm, d) {
    alpha <- parm[d$pos.alpha]

    beta <- parm[d$pos.beta]

    phi <- LaplacesDemon::interval(parm[d$pos.phi], 0.0001, Inf)
    parm[d$pos.phi] <- phi

    # log-prior
    alpha_prior <- 0
    beta_prior <- 0

    if(is.numeric(model$phi_prior_mu) & is.numeric(model$phi_prior_phi)){
      phi_prior <- dgamma(phi, shape=model$phi_prior_mu*model$phi_prior_phi, scale=1/model$phi_prior_phi, log=T)
    } else {
      phi_prior <- 0
    }

    # log-likelihood
    mu <- LaplacesDemon::interval(pnorm(alpha + beta*d$x), 0.0001, 0.9999)

    a <- mu*phi
    b <- (1 - mu)*phi

    LL <- sum(dbeta(d$y, a, b, log=T), na.rm=T)

    # log-posterior
    LP <- LL + alpha_prior + beta_prior + phi_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rbeta(length(d$x), a, b), parm=parm))
  }

  ld_initial_values <- c(
    qnorm(y[pre_interval][1]), # alpha: qnorm first y
    0, # beta: 0
    ((mean(y[pre_interval])*(1-mean(y[pre_interval])))/var(y[pre_interval]))-1 # phi: sample size from variance y
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior))
  )

  if(param_diagnostic) {
    full_samples$param_alpha <- rep(ld_posterior$alpha, each=length(y))
    full_samples$param_beta <- rep(ld_posterior$beta, each=length(y))
    full_samples$param_phi <- rep(ld_posterior$phi, each=length(y))
  }

  full_samples$y_model_mean <- LaplacesDemon::interval(pnorm(rep(ld_posterior$alpha, each=length(y)) + rep(ld_posterior$beta, each=length(y))*full_samples$x), 0.0001, 0.9999)
  full_samples$y_model_concentration <- rep(ld_posterior$phi, each=length(y))
  full_samples$y_model <- rbeta(
    nrow(ld_posterior)*length(y),
    full_samples$y_model_mean*full_samples$y_model_concentration,
    (1 - full_samples$y_model_mean)*full_samples$y_model_concentration
  )

  full_samples$model_iter <- rep(1:nrow(ld_posterior), each=length(y))

  return(full_samples)
}

dcc_beta_control <- function(y, x, intervention_start, ignored, covariates, model, nsamples, nburnin, param_diagnostic) {
  y <- dplyr::case_when(y < 0.0001 ~ 0.0001, y > 0.9999 ~ 0.9999, TRUE ~ y)

  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- c()
  }

  if(!is.null(ignored)) {
    pre_interval <- setdiff(pre_interval, ignored)
  }

  s <- scale(covariates)
  M <- dim(s)[[2]] # number of columns
  parm_names <- LaplacesDemon::as.parm.names(list(alpha=0, beta=rep(0, M), gamma=rep(0, M), rho=0, phi=0))

  ld_data <- list(
    y = y[pre_interval],
    x = x[pre_interval],
    s = s[pre_interval,],
    mon.names = "LP",
    parm.names = parm_names,
    pos.alpha = grep("alpha", parm_names),
    pos.beta = grep("beta", parm_names),
    pos.gamma = grep("gamma", parm_names),
    pos.rho = grep("rho", parm_names),
    pos.phi = grep("phi", parm_names)
  )

  ld_model <- function(parm, d) {
    alpha <- parm[d$pos.alpha]

    beta <- LaplacesDemon::interval(parm[d$pos.beta], 0, Inf) # no inverse covariates
    parm[d$pos.beta] <- beta

    gamma <- LaplacesDemon::interval(parm[d$pos.gamma], 0.0001, 0.9999)
    parm[d$pos.gamma] <- gamma

    rho <- LaplacesDemon::interval(parm[d$pos.rho], 0.0001, Inf)
    parm[d$pos.rho] <- rho
    rhos <- rep(rho, M)
    rhos[round(gamma) == 0] <- 0.1 # precision = 10

    phi <- LaplacesDemon::interval(parm[d$pos.phi], 0.0001, Inf)
    parm[d$pos.phi] <- phi

    # log-prior
    alpha_prior <- 0
    beta_prior <- sum(dnorm(beta, 0, rhos, log=T))

    if(!is.null(model$inclusion_probabilities)) {
      gamma_prior <- sum(dbinom(round(gamma), 1, model$inclusion_probabilities, log=T))
    } else {
      gamma_prior <- sum(dbinom(round(gamma), 1, model$nmatches/M, log=T))
    }

    rho_prior <- 0

    if(is.numeric(model$phi_prior_mu) & is.numeric(model$phi_prior_phi)){
      phi_prior <- dgamma(phi, shape=model$phi_prior_mu*model$phi_prior_phi, scale=1/model$phi_prior_phi, log=T)
    } else {
      phi_prior <- 0
    }

    # log-likelihood
    betas <- as.vector(tcrossprod(d$s, t(round(gamma)*beta)))
    mu <- LaplacesDemon::interval(pnorm(alpha + betas), 0.0001, 0.9999)

    a <- mu*phi
    b <- (1 - mu)*phi

    LL <- sum(dbeta(d$y, a, b, log=T), na.rm=T)

    # log-posterior
    LP <- LL + alpha_prior + beta_prior + gamma_prior + rho_prior + phi_prior

    return(list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rbeta(length(d$x), a, b), parm=parm))
  }

  ld_initial_values <- c(
    qnorm(mean(y[pre_interval])), # alpha: qnorm mean y
    rep(0, M), # betas: 0
    rep(1, M), # gammas: 1
    3, # rho: 3
    ((mean(y[pre_interval])*(1-mean(y[pre_interval])))/var(y[pre_interval]))-1 # phi: sample size from variance y
  )

  niter <- nburnin + nsamples
  ld_fit <- LaplacesDemon::LaplacesDemon(
    ld_model, Data = ld_data, Initial.Values = ld_initial_values,
    Iterations = niter, Status = 500, Thinning = 1, LogFile = "ld_dump.log", Algorithm = "AMWG", Specs = list(B=NULL, n=0, Periodicity = 30)
  )

  ld_posterior <- as.data.frame(ld_fit$Posterior1[(nburnin+1):niter,])
  full_samples <- data.frame(
    y = rep(y, times=nrow(ld_posterior)),
    x = rep(x, times=nrow(ld_posterior))
  )
  full_samples <- cbind(full_samples, magrittr::set_colnames(matrix(rep(t(s), nrow(ld_posterior)), ncol=M, byrow=T), paste0("s_", 1:M)))

  if(param_diagnostic) {
    full_samples$param_alpha <- rep(ld_posterior$alpha, each=length(y))
    full_samples <- cbind(full_samples, magrittr::set_colnames(matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("beta"))), each=length(y)), ncol=M, byrow=F), paste0("param_beta_", 1:M)))
    full_samples <- cbind(full_samples, magrittr::set_colnames(matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("gamma"))), each=length(y)), ncol=M, byrow=F), paste0("param_gamma_", 1:M)))
    full_samples$param_rho <- rep(ld_posterior$rho, each=length(y))
    full_samples$param_phi <- rep(ld_posterior$phi, each=length(y))
  }

  full_samples$y_model_intercept <- rep(ld_posterior$alpha, each=length(y))
  full_samples$y_model_covariates <- rowSums( matrix(rep(t(s), nrow(ld_posterior)), ncol=M, byrow=T) * ( round(matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("gamma"))), each=length(y)), ncol=M, byrow=F)) * matrix(rep(as.matrix(dplyr::select(ld_posterior, starts_with("beta"))), each=length(y)), ncol=M, byrow=F) ) )
  full_samples$y_model_mean <- LaplacesDemon::interval(pnorm(full_samples$y_model_intercept + full_samples$y_model_covariates), 0.0001, 0.9999)
  full_samples$y_model_concentration <- rep(ld_posterior$phi, each=length(y))
  full_samples$y_model <- rbeta(
    nrow(ld_posterior)*length(y),
    full_samples$y_model_mean*full_samples$y_model_concentration,
    (1 - full_samples$y_model_mean)*full_samples$y_model_concentration
  )

  full_samples$model_iter <- rep(1:nrow(ld_posterior), each=length(y))

  return(full_samples)
}

calculate_p_value <- function(y, intervention_start, ignored, full_samples) {
  if(!is.null(intervention_start)) {
    pre_interval <- 1:(intervention_start-1)
    post_interval <- intervention_start:length(y)
  } else {
    pre_interval <- 1:length(y)
    post_interval <- c()
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
