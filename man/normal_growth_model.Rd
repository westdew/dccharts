% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCC.R
\name{normal_growth_model}
\alias{normal_growth_model}
\title{Normal Growth DCC Model Specification}
\usage{
normal_growth_model(
  alpha_prior_mean = NULL,
  alpha_prior_sd = NULL,
  beta_prior_mean = NULL,
  beta_prior_sd = NULL,
  sigma_prior_mode = NULL,
  sigma_prior_concentration = NULL
)
}
\arguments{
\item{alpha_prior_mean}{The prior intercept of the normal growth data model.}

\item{alpha_prior_sd}{The standard deviation as confidence of the prior intercept.}

\item{beta_prior_mean}{The prior slope of the normal growth data model.}

\item{beta_prior_sd}{The standard deviation as confidence of the prior slope.}

\item{sigma_prior_mode}{The prior standard deviation of the normal growth data model.}

\item{sigma_prior_concentration}{The confidence of the prior standard deviation.}
}
\value{
A Normal Growth DCC model specification with specified hyperparameters.
}
\description{
This function specifies a normal growth data model for a DCC.
}
