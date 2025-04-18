% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCC.R
\name{normal_model}
\alias{normal_model}
\title{Normal DCC Model Specification}
\usage{
normal_model(
  alpha_prior = NULL,
  alpha_prior_confidence = NULL,
  sigma_prior = NULL,
  sigma_prior_confidence = NULL
)
}
\arguments{
\item{alpha_prior}{The prior mean of the normal data model}

\item{alpha_prior_confidence}{The confidence of the prior mean (scale: standard deviation)}

\item{sigma_prior}{The prior standard deviation of the normal data model}

\item{sigma_prior_confidence}{The confidence of the prior standard deviation (rate: sample size)}
}
\value{
A Normal DCC model specification with specified hyperparameters
}
\description{
This function specifies a stationary normal data model for a DCC.
}
