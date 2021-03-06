% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM.R
\name{CalcInitialLin}
\alias{CalcInitialLin}
\title{Estimates Initial Probabilities}
\usage{
CalcInitialLin(data, freq_vec, pi_0, forw, backw, likelihoods, max2, max3)
}
\arguments{
\item{data}{two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data}

\item{freq_vec}{frequency vector to go along with data}

\item{pi_0}{probability of being stayer}

\item{forw}{forward quantity for each individual}

\item{backw}{backward quantity for each individual}

\item{likelihoods}{list of likelihoods for each individual}

\item{max2}{number of values the second test can take on}

\item{max3}{number of values the third test can take on}
}
\value{
estimates vector of initial probabilities
}
\description{
Calculates second expectation in section 2.6
}
