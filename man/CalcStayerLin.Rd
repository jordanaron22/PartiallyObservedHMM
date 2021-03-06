% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM.R
\name{CalcStayerLin}
\alias{CalcStayerLin}
\title{Estimates Stayer Proportion}
\usage{
CalcStayerLin(data, init, tran, class, freq_vec, pi_0, likelihoods)
}
\arguments{
\item{data}{two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data}

\item{init}{initial probabilities}

\item{tran}{transition probabilities}

\item{class}{classification probabilities}

\item{freq_vec}{frequency vector to go along with data}

\item{pi_0}{probability of being stayer}

\item{likelihoods}{list of likelihoods for each individual}
}
\value{
estimates stayer probability
}
\description{
Calculates first expectation in supplementary materials section 1.2
}
