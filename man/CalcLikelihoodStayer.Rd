% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM.R
\name{CalcLikelihoodStayer}
\alias{CalcLikelihoodStayer}
\title{Calculates Likelihood of Stayers}
\usage{
CalcLikelihoodStayer(data, class, pi)
}
\arguments{
\item{data}{two dimensional matrix for individual.  1st dimension is time 2nd is partially observed data}

\item{class}{classification matrix}

\item{pi}{probability of being a stayer}
}
\value{
stayer likelihood
}
\description{
Calculates Likelihood of Stayers
}
