% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transforms.R
\name{kTransform}
\alias{kTransform}
\title{Kendall transformation}
\usage{
kTransform(x)
}
\arguments{
\item{x}{Vector or data frame to be Kendall-transformed; allowed feature types are numeric, integer (treated as numeric), ordered factor, logical and unordered factor with two or less levels.
\code{NA} and non-finite values are allowed; \code{NaN} is treated as \code{NA}.}
}
\value{
A transformed vector or data frame with transformed columns.
}
\description{
Kendall transformation
}
\examples{
kTransform(data.frame(Asc=1:3,Desc=3:1,Vsh=c(2,1,2)))
}
\references{
"Kendall transformation brings a robust categorical representation of ordinal data" M.B. Kursa. SciRep 12, 8341 (2022).
}
