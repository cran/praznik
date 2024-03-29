% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithms.R
\name{JMI}
\alias{JMI}
\title{Joint mutual information filter}
\usage{
JMI(X, Y, k = 3, threads = 0)
}
\arguments{
\item{X}{Attribute table, given as a data frame with either factors (preferred), booleans, integers (treated as categorical) or reals (which undergo automatic categorisation; see below for details).
Single vector will be interpreted as a data.frame with one column.
\code{NA}s are not allowed.}

\item{Y}{Decision attribute; should be given as a factor, but other options are accepted, exactly like for attributes.
\code{NA}s are not allowed.}

\item{k}{Number of attributes to select.
Must not exceed \code{ncol(X)}.}

\item{threads}{Number of threads to use; default value, 0, means all available to OpenMP.}
}
\value{
A list with two elements: \code{selection}, a vector of indices of the selected features in the selection order, and \code{score}, a vector of corresponding feature scores.
Names of both vectors will correspond to the names of features in \code{X}.
Both vectors will be at most of a length \code{k}, as the selection may stop sooner, even during initial selection, in which case both vectors will be empty.
}
\description{
The method starts with a feature of a maximal mutual information with the decision \eqn{Y}.
Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
\deqn{J(X)=\sum_{W\in S} I(X,W;Y),}
where \eqn{S} is the set of already selected features.
}
\note{
\code{\link{DISR}} is a normalised version of JMI; \code{\link{JMIM}} and \code{\link{NJMIM}} are modifications of JMI and DISR in which minimal joint information over already selected features is used instead of a sum.

The method requires input to be discrete to use empirical estimators of distribution, and, consequently, information gain or entropy.
To allow smoother user experience, praznik automatically coerces non-factor vectors in inputs, which requires additional time, memory and may yield confusing results -- the best practice is to convert data to factors prior to feeding them in this function.
Real attributes are cut into about 10 equally-spaced bins, following the heuristic often used in literature.
Precise number of cuts depends on the number of objects; namely, it is \eqn{n/3}, but never less than 2 and never more than 10.
Integers (which technically are also numeric) are treated as categorical variables (for compatibility with similar software), so in a very different way -- one should be aware that an actually numeric attribute which happens to be an integer could be coerced into a \eqn{n}-level categorical, which would have a perfect mutual information score and would likely become a very disruptive false positive.
}
\examples{
data(MadelonD)
JMI(MadelonD$X,MadelonD$Y,20)
}
\references{
"Data Visualization and Feature Selection: New Algorithms for Nongaussian Data" H. Yang and J. Moody, NIPS (1999)
}
