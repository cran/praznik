# Feature selection algorithms

#' Mutual information maximisation filter
#'
#' Calculates mutual information between all features and the decision, then returns top k.
#' @template input
#' @template y
#' @template k
#' @template output
#' @examples data(MadelonD)
#' MIM(MadelonD$X,MadelonD$Y,20)
#' @export
MIM<-function(X,Y,k=3,threads=0)
 .Call(C_MIM,X,Y,as.integer(k),as.integer(threads))

#' Minimal conditional mutual information maximisation filter
#'
#' The method starts with a feature of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\min(I(X;Y),\min_{W\in S} I(X;Y|W)),}
#' where \eqn{S} is the set of already selected features.
#' @references "Fast Binary Feature Selection using Conditional Mutual Information Maximisation" F. Fleuret, JMLR (2004)
#' @references "Object recognition with informative features and linear classification" M. Vidal-Naquet and S. Ullman, IEEE Conference on Computer Vision and Pattern Recognition (2003).
#' @template input
#' @template y
#' @template k
#' @template output
#' @examples data(MadelonD)
#' CMIM(MadelonD$X,MadelonD$Y,20)
#' @export
CMIM<-function(X,Y,k=3,threads=0)
 .Call(C_CMIM,X,Y,as.integer(k),as.integer(threads))

#' Conditional mutual information maximisation filter
#'
#' The method starts with a feature of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=I(X;Y|S),}
#' where \eqn{S} is the set of already selected features.
#' @template input
#' @template y
#' @template k
#' @template output
#' @examples data(MadelonD)
#' CMI(MadelonD$X,MadelonD$Y,20)
#' @export
CMI<-function(X,Y,k=3,threads=0)
 .Call(C_CMI,X,Y,as.integer(k),as.integer(threads))

#' Minimum redundancy maximal relevancy filter
#'
#' The method starts with a feature of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=I(X;Y)-\frac{1}{|S|}\sum_{W\in S} I(X;W),}
#' where \eqn{S} is the set of already selected features.
#' @references "Feature Selection Based on Mutual Information: Criteria of Max-Dependency, Max-Relevance, and Min-Redundancy" H. Peng et al. IEEE Pattern Analysis and Machine Intelligence (PAMI) (2005)
#' @template input
#' @template y
#' @template k
#' @param positive If true, algorithm won't return features with negative scores (i.e., with redundancy term higher than the relevance therm).
#'  In that case, \code{k} controls the maximal number of returned features, and is set to `ncol(X)` by default.
#' @template output
#' @examples data(MadelonD)
#' MRMR(MadelonD$X,MadelonD$Y,20)
#' @export
MRMR<-function(X,Y,k=if(positive) ncol(X) else 3,positive=FALSE,threads=0)
 .Call(C_MRMR,X,Y,as.integer(k),as.integer(threads),as.logical(positive))

#' Joint mutual information filter
#'
#' The method starts with a feature of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\sum_{W\in S} I(X,W;Y),}
#' where \eqn{S} is the set of already selected features.
#' @note \code{\link{DISR}} is a normalised version of JMI; \code{\link{JMIM}} and \code{\link{NJMIM}} are modifications of JMI and DISR in which minimal joint information over already selected features is used instead of a sum.
#' @references "Data Visualization and Feature Selection: New Algorithms for Nongaussian Data" H. Yang and J. Moody, NIPS (1999)
#' @template input
#' @template y
#' @template k
#' @template output
#' @examples data(MadelonD)
#' JMI(MadelonD$X,MadelonD$Y,20)
#' @export
JMI<-function(X,Y,k=3,threads=0)
 .Call(C_JMI,X,Y,as.integer(k),as.integer(threads))

#' Third-order joint mutual information filter
#'
#' The method starts with two features: \eqn{X_1} of a maximal mutual information with the decision \eqn{Y}, and \eqn{X_2} of a maximal value of \eqn{I(X_1,X_2;Y)}, as would be selected second by a regular \code{\link{JMI}}.
#' Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\frac{1}{2}\sum_{(U,W)\in S^2; U\neq W} I(X,U,W;Y),}
#' where \eqn{S} is the set of already selected features.
#' @note This method has a complexity of \eqn{O(k^2\cdot m \cdot n)}, while other filters have \eqn{O(k\cdot m \cdot n)} --- for larger \eqn{k}, it will be substantially slower.
#' In the original paper, special shrinkage estimator of MI is used; in praznik, all algorithms use ML estimators, so is \code{JMI3}.
#' @references "Efficient feature selection using shrinkage estimators" K. Sechidis, L. Azzimonti, A. Pocock, G. Corani, J. Weatherall and G. Brown. Machine Learning, 108 (8-9), pp. 1261-1286 (2019)
#' @template input
#' @template y
#' @template k
#' @template output
#' @examples \dontrun{data(MadelonD)
#' JMI3(MadelonD$X,MadelonD$Y,20)}
#' @export
JMI3<-function(X,Y,k=3,threads=0)
 .Call(C_JMI3,X,Y,as.integer(k),as.integer(threads))

#' Double input symmetrical relevance filter
#'
#' The method starts with a feature of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\sum_{W\in S} \frac{I(X,W;Y)}{H(X,W,Y)},}
#' where \eqn{S} is the set of already selected features.
#' @note DISR is a normalised version of \code{\link{JMI}}; \code{\link{JMIM}} and \code{\link{NJMIM}} are modifications of JMI and DISR in which minimal joint information over already selected features is used instead of a sum.
#' @references "On the Use of Variable Complementarity for Feature Selection in Cancer Classification" P. Meyer and G. Bontempi, (2006)
#' @template input
#' @template y
#' @template k
#' @template output
#' @examples data(MadelonD)
#' DISR(MadelonD$X,MadelonD$Y,20)
#' @export
DISR<-function(X,Y,k=3,threads=0)
 .Call(C_DISR,X,Y,as.integer(k),as.integer(threads))

#' Minimal joint mutual information maximisation filter
#'
#' The method starts with a feature of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\min_{W\in S} I(X,W;Y),}
#' where \eqn{S} is the set of already selected features.
#' @note \code{\link{NJMIM}} is a normalised version of JMIM; \code{\link{JMI}} and \code{\link{DISR}} are modifications of JMIM and NJMIM in which a sum of joint information over already selected features is used instead of a minimum.
#' @template input
#' @template y
#' @template k
#' @template output
#' @examples data(MadelonD)
#' JMIM(MadelonD$X,MadelonD$Y,20)
#' @references "Feature selection using Joint Mutual Information Maximisation" M. Bennasar, Y. Hicks and R. Setchi, (2015)
#' @export
JMIM<-function(X,Y,k=3,threads=0)
 .Call(C_JMIM,X,Y,as.integer(k),as.integer(threads))

#' Minimal normalised joint mutual information maximisation filter
#'
#' The method starts with a feature of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\min_{W\in S} \frac{I(X,W;Y)}{H(X,W,Y)},}
#' where \eqn{S} is the set of already selected features.
#' @note NJMIM is a normalised version of \code{\link{JMIM}}; \code{\link{JMI}} and \code{\link{DISR}} are modifications of JMIM and NJMIM in which a sum of joint information over already selected features is used instead of a minimum.
#' It stops returning features when the best score reaches 0.
#' @template input
#' @template y
#' @template k
#' @template output
#' @examples data(MadelonD)
#' NJMIM(MadelonD$X,MadelonD$Y,20)
#' @references "Feature selection using Joint Mutual Information Maximisation" M. Bennasar, Y. Hicks and R. Setchi, (2015)
#' @export
NJMIM<-function(X,Y,k=3,threads=0)
 .Call(C_NJMIM,X,Y,as.integer(k),as.integer(threads))

#' Joint impurity filter
#'
#' The method starts with a feature of a maximal impurity gain with the decision \eqn{Y}.
#' Then, it greedily adds feature \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\sum_{W\in S} G(X,W;Y),}
#' where \eqn{S} is the set of already selected features, and
#' \deqn{G(X;Y)=\sum_{xy}\frac{p_{xy}^2}{p_x}-\sum_{y} p_y^2}
#' is the Gini impurity gain from partitioning \eqn{Y} according to \eqn{X}.
#' @note This is an impurity-based version of \code{\link{JMI}}; expect similar results in slightly shorter time.
#' @template input
#' @template y
#' @template k
#' @template output
#' @examples data(MadelonD)
#' JIM(MadelonD$X,MadelonD$Y,20)
#' @export
JIM<-function(X,Y,k=3,threads=0)
 .Call(C_JIM,X,Y,as.integer(k),as.integer(threads))
