% Please edit documentation in R/msmediate.R
\name{vartest.msmediate}
\alias{vartest.msmediate}
\title{Include Variance Testing for Multisite Causal Mediation Analysis}
\usage{
vartest.msmediate(data, y, treatment, mediator, X, site, npermute = 200)
}
\arguments{
\item{data}{The data set for analysis.}

\item{y}{The name of the outcome variable (string).}

\item{treatment}{The name of the treatment variable (string).}

\item{mediator}{The name of the mediator variable (string).}

\item{X}{A vector of variable names (string) of pretreatment confounders, which will be included in the propensity score model. For now, the multilevel propensity score model only allows for one random intercept.}

\item{site}{The variable name for the site ID (string).}

\item{npermute}{The number of permutations for the permutation test. The default value is 200. It may take a long time, depending on the sample size and the length of X.}
}
\value{
A list contains the hypothesis testing results of the between-site variance of the causal effects, besides the same output as given by the function msmediate().
}
\description{
Performs hypothesis testing for the between-site variance of direct effect and that of indirect effect, besides providing the same output as given by the function msmediate().
}
\examples{
data(sim)

vartest.msmediate(data = sim, y = "y", treatment = "tr", mediator = "me", X = c("x1", 
    "x2", "x3"), site = "site", npermute = 2)

}
\references{
Qin, X., & Hong, G (in press). A weighting method for assessing between-site heterogeneity in causal mediation mechanism. Journal of Educational and Behavioral Statistics.
}
\author{
Xu Qin and Guanglei Hong
}
