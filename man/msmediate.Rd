% Please edit documentation in R/msmediate.R
\name{msmediate}
\alias{msmediate}
\title{Multisite Causal Mediation Analysis}
\usage{
msmediate(data, y, treatment, mediator, X, site)
}
\arguments{
\item{data}{The data set for analysis.}

\item{y}{The name of the outcome variable (string).}

\item{treatment}{The name of the treatment variable (string).}

\item{mediator}{The name of the mediator variable (string).}

\item{X}{A vector of variable names (string) of pretreatment confounders, which will be included in the propensity score model. For now, the multilevel propensity score model only allows for one random intercept.}

\item{site}{The variable name for the site ID (string).}
}
\value{
A list contains the estimates of the between-site variance of direct effect, that of indirect effect, and the correlation between the direct and indirect effects across sites ($Random_effects), and the population average direct and indirect effect estimates along with their hypothesis testing results ($Fixed_effects).
}
\description{
Performs causal mediation analysis in multisite trials. It is used to estimate both the population average and between-site variance of direct and indirect effects.
}
\examples{
data(sim)

msmediate(data = sim, y = "y", treatment = "tr", mediator = "me", X = c("x1", "x2", 
    "x3"), site = "site")
}
\references{
Qin, X., & Hong, G (in press). A weighting method for assessing between-site heterogeneity in causal mediation mechanism. Journal of Educational and Behavioral Statistics.
}
\author{
Xu Qin and Guanglei Hong
}
