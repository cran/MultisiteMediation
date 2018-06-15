% Please edit documentation in R/msmediate.R
\name{msmediate.weights}
\alias{msmediate.weights}
\title{Causal mediation analysis in multisite trials in the presence of complex sample and survey designs and non-random nonresponse}
\usage{
msmediate.weights(data, y, treatment, mediator, response, XR1, XR0, XM1, XM0, site, 
    sample.weight)
}
\arguments{
\item{data}{The data set for analysis.}

\item{y}{The name of the outcome variable (string).}

\item{treatment}{The name of the treatment variable (string).}

\item{mediator}{The name of the mediator variable (string).}

\item{response}{The name of the response variable (string), which is equal to 1 if the individual responded and 0 otherwise.}

\item{XR1}{A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the treatment condition. The function will treat covariates with fewer than 10 unique values as categorical variables. For now, the multilevel propensity score model only allows for one random intercept.}

\item{XR0}{A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the control condition. The function will treat covariates with fewer than 10 unique values as categorical variables. For now, the multilevel propensity score model only allows for one random intercept.}

\item{XM1}{A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the treatment condition. The function will treat covariates with fewer than 10 unique values as categorical variables. For now, the multilevel propensity score model only allows for one random intercept.}

\item{XM0}{A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the control condition. The function will treat covariates with fewer than 10 unique values as categorical variables. For now, the multilevel propensity score model only allows for one random intercept.}

\item{site}{The variable name for the site ID (string).}

\item{sample.weight}{The variable name for the sample weight given by design (string).}
}
\value{
A list contains the estimates of the between-site variances of natural direct effect, natural indirect effect, pure indirect effect, and treatment-by-mediator interaction effect, and the correlations between the effects across sites ($Random_effects), and the population average effect estimates along with their hypothesis testing results ($Fixed_effects).
}
\description{
This function is used to estimate both the population average and between-site variance of natural direct effect, natural indirect effect, pure indirect effect, and treatment-by-mediator interaction effect. 
It incorporates a sample weight to adjust for complex sample and survey designs and employs an estimated nonresponse weight to account for non-random nonresponse.
}
\examples{
data(sim.weights)

msmediate.weights(data = sim.weights, y = "y", treatment = "tr", mediator = "me", 
    response = "R", XR1 = c("x1", "x2", "x3"), XR0 = c("x1", "x2", "x3"), XM1 = c("x1", 
        "x2", "x3"), XM0 = c("x1", "x2", "x3"), site = "site", sample.weight = "WD")
}
\references{
Qin, X., Hong, G ., Deutsch, J., & Bein, E. (under review). A Multisite causal mediation analysis in the presence of complex sample and survey designs and non-random nonresponse.
}
\author{
Xu Qin, Guanglei Hong, Jonah Deutsch, and Edward Bein
}
