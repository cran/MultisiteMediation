% Please edit documentation in R/msmediate.R
\name{balance}
\alias{balance}
\title{Balance checking for causal mediation analysis in multisite trials}
\usage{
balance(data, y, treatment, mediator, response, XR1, XR0, XM1, XM0, X, site)
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

\item{X}{A vector of variable names (string) of all the pretreatment covariates to be checked balance for.}

\item{site}{The variable name for the site ID (string).}
}
\value{
A list of tables containing the balance checking results for the response before weighting ($balance.R$balance1 under the treatment condition and $balance.R$balance0 under the control condition) and after weighting ($balance.R$balance1.adj under the treatment condition and $balance.R$balance0.adj under the control condition); and the balance checking results for the mediator before weighting ($balance.M$balance1 under the treatment condition and $balance.M$balance0 under the control condition) and after weighting ($balance.M$balance1.adj under the treatment condition and $balance.M$balance0.adj under the control condition). It also contains a set of balance checking plots corresponding to the tables.
\item{balance.mean}{Population average of standardized bias. The standardized bias is calculated by dividing the unweighted (before weighting) or weighted (after weighting) mean difference between response or mediator levels in each covariate by the standard deviation of the covariate }
\item{balance.sd}{Between-site standard deviation of standardized bias.}
\item{balance.lower}{Lower bound of the 95\% plausible value range of the site-specific standardized bias.}
\item{balance.upper}{Upper bound of the 95\% plausible value range of the site-specific standardized bias.}
}
\description{
This function is used to check if, within a treatment group, the estimated nonresponse weight balances the distribution of the observed covariates between the respondents and the nonrespondents, or if the estimated RMPW weight balances the distribution of the observed covariates between those whose mediator takes value 1 and those whose mediator takes value 0.
}
\examples{
data(sim.weights)

balance(data = sim.weights, y = "y", treatment = "tr", mediator = "me", response = "R", 
    XR1 = c("x1", "x2", "x3"), XR0 = c("x1", "x2", "x3"), XM1 = c("x1", "x2", "x3"), 
    XM0 = c("x1", "x2", "x3"), X = c("x1", "x2", "x3"), site = "site")
}
\references{
Qin, X., Hong, G ., Deutsch, J., & Bein, E. (under review). A Multisite causal mediation analysis in the presence of complex sample and survey designs and non-random nonresponse.
}
\author{
Xu Qin, Guanglei Hong, Jonah Deutsch, and Edward Bein
}
