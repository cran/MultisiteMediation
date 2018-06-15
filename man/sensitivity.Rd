% Please edit documentation in R/msmediate.R
\name{sensitivity}
\alias{sensitivity}
\title{Sensitivity analysis for causal mediation analysis in multisite trials}
\usage{
sensitivity(data, y, treatment, mediator, response, XR1, XR0, XM1, XM0, omit.X = NULL, 
    ran.omit.X = 1, site)
}
\arguments{
\item{data}{The data set for analysis.}

\item{y}{The name of the outcome variable (string).}

\item{treatment}{The name of the treatment variable (string).}

\item{mediator}{The name of the mediator variable (string).}

\item{response}{The name of the response variable (string), which is equal to 1 if the individual responded and 0 otherwise.}

\item{XR1}{A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the treatment condition in the original analysis. The function will treat covariates with fewer than 10 unique values as categorical variables. For now, the multilevel propensity score model only allows for one random intercept.}

\item{XR0}{A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the control condition in the original analysis. The function will treat covariates with fewer than 10 unique values as categorical variables. For now, the multilevel propensity score model only allows for one random intercept.}

\item{XM1}{A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the treatment condition in the original analysis. The function will treat covariates with fewer than 10 unique values as categorical variables. For now, the multilevel propensity score model only allows for one random intercept.}

\item{XM0}{A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the control condition in the original analysis. The function will treat covariates with fewer than 10 unique values as categorical variables. For now, the multilevel propensity score model only allows for one random intercept.}

\item{omit.X}{A vector of variable names (string) of omitted pretreatment confounders of the response-mediator, response-outcome, or mediator-outcome relationships.}

\item{ran.omit.X}{A vector of variable names (string) of pretreatment covariates that have been included in the original analysis but whose random slopes are omitted from response models or mediator models in the original analysis. It takes value 1 if there are no omitted random slopes.}

\item{site}{The variable name for the site ID (string).}
}
\value{
A list contains sensitivity parameters and the effect size of bias due to omission of pretreatment confounders ($`Bias due to omission of confounders) or random slopes ($`Bias due to omission of random slopes) for each causal parameter ($ITT for the population average ITT effect, $var.ITT for the between-site variance of ITT effect; $NIE for the population average natural indirect effect, $var.NIE for the between-site variance of natural indirect effect; $NDE for the population average natural direct effect, $var.NDE for the between-site variance of natural direct effect; $cov.NIE.NDE for the between-site covariance between natural indirect and direct effects; $PIE for the population average pure indirect effect; and $ INT for the population average natural treatment-by-mediator interaction effect.
\item{c1}{Average conversion coefficient under the experimental condition.}
\item{sigma1.ITT}{Average standard deviation of the ITT weight discrepancy in the experimental group at the same site. It is associated with the degree to which the omissions predict response under either treatment condition.}
\item{rho1.ITT}{Average correlation between the ITT weight discrepancy and the outcome in the experimental group at the same site. It is related to the degree to which the omissions predict the outcome within levels of the response in the experimental group.}
\item{cov.sigma1.rho1.ITT}{Covariance between site-specific sigma1.ITT and rho1.ITT.}
\item{c0}{Average conversion coefficient under the control condition.}
\item{sigma0.ITT}{Average standard deviation of the ITT weight discrepancy in the control group at the same site.}
\item{rho0.ITT}{Average correlation between the ITT weight discrepancy and the outcome in the control group at the same site.}
\item{cov.sigma0.rho0.ITT}{Covariance between site-specific sigma0.ITT and rho0.ITT.}
\item{bias.ITT}{Effect size of bias of the population average ITT effect due to the omissions.}
\item{var.c.sigma.rho}{Between-site variance of the effect size of bias for site-specific ITT.}
\item{cov.b.ITT.c.sigma.rho}{Covariance between the effect size of site-specific ITT and the effect size of bias for site-specific ITT.}
\item{bias.var.ITT}{Effect size of bias of the between-site variance of the ITT effect due to the omissions.}
\item{sigma.counter}{Average standard deviation of the overal weight (product of the ITT weight and the RMPW weight) discrepancy in the experimental group at the same site. It is associated with the degree to which the omissions predict both response and mediator value assignment under either treatment condition.}
\item{rho.counter}{Average correlation between the overall weight (product of the ITT weight and the RMPW weight) discrepancy and the outcome in the experimental group at the same site. It is related to the degree to which the omissions predict the outcome within levels of the response and mediator in the experimental group.}
\item{cov.sigma.rho.counter}{Covariance between site-specific sigma.counter and rho.counter.}
\item{bias.NIE}{Effect size of bias of the population average natural indirect effect due to the omissions.}
\item{var.c.sigma.rho.NIE}{Between-site variance of the effect size of bias for site-specific NIE.}
\item{cov.b.NIE.c.sigma.rho}{Covariance between the effect size of site-specific NIE and the effect size of bias for site-specific NIE.}
\item{bias.var.NIE}{Effect size of bias of the between-site variance of the natural indirect effect due to the omissions.}
\item{bias.NDE}{Effect size of bias of the population average natural direct effect due to the omissions.}
\item{var.c.sigma.rho.NDE}{Between-site variance of the effect size of bias for site-specific NDE.}
\item{cov.b.NDE.c.sigma.rho}{Covariance between the effect size of site-specific NDE and the effect size of bias for site-specific NDE.}
\item{bias.var.NDE}{Effect size of bias of the between-site variance of the natural direct effect due to the omissions.}
\item{bias.cov.NIE.NDE}{Effect size of bias of the between-site covariance between the natural direct and indirect effects due to the omissions.}
\item{sigma.PIE}{Average standard deviation of the overal weight (product of the ITT weight and the RMPW weight) discrepancy in the control group at the same site.}
\item{rho.PIE}{Average correlation between the overall weight (product of the ITT weight and the RMPW weight) discrepancy and the outcome in the control group at the same site.}
\item{cov.sigma.rho.PIE}{Covariance between site-specific sigma.PIE and rho.PIE.}
\item{bias.PIE}{Effect size of bias of the population average pure indirect effect due to the omissions.}
\item{bias.INT}{Effect size of bias of the population average natural treatment-by-mediator interactioin effect due to the omissions.}
}
\description{
This function is used to calculate the effect size of the hidden bias associated with one or more omitted confounders or omitted random slopes of existing confounders, based on a small number of weighting-based sensitivity parameters.
}
\examples{
data(sim.weights)

sensitivity(data = sim.weights, y = "y", treatment = "tr", mediator = "me", 
    response = "R", XR1 = "x2", XR0 = "x2", XM1 = c("x1", "x2"), XM0 = "x2", 
    omit.X = c("x1", "x3"), ran.omit.X = "x2", site = "site")
}
\references{
Qin, X., Hong, G ., Deutsch, J., & Bein, E. (under review). A Multisite causal mediation analysis in the presence of complex sample and survey designs and non-random nonresponse.
}
\author{
Xu Qin, Guanglei Hong, Jonah Deutsch, and Edward Bein
}
