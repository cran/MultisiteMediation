#' A simulated example data
#' 
#' This simulated data list is for demonstration.
#' 
#' @docType data
#' @name sim
#' @return A list containing
#' \item{site}{Site ID}
#' \item{y}{Outcome}
#' \item{tr}{Treatment}
#' \item{me}{Mediator}
#' \item{x1}{Pretreatment covariate}
#' \item{x2}{Pretreatment covariate}
#' \item{x3}{Pretreatment covariate}
#' \item{x4}{Pretreatment covariate}
NULL

#' Causal mediation analysis in multisite trials
#' 
#' This function is used to estimate both the population average and between-site variance of direct and indirect effects. 
#'
#' @param data The data set for analysis.
#' @param y The name of the outcome variable (string).
#' @param treatment The name of the treatment variable (string).
#' @param mediator The name of the mediator variable (string).
#' @param X A vector of variable names (string) of pretreatment covariates, which will be included in the propensity score model. For now, the multilevel propensity score model only allows for one random intercept.
#' @param site The variable name for the site ID (string).
#' @return A list contains the estimates of the between-site variance of direct effect, that of indirect effect, and the correlation between the direct and indirect effects across sites ($Random_effects), and the population average direct and indirect effect estimates along with their hypothesis testing results ($Fixed_effects).
#' @author Xu Qin and Guanglei Hong
#' @references Qin, X., & Hong, G (2017). A weighting method for assessing between-site heterogeneity in causal mediation mechanism. Journal of Educational and Behavioral Statistics. Journal of Educational and Behavioral Statistics, 42(3), 308-340. \doi{10.3102/1076998617694879}
#' @export
#' @importFrom stats as.formula binomial coef fitted lm model.matrix pnorm
#' @importFrom lme4 VarCorr fixef glmer ranef
#' @importFrom statmod gauss.quad.prob
#' @examples 
#' data(sim)
#'
#' msmediate(data = sim, y = "y", treatment = "tr", mediator = "me", X = c("x1", "x2", "x3", "x4"), site = "site")
#'                     
msmediate = function(data, y, treatment, mediator, X, site) {
  data = as.data.frame(data)
  colnames(data)[which(colnames(data) == treatment)] = "tr"
  colnames(data)[which(colnames(data) == mediator)] = "me"
  colnames(data)[which(colnames(data) == site)] = "site"
  colnames(data)[which(colnames(data) == y)] = "y"
  ranX = 1
  data = data[order(data$site), ]
  # # Factorize categorical covariates (with fewer than 10 categories)
  # for(i in 1:length(X)){
  #   if(length(unique(data[, X[i]])) > 2 & length(unique(data[, X[i]])) < 10){
  #     data[, X[i]] = as.factor(data[, X[i]])
  #   }
  # }
  # covariates = model.matrix(as.formula(paste("~", paste(X, collapse = "+"))), data)[, -1]
  # data = data[, -which(colnames(data) %in% X)]
  # data = cbind(data, covariates)
  data1 = data[data$tr == 1, ]
  data0 = data[data$tr == 0, ]
  # X = colnames(covariates)
  
  l_tr = suppressWarnings({glmer(as.formula(paste("me", "~", paste(X, collapse="+"), "+(", ranX, "|site)")), data = data1, family = binomial, nAGQ = 10)})
  data$p1[data$tr == 1] = fitted(l_tr)
  l_ctrl = suppressWarnings({glmer(as.formula(paste("me", "~", paste(X, collapse="+"), "+(", ranX, "|site)")), data = data0, family = binomial, nAGQ = 10)})
  data$p0[data$tr == 0] = fitted(l_ctrl)
  
  predict.lmer = function(l, data, X, ranX) {
    nj = NULL
    for (j in as.numeric(rownames(ranef(l)$site))) {
      nj = c(nj, sum(data$site == j))
    }
    ranef = NULL
    for (k in 1:dim(ranef(l)$site)[2]) {
      ranef = cbind(ranef, rep(ranef(l)$site[, k], nj))
    }
    pred_logit = X %*% fixef(l) + apply(ranX * ranef, 1, sum)
    pred = exp(pred_logit)/(exp(pred_logit) + 1)
  }
  data$p0[data$tr == 1] = predict.lmer(l_ctrl, data1, model.matrix(l_tr), rep(1, nrow(data1)))
  data$p1[data$tr == 0] = predict.lmer(l_tr, data0, model.matrix(l_ctrl), rep(1, nrow(data0)))
  
  data$rmpw[data$tr == 1 & data$me == 1] = (data$p0/data$p1)[data$tr == 1 & data$me == 1]
  data$rmpw[data$tr == 1 & data$me == 0] = ((1 - data$p0)/(1 - data$p1))[data$tr == 1 & data$me == 0]
  data$rmpw[data$tr == 0] = 1

  data = data[order(data$site),]
  data0 = data[data$tr == 0, ]
  data1 = data[data$tr == 1, ]
  
  f = function(x, alpha, s, F, B){
    return(1 / (1 + exp(-(x %*% alpha + s * F * B)))) # When the estimated variance of both two or one of the two random effects is 1
  }
  
  nnodes = 10
  temp = gauss.quad.prob(nnodes, "normal") # Here, random effects are assumed normal
  nodes = temp$nodes
  weights = temp$weights
  # When there is one random effect, i.e. when the estimated variance of one random effect is 0
  B_1 = t(as.matrix(nodes))
  A_1 = weights
  
  x = model.matrix(as.formula(paste("me", "~", paste(X, collapse="+"))), data = data) # Here we have the same covariates in the models for both the treatment and control group
  h1_fun = function(data0, data){
    if(data0$tr[1] == 0){
      l0 = l_ctrl
    } else{
      l0 = l_tr
    }        
    alpha0 = fixef(l0)
    if(round(VarCorr(l0)$site[1], 5) == 0){
      F0 = 0
      v0 = 0
      s0 = matrix(0, nrow(x), 1)
      J0 = 0
      A0 = 1
      B0 = as.matrix(0)
    } else {
      F0 = sqrt(VarCorr(l0)$site[1])
      v0 = ranef(l0)$site
      s0 = matrix(1, nrow(x), 1)
      J0 = 1
      A0 = A_1
      B0 = B_1
    }
    p0 = NULL
    for(q in 1:length(A0)){
      p0 = cbind(p0, f(x = x, alpha = alpha0, s = s0, F = F0, B = B0[, q]))
    }
    p0 = as.matrix(p0)
    
    h_alpha0 = NULL
    h_F0 = NULL
    for(j in unique(data$site)){
      l0j = 0
      h_alpha0j = 0
      h_F0j = 0
      for(q in 1:length(A0)){
        p0q = p0[, q]
        if(data0$tr[1] == 0){
          f0jq = 1
          for(i in which(data$site == j)){
            f0jq = f0jq * (p0q[i]^data$me[i] * (1 - p0q[i])^(1 - data$me[i]))^(1 - data$tr[i])
          }
          l0j = l0j + f0jq * A0[q]
          h_alpha0j = h_alpha0j + as.numeric(f0jq * A0[q] * (1 - data$tr[data$site == j]) * (data$me[data$site == j] - p0q[data$site == j])) * x[data$site == j,]      
          h_F0j =  h_F0j + as.numeric(f0jq * A0[q] * (1 - data$tr[data$site == j]) * (data$me[data$site == j] - p0q[data$site == j])) * kronecker(t(B0[, q]), s0[data$site == j,]) %*% t(J0)    
        } else {
          f0jq = 1
          for(i in which(data$site == j)){
            f0jq = f0jq * (p0q[i]^data$me[i] * (1 - p0q[i])^(1 - data$me[i]))^(data$tr[i])
          }
          l0j = l0j + f0jq * A0[q]
          h_alpha0j = h_alpha0j + as.numeric(f0jq * A0[q] * data$tr[data$site == j] * (data$me[data$site == j] - p0q[data$site == j])) * x[data$site == j,]              
          h_F0j =  h_F0j + as.numeric(f0jq * A0[q] * data$tr[data$site == j] * (data$me[data$site == j] - p0q[data$site == j])) * kronecker(t(B0[, q]), s0[data$site == j,]) %*% t(J0)    
        }        
      }
      h_alpha0 = rbind(h_alpha0, h_alpha0j/l0j)
      h_F0 = rbind(h_F0, h_F0j/l0j)
    }
    
    if(round(VarCorr(l0)$site[1], 5) == 0){
      h1_0 = h_alpha0
    } else {
      h1_0 = cbind(h_alpha0, h_F0)
    }
    
    if(round(VarCorr(l0)$site[1], 5) == 0){
      theta0 = matrix(0, 1, length(unique(data0$site)))
    } else {
      theta0 = solve(F0) %*% t(v0)
    }
    
    return(list(h1 = h1_0, theta = theta0, s = s0, J = J0))
  }
  
  h1_list0 = h1_fun(data0, data)
  h1_list1 = h1_fun(data1, data)
  h1_0 = h1_list0$h1 
  h1_1 = h1_list1$h1
  h1 = cbind(h1_0, h1_1)
  p = ncol(h1)
  theta0 = h1_list0$theta
  theta1 = h1_list1$theta
  s0 = h1_list0$s
  s1 = h1_list1$s
  J0 = h1_list0$J
  J1 = h1_list1$J
  
  N = nrow(data)
  J = length(unique(data$site))
  h2 = matrix(0, N, J * 3)
  mu = NULL
  nj = NULL
  for (j in unique(data$site)) {
    dataj = data[data$site == j, ]
    nj = c(nj, nrow(dataj))
    mu_0 = mean(dataj$y[dataj$tr == 0])
    mu_counter = sum(dataj$y * dataj$rmpw * dataj$tr)/sum(dataj$rmpw * dataj$tr)
    mu_1 = mean(dataj$y[dataj$tr == 1])
    h_0 = (data$y - mu_0) * (1 - data$tr) * (data$site == j)
    h_counter = (data$y - mu_counter) * data$rmpw * data$tr * (data$site == j)
    h_1 = (data$y - mu_1) * data$tr * (data$site == j)
    h2[, (3 * (which(unique(data$site) == j) - 1) + 1):(3 * which(unique(data$site) == j))] = cbind(h_0, 
                                                                                                    h_counter, h_1)
    mu = c(mu, mu_0, mu_counter, mu_1)
  }
  
  h = cbind(h1, h2)  # Stack the moment functions from both steps
  H = 1/N * t(as.matrix(h)) %*% as.matrix(h)  # See Equation S13
  ## The derivation of R (See Equation S14)
  bdiag = function(A, B){
    C = matrix(0, nrow = nrow(A) + nrow(B), ncol = ncol(A) + ncol(B))
    C[1:nrow(A),1:ncol(A)] = A
    C[(1 + nrow(A)):nrow(C),(1 + ncol(A)):ncol(C)] = B
    return(C)
  }
  R11 = as.matrix(bdiag(-1/N * t(h1_0) %*% as.matrix(h1_0), -1/N * t(h1_1) %*% as.matrix(h1_1)))
  
  R22_diag = NULL
  for (j in unique(data$site)) {
    R22_diag = c(R22_diag, mean(-(1 - data$tr) * (data$site == j)), mean(-data$rmpw * data$tr * (data$site == j)), mean(-data$tr * (data$site == j)))
  }
  R22 = diag(R22_diag)
  
  p0 = data$p0
  p1 = data$p1
  R21 = NULL
  for (j in 1:J) {
    dh_counter_dalpha0 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (data$me/p1 - (1 - data$me)/(1 - p1)) *   
                                 p0 * (1 - p0) * x, 2, mean)
    dh_counter_dF0 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (data$me/p1 - (1 - data$me)/(1 - p1)) * p0 
                           * (1 - p0) * kronecker(t(theta0[, j]), s0) %*% t(J0), 2, mean)
    dh_counter_dalpha1 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (-data$me * p0/p1 * (1 - p1) + (1 - data$me) 
                                                                  * (1 - p0)/(1 - p1) * p1) * x, 2, mean)
    dh_counter_dF1 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (-data$me * p0/p1 * (1 - p1) + (1 - data$me) 
                                                              * (1 - p0)/(1 - p1) * p1) * kronecker(t(theta1[, j]), s1) %*% t(J1), 2, mean)
    if (ncol(h1_0) - ncol(x) == 0 & ncol(h1_1) - ncol(x) == 0) {
      R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dalpha1), rep(0, p)))
    } else if (ncol(h1_0) - ncol(x) == 0) {
      R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dalpha1, dh_counter_dF1), 
                             rep(0, p)))
    } else if (ncol(h1_1) - ncol(x) == 0) {
      R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dF0, dh_counter_dalpha1), 
                             rep(0, p)))
    } else {
      R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dF0, dh_counter_dalpha1, 
                                          dh_counter_dF1), rep(0, p)))
    }
  }
  
  R12 = matrix(0, p, 3 * J)
  
  R = matrix(0, p + 3 * J, p + 3 * J)
  R[1:p, 1:p] = R11
  R[1:p, (p + 1):(p + 3 * J)] = R12
  R[(p + 1):(p + 3 * J), 1:p] = R21
  R[(p + 1):(p + 3 * J), (p + 1):(p + 3 * J)] = R22
  
  V_all = 1/N * solve(R) %*% H %*% t(solve(R))  # See Appendix A
  V_mu = V_all[(ncol(V_all) - 3 * J + 1):ncol(V_all), (ncol(V_all) - 3 * J + 1):ncol(V_all)]  
  
  Phi = matrix(c(-1, 0, 1, -1, 0, 1), 2, 3)
  for (j in 1:(J - 1)) {
    Phi = as.matrix(bdiag(Phi, matrix(c(-1, 0, 1, -1, 0, 1), 2, 3)))
  }
  beta = Phi %*% mu # Estimated site-specific direct and indirect effects
  V_beta = Phi %*% V_mu %*% t(Phi)
  dej = beta[seq(1, 2 * J, 2)]
  iej = beta[seq(2, 2 * J, 2)]
  
  Psi = cbind(rep(c(1, 0), J), rep(c(0, 1), J))
  gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% beta
  de = gamma[1]
  ie = gamma[2]
  G = 0
  V = 0
  for(j in 1:J){
    G = G + (beta[(2 * j - 1): (2 * j)] - gamma) %*% t(beta[(2 * j - 1): (2 * j)] - gamma)    
    Vj = V_beta[(2 * j - 1):(2 * j), (2 * j - 1):(2 * j)]
    V = V + Vj
  }
  tau = 1/(J - 1) * (G - V + 1/J * t(Psi) %*% V_beta %*% Psi)
  
  V_gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% (V_beta + kronecker(diag(rep(1,J)), tau)) %*% Psi %*%   
    solve(t(Psi) %*% Psi)
  V_de = V_gamma[1, 1]
  V_ie = V_gamma[2, 2]   
  
  se_de = sqrt(V_de)
  se_ie = sqrt(V_ie)
  
  z_de = de/se_de
  z_ie = ie/se_ie
  
  p_de = (1 - pnorm(abs(z_de))) * 2
  p_ie = (1 - pnorm(abs(z_ie))) * 2
  
  est = cbind(round(c(de, ie), 5), round(c(se_de, se_ie), 5), round(c(z_de, z_ie), 3), round(c(p_de, p_ie), 3))
  
  sig = NULL
  sig[est[, 4] <= 0.001] = "**"
  sig[est[, 4] > 0.001 & est[, 4] <= 0.01] = "*"
  sig[est[, 4] > 0.01 & est[, 4] <= 0.05] = "."
  sig[est[, 4] > 0.05] = ""
  
  est = cbind(est, sig)
  est[est[, 4] < 0.001, 4] = "<0.001"
  est = as.data.frame(est)
  colnames(est) = c("Estimate", "Std.Error", "t value", "Pr(>|t|)", "")
  rownames(est) = c("Natrual Direct Effect", "Natrual Indirect Effect")
  
  tau.ori = tau
  tau.sub = tau
  
  if("TRUE" %in% (diag(tau.ori) < 0)){
    tau[which(diag(tau < 0)), ] = 0
    tau[, which(diag(tau < 0))] = 0
    tau.sub = tau.ori[-which(diag(tau.ori < 0)), ]
    tau.sub = tau.sub[, -which(diag(tau.sub < 0))]    
  }
  
  cor.sub = cov2cor(tau.sub)
  if(max(na.omit(cor.sub)) > 1)
    cor.sub = suppressWarnings(cor.smooth(cor.sub))
  cor = cor.sub
  
  if("TRUE" %in% (diag(tau.ori) < 0)){
    cor = tau.ori
    cor[which(diag(tau.ori < 0)), ] = NaN
    cor[, which(diag(tau.ori < 0))] = NaN
    cor[which(diag(tau.ori >= 0)), which(diag(tau.ori >= 0))] = cor.sub
  }
  
  cor_de_ie = cor[2, 1]
  
  var = cbind(round(c(tau[1, 1], tau[2, 2]), 3), round(c(sqrt(tau[1, 1]), sqrt(tau[2, 2])), 3), c("", round(cor_de_ie, 3)))
  var = as.data.frame(var)
  colnames(var) = c("Variance", "Std.Dev.", "Corr")
  rownames(var) = c("Natrual Direct Effect", "Natrual Indirect Effect")
  
  result = list(Random_effects = var, Fixed_effects = est)
  
  return(result)
}

#' A simulated example data
#' 
#' This simulated data list is for demonstration.
#' 
#' @docType data
#' @name sim.weights
#' @return A list containing
#' \item{site}{Site ID}
#' \item{y}{Outcome}
#' \item{tr}{Treatment}
#' \item{me}{Mediator}
#' \item{x1}{Pretreatment covariate}
#' \item{x2}{Pretreatment covariate}
#' \item{x3}{Pretreatment covariate}
#' \item{R}{Response indicator}
#' \item{WD}{Sample weight}
NULL

#' Causal mediation analysis in multisite trials in the presence of complex sample and survey designs and non-random nonresponse
#' 
#' This function is used to estimate both the population average and between-site variance of natural direct effect, natural indirect effect, pure indirect effect, and treatment-by-mediator interaction effect. 
#' It incorporates a sample weight to adjust for complex sample and survey designs and employs an estimated nonresponse weight to account for non-random nonresponse.
#'
#' @param data The data set for analysis.
#' @param y The name of the outcome variable (string).
#' @param treatment The name of the treatment variable (string).
#' @param mediator The name of the mediator variable (string).
#' @param response The name of the response variable (string), which is equal to 1 if the individual responded and 0 otherwise.
#' @param XR1 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XR0 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM1 A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM0 A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param site The variable name for the site ID (string).
#' @param sample.weight The variable name for the sample weight given by design (string).
#' @return A list contains the estimates of the between-site variances of natural direct effect, natural indirect effect, pure indirect effect, and treatment-by-mediator interaction effect, and the correlations between the effects across sites ($Random_effects), and the population average effect estimates along with their hypothesis testing results ($Fixed_effects).
#' @author Xu Qin, Guanglei Hong, Jonah Deutsch, and Edward Bein
#' @references Qin, X., Hong, G., Deutsch, J., & Bein, E. (2019). Multisite causal mediation analysis in the presence of complex sample and survey designs and non-random non-response. Journal of the Royal Statistical Society: Series A (Statistics in Society), 182(4), 1343-1370. \doi{10.1111/rssa.12446}
#' @export
#' @importFrom stats as.formula binomial coef fitted lm model.matrix pnorm na.omit cov2cor
#' @importFrom lme4 VarCorr fixef glmer ranef glmerControl
#' @importFrom statmod gauss.quad.prob
#' @importFrom psych cor.smooth
#' @importFrom MASS ginv
#' @examples 
#' data(sim.weights)
#'
#' msmediate.weights(data = sim.weights, y = "y", treatment = "tr", mediator = "me", response = "R", XR1 = c("x1", "x2", "x3"), XR0 = c("x1", "x2", "x3"), XM1 = c("x1", "x2", "x3"), XM0 = c("x1", "x2", "x3"), site = "site", sample.weight = "WD")
#'   
msmediate.weights = function(data, y, treatment, mediator, response, XR1, XR0, XM1, XM0, site, sample.weight){
  data = as.data.frame(data)
  colnames(data)[which(colnames(data) == y)] = "y"
  colnames(data)[which(colnames(data) == treatment)] = "tr"
  colnames(data)[which(colnames(data) == mediator)] = "M"
  colnames(data)[which(colnames(data) == response)] = "R"
  colnames(data)[which(colnames(data) == site)] = "site"
  colnames(data)[which(colnames(data) == sample.weight)] = "WD"
  data = data[order(data$site), ]
  # # Factorize categorical covariates (with fewer than 10 categories)
  # transform = function(X){
  #   for(i in 1:length(X)){
  #     if(length(unique(data[, X[i]])) > 2 & length(unique(data[, X[i]])) < 10){
  #       data[, X[i]] = as.factor(data[, X[i]])
  #     }
  #   }
  #   covariates = model.matrix(as.formula(paste("~", paste(X, collapse = "+"))), data)
  #   X = colnames(covariates)
  #   return(list(covariates = covariates[, -1], X = X[-1]))
  # }
  # transform.XR1 = transform(XR1)
  # transform.XR0 = transform(XR0)
  # transform.XM1 = transform(XM1)
  # transform.XM0 = transform(XM0)
  # data = data[, -which(colnames(data) %in% unique(c(XR1, XR0, XM1, XM0)))]
  # XR1 = transform.XR1$X
  # XR0 = transform.XR0$X
  # XM1 = transform.XM1$X
  # XM0 = transform.XM0$X
  # covariates = cbind(transform.XR1$covariates, transform.XR0$covariates, transform.XM1$covariates, transform.XM0$covariates)
  # colnames(covariates) = c(XR1, XR0, XM1, XM0)
  # data = cbind(data, covariates)
  # data = data[, colnames(unique(as.matrix(data), MARGIN = 2))] 
  # 
  ######## Separate the data set into two, one for the treatment group and the other for the control group.
  data1 = data[data$tr == 1, ]
  data0 = data[data$tr == 0, ]
  ######## Nonresponse Weight Estimation in Step 1
  #### Fit multilevel logistic regressions of the response indicator
  ## Model fitted to the treatment group:
  lR1 = glmer(as.formula(paste("R", "~", paste(XR1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  ## Model fitted to the control group:
  lR0 = glmer(as.formula(paste("R", "~", paste(XR0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  #### Predict the response probabilities
  data$pR[data$tr == 1 & data$R == 1] = fitted(lR1)[data1$R == 1]
  data$pR[data$tr == 0 & data$R == 1] = fitted(lR0)[data0$R == 1]
  #### Numerator of the nonresponse weight
  lR.nu1 = glmer(R ~ 1 + (1|site), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  lR.nu0 = glmer(R ~ 1 + (1|site), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  data$pR.nu[data$tr == 1 & data$R == 1] = fitted(lR.nu1)[data1$R == 1]
  data$pR.nu[data$tr == 0 & data$R == 1] = fitted(lR.nu0)[data0$R == 1]  
  #### Construct the nonresponse weight
  data$WR[data$R == 1] = data$pR.nu[data$R == 1]/data$pR[data$R == 1]
  data$WR[data$R == 0] = 1 # Otherwise G is not invertible due to NA's. Because WR will not be used for nonrespondents, this will not affect the final results.
  data$pR[data$R == 0] = 1
  data$pR.nu[data$R == 0] = 1
  
  #### Moment functions for nonresponse weight estimation in step 1
  ## Generate G-H quadrature points and weights
  nnodes = 10
  temp = gauss.quad.prob(nnodes, "normal") # Here, random effects are assumed normal
  nodes = temp$nodes
  weights = temp$weights
  B_1 = t(as.matrix(nodes)) # When there is one random effect, i.e. when the estimated variance of one random effect is 0
  A_1 = weights # When there is one random effect, i.e. when the estimated variance of one random effect is 0
  ## Moment function
  data$R.nu = data$R # This is for the use in fjq, h_pij, h_sigmaj
  h1_fun = function(model, group, data){ # model = "M" for mediator model; model = "R" for response model; group = 0 for control group; group = 1 for treatment group.
    l = get(paste0("l", model, group))
    if(model == "R.nu"){
      x = as.matrix(cbind(rep(1, nrow(data))))
    } else {
      x = as.matrix(cbind(rep(1, nrow(data)), data[, get(paste0("X", model, group))]))   
    }	  
    pi = fixef(l)
    if(round(VarCorr(l)$site[1], 5) == 0){
      sigma = 0
      theta = matrix(0, length(unique(data$site)), 1)
      A = 1
      B = as.matrix(0)
    } else {
      sigma = sqrt(VarCorr(l)$site[1])
      theta = as.matrix(ranef(l)$site)/sigma
      A = A_1
      B = B_1
    }
    p = NULL
    for(q in 1:length(A)){
      p = cbind(p, 1 / (1 + exp(-(x %*% pi + sigma * B[, q]))))
    }
    p = as.matrix(p)
    
    h_pi = NULL
    h_sigma = NULL
    for(j in unique(data$site)){
      data$theta[data$site == j] = theta[as.numeric(rownames(theta)) == j, ]
      lj = 0
      h_pij = 0
      h_sigmaj = 0
      for(q in 1:length(A)){
        pq = p[, q]
        fjq = 1
        for(i in which(data$site == j)){
          fjq = fjq * (pq[i]^data[i, model] * (1 - pq[i])^(1 - data[i, model]))^(data$tr[i] == group)
        }
        lj = lj + fjq * A[q]
        h_pij = h_pij + fjq * A[q] * (data$tr[data$site == j] == group) * (data[data$site == j, model] - pq[data$site == j]) * as.matrix(x[data$site == j,])      
        h_sigmaj =  h_sigmaj + as.matrix(fjq * A[q] * (data$tr[data$site == j] == group) * (data[data$site == j, model] - pq[data$site == j])) * B[, q]
      }
      h_pi = rbind(h_pi, h_pij/lj)
      h_sigma = rbind(h_sigma, h_sigmaj/lj)
    }
    
    if(round(VarCorr(l)$site[1], 5) == 0){
      h1 = h_pi
      V = x
    } else {
      h1 = cbind(h_pi, h_sigma)
      V = cbind(x, data$theta)
    }
    
    return(list(h1 = h1, V = V, sigma = sigma))
  }
  
  h1_R1_list = h1_fun(model = "R", group = 1, data = data) #denominator for treatment group 
  h1_R0_list = h1_fun(model = "R", group = 0, data = data) #denominator for control group
  h1_R1 = h1_R1_list$h1
  h1_R0 = h1_R0_list$h1
  V_R1 = h1_R1_list$V
  V_R0 = h1_R0_list$V
  sigma_R1 = h1_R1_list$sigma
  sigma_R0 = h1_R0_list$sigma
  
  h1_R.nu1_list = h1_fun(model = "R.nu", group = 1, data = data) #numerator for treatment group 
  h1_R.nu0_list = h1_fun(model = "R.nu", group = 0, data = data) #numerator for control group
  h1_R.nu1 = h1_R.nu1_list$h1
  h1_R.nu0 = h1_R.nu0_list$h1
  V_R.nu1 = h1_R.nu1_list$V
  V_R.nu0 = h1_R.nu0_list$V
  sigma_R.nu1 = h1_R.nu1_list$sigma
  sigma_R.nu0 = h1_R.nu0_list$sigma
  
  ######## RMPW Estimation in Step 1
  #### Fit multilevel logistic regressions of the mediator
  ## Sample: respondents in the full sample
  data1 = data[data$tr == 1 & data$R == 1, ]
  data0 = data[data$tr == 0 & data$R == 1, ]
  ## Model fitted to the treatment group:
  lM1 = glmer(as.formula(paste("M", "~", paste(XM1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  ## Model fitted to the control group:
  lM0 = glmer(as.formula(paste("M", "~", paste(XM0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  
  XM1 = colnames(model.matrix(lM1))[-1]
  XM0 = colnames(model.matrix(lM0))[-1]
  
  #### Predict the mediator probabilities
  predict.lmer = function(l, data, X, ranX) {
    nj = NULL
    for (j in as.numeric(rownames(ranef(l)$site))) {
      nj = c(nj, sum(data$site == j))
    }
    ranef = NULL
    for (k in 1:dim(ranef(l)$site)[2]) {
      ranef = cbind(ranef, rep(ranef(l)$site[, k], nj))
    }
    X = cbind(rep(1, nrow(data)), data[, X])
    pred_logit = as.matrix(X) %*% fixef(l) + apply(ranX * ranef, 1, sum)
    pred = exp(pred_logit)/(exp(pred_logit) + 1)
  }
  data$pM1[data$tr == 1 & data$R == 1] = fitted(lM1)
  data$pM0[data$tr == 0 & data$R == 1] = fitted(lM0)
  data$pM1[data$tr == 0 & data$R == 1] = predict.lmer(lM1, data0, XM1, 1) 
  data$pM0[data$tr == 1 & data$R == 1] = predict.lmer(lM0, data1, XM0, 1) 
  
  #### Construct the RMPW weight
  data$rmpw[data$tr == 1 & data$M == 1 & data$R == 1] = data$pM0[data$tr == 1 & data$M == 1 & data$R == 1]/data$pM1[data$tr == 1 & data$M == 1 & data$R == 1]
  data$rmpw[data$tr == 1 & data$M == 0 & data$R == 1] = (1 - data$pM0[data$tr == 1 & data$M == 0 & data$R == 1])/(1 - data$pM1[data$tr == 1 & data$M == 0 & data$R == 1])
  data$rmpw[data$tr == 0 & data$M == 1 & data$R == 1] = data$pM1[data$tr == 0 & data$M == 1 & data$R == 1]/data$pM0[data$tr == 0 & data$M == 1 & data$R == 1]
  data$rmpw[data$tr == 0 & data$M == 0 & data$R == 1] = (1 - data$pM1[data$tr == 0 & data$M == 0 & data$R == 1])/(1 - data$pM0[data$tr == 0 & data$M == 0 & data$R == 1])
  data$rmpw[data$R == 0] = 0
  
  #### Moment functions for RMPW weight estimation in step 1
  h1_M1_list = h1_fun(model = "M", group = 1, data = data[data$R == 1, ]) #Mediator model for treatment group 
  h1_M0_list = h1_fun(model = "M", group = 0, data = data[data$R == 1, ]) #Mediator model for control group
  h1_M1 = matrix(0, nrow(data), ncol(h1_M1_list$h1))
  h1_M0 = matrix(0, nrow(data), ncol(h1_M0_list$h1))
  h1_M1[data$R == 1, ] = h1_M1_list$h1 #Mediator is modeled based on respondents. Therefore, the moment function of the mediator model for nonrespondents is 0.
  h1_M0[data$R == 1, ] = h1_M0_list$h1
  V_M1 = matrix(0, nrow(data), ncol(h1_M1_list$V))
  V_M0 = matrix(0, nrow(data), ncol(h1_M0_list$V))
  V_M1[data$R == 1, ] = h1_M1_list$V
  V_M0[data$R == 1, ] = h1_M0_list$V
  sigma_M1 = h1_M1_list$sigma
  sigma_M0 = h1_M0_list$sigma
  
  h1 = cbind(h1_R0, h1_R1, h1_R.nu0, h1_R.nu1, h1_M0, h1_M1)
  
  ######## Site-Specific Mean Potential Outcome Estimation in Step 2
  N = nrow(data)
  J = length(unique(data$site))
  h2 = matrix(0, N, J * 4)
  mu = NULL
  data$y[data$R == 0] = 0
  data$WD[data$R == 0] = 0
  for (j in unique(data$site)) {
    mu_0 = sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$y)/sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR)
    mu_counter0 = sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw * data$y)/sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw)
    mu_counter1 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw)
    mu_1 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR)
    
    h_0 = data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * (data$y - mu_0)
    h_counter0 = data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw * (data$y - mu_counter0)
    h_counter1 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw * (data$y - mu_counter1)
    h_1 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * (data$y - mu_1)
    
    mu = c(mu, mu_0, mu_counter0, mu_counter1, mu_1)
    h2[, (4 * (which(unique(data$site) == j) - 1) + 1):(4 * which(unique(data$site) == j))] = cbind(h_0, h_counter0, h_counter1, h_1)
  }
  
  if(sum(is.na(mu)) > 0){
    site.omit = which(is.na(mu))[4]/4
    J = J - sum(is.na(mu))/4
    h2 = h2[, -which(is.na(mu))]
    mu = na.omit(mu)
  }
  
  ######## Asymptotic Sampling Variance of the Two-Step Estimators
  #### Stack the moment functions from both steps
  h = cbind(h1, h2)
  H = 1/N * t(as.matrix(h)) %*% as.matrix(h)
  bdiag = function(A, B){
    C = matrix(0, nrow = nrow(A) + nrow(B), ncol = ncol(A) + ncol(B))
    C[1:nrow(A),1:ncol(A)] = A
    C[(1 + nrow(A)):nrow(C),(1 + ncol(A)):ncol(C)] = B
    return(C)
  }
  G11 = as.matrix(bdiag(bdiag(bdiag(-1/N * t(h1_R0) %*% as.matrix(h1_R0), -1/N * t(h1_R1) %*% as.matrix(h1_R1)), bdiag(-1/N * t(h1_R.nu0) %*% as.matrix(h1_R.nu0), -1/N * t(h1_R.nu1) %*% as.matrix(h1_R.nu1))), bdiag(-1/N * t(h1_M0) %*% as.matrix(h1_M0), -1/N * t(h1_M1) %*% as.matrix(h1_M1))))
  G22_diag = NULL
  if(length(unique(data$site)) != J){
    for (j in unique(data$site)[-site.omit]) {
      G22_diag = c(G22_diag, mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR),
                   mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR))
    }
  } else {
    for (j in unique(data$site)) {
      G22_diag = c(G22_diag, mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR),
                   mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR))
    }
  }
  G22 = diag(G22_diag)
  G21 = matrix(0, 4 * J, ncol(G11))
  for(j in 1:J){
    G21[4 * j - 3, 1:ncol(h1_R0)] = apply(h2[, (4 * j - 3)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R0, 2, mean)
    G21[4 * j - 2, 1:ncol(h1_R0)] = apply(h2[, (4 * j - 2)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R0, 2, mean)
    G21[4 * j - 3, (ncol(h1_R0) + ncol(h1_R1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0))] = apply(h2[, (4 * j - 3)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu0, 2, mean)
    G21[4 * j - 2, (ncol(h1_R0) + ncol(h1_R1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0))] = apply(h2[, (4 * j - 2)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu0, 2, mean)
    G21[4 * j - 1, (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (4 * j - 1)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
    G21[4 * j, (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (4 * j)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
    G21[4 * j - 1, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (4 * j - 1)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
    G21[4 * j, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (4 * j)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
    G21[4 * j - 2, (ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + 1):(ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_M0))] = apply((h2[, (4 * j - 2)]/data$rmpw * (-data$M * data$pM1/data$pM0 * (1 - data$pM0) + (1 - data$M) * (1 - data$pM1)/(1 - data$pM0) * data$pM0) * V_M0)[data$R==1, ], 2, sum)/nrow(data)
    G21[4 * j - 2, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1)+ ncol(h1_M0) + 1):ncol(G11)] = apply((h2[, (4 * j - 2)]/data$rmpw * (data$M/data$pM0 - (1 - data$M)/(1 - data$pM0)) * data$pM1 * (1 - data$pM1) * V_M1)[data$R==1, ], 2, sum)/nrow(data)
    G21[4 * j - 1, (ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + 1):(ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_M0))] = apply((h2[, (4 * j - 1)]/data$rmpw * (data$M/data$pM1 - (1 - data$M)/(1 - data$pM1)) * data$pM0 * (1 - data$pM0) * V_M0)[data$R==1, ], 2, sum)/nrow(data)
    G21[4 * j - 1, (ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_M0) + 1):ncol(G11)] = apply((h2[, (4 * j - 1)]/data$rmpw * (-data$M * data$pM0/data$pM1 * (1 - data$pM1) + (1 - data$M) * (1 - data$pM0)/(1 - data$pM1) * data$pM1) * V_M1)[data$R==1, ], 2, sum)/nrow(data)
  }
  G12 = matrix(0, nrow(G11), ncol(G22))
  G= cbind(rbind(G11, G21), rbind(G12, G22))
  V = 1/N * ginv(G) %*% H %*% t(ginv(G))
  V_mu = V[(ncol(V) - 4 * J + 1):ncol(V), (ncol(V) - 4 * J + 1):ncol(V)]  
  
  ######## Population Average Effect Estimation
  Phi = matrix(c(-1,-1,0,-1,1,0,0,0,1,-1,0,1,-1,0,-1,1,0,1,0,1),5, 4)
  beta = kronecker(diag(1, J), Phi) %*% mu
  V_beta = kronecker(diag(1, J), Phi) %*% V_mu %*% kronecker(diag(1, J), t(Phi))
  
  Psi = kronecker(rep(1, J), diag(1, 5))
  gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% beta
  
  ######## Between-Site Variance Estimation
  B = 0
  W = 1/(J * (J - 1)) * t(Psi) %*% V_beta %*% Psi
  for(j in 1:J){
    B = B + 1/(J - 1) * (beta[(5 * j - 4):(5 * j)] - gamma) %*% t((beta[(5 * j - 4):(5 * j)] - gamma))
    W = W - 1/(J - 1) * V_beta[(5 * j - 4):(5 * j), (5 * j - 4):(5 * j)]
  }
  tau = B + W
  
  ######## Standard Error of Population Average Effects
  V_gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% (V_beta + kronecker(diag(1, J), tau)) %*% Psi %*% solve(t(Psi) %*% Psi)
  SE_gamma = sqrt(diag(V_gamma))
  t_gamma = gamma/SE_gamma
  p_gamma = (1 - pnorm(abs(t_gamma))) * 2
  est_gamma = cbind(gamma, SE_gamma, t_gamma, p_gamma)
  est_gamma = as.data.frame(round(est_gamma, 3))

  sig = NULL
  sig[est_gamma[, 4] <= 0.001] = "**"
  sig[est_gamma[, 4] > 0.001 & est_gamma[, 4] <= 0.01] = "*"
  sig[est_gamma[, 4] > 0.01 & est_gamma[, 4] <= 0.05] = "."
  sig[est_gamma[, 4] > 0.05] = ""
  est_gamma = cbind(est_gamma, sig)
  est_gamma[est_gamma[, 4] < 0.001, 4] = "<0.001"
  est_gamma = as.data.frame(est_gamma)
  colnames(est_gamma) = c("Estimate", "Std.Error", "t value", "Pr(>|t|)", "")
  rownames(est_gamma) = c("ITT Effect on Outcome","Natural Direct Effect", "Natural Indirect Effect", "Pure Indirect Effect", "T-M Interaction Effect")
  
  ######## Between-Site Variance Output
  tau.ori = tau
  tau.sub = tau
  
  if("TRUE" %in% (diag(tau.ori) < 0)){
    tau[which(diag(tau < 0)), ] = 0
    tau[, which(diag(tau < 0))] = 0
    tau.sub = tau.ori[-which(diag(tau.ori < 0)), ]
    tau.sub = tau.sub[, -which(diag(tau.sub < 0))]    
  }
  
  cor.sub = cov2cor(tau.sub)
  if(max(na.omit(cor.sub)) > 1)
    cor.sub = suppressWarnings(cor.smooth(cor.sub))
  cor = cor.sub
  
  if("TRUE" %in% (diag(tau.ori) < 0)){
    cor = tau.ori
    cor[which(diag(tau.ori < 0)), ] = NaN
    cor[, which(diag(tau.ori < 0))] = NaN
    cor[which(diag(tau.ori >= 0)), which(diag(tau.ori >= 0))] = cor.sub
  }
  
  var = round(cbind(diag(tau), sqrt(diag(tau)), cor[, 1:4]), 3)
  var[1, 3:6] = ""
  var[2, 4:6] = ""
  var[3, 5:6] = ""
  var[4, 6] = ""
  var = as.data.frame(var)
  colnames(var) = c("Variance", "Std.Dev.", "Corr", "", "", "")
  rownames(var) = c("ITT Effect on Outcome","Natural Direct Effect", "Natural Indirect Effect", "Pure Indirect Effect", "T-M Interaction Effect")
  
  return(list(Random_effects = var, Fixed_effects = est_gamma))
}

#' Variance testing for multisite causal mediation analysis
#'
#' This function performs hypothesis testing for the between-site variance of direct effect and that of indirect effect, besides providing the same output as given by the function msmediate().
#'
#' @param data The data set for analysis.
#' @param y The name of the outcome variable (string).
#' @param treatment The name of the treatment variable (string).
#' @param mediator The name of the mediator variable (string).
#' @param X A vector of variable names (string) of pretreatment covariates, which will be included in the propensity score model. For now, the multilevel propensity score model only allows for one random intercept.
#' @param site The variable name for the site ID (string).
#' @param npermute The number of permutations for the permutation test. The default value is 200. It may take a long time, depending on the sample size and the length of X.
#' @return A list contains the hypothesis testing results of the between-site variance of the causal effects, besides the same output as given by the function msmediate().
#' @author Xu Qin and Guanglei Hong
#' @references Qin, X., & Hong, G (2017). A weighting method for assessing between-site heterogeneity in causal mediation mechanism. Journal of Educational and Behavioral Statistics. Journal of Educational and Behavioral Statistics. Journal of Educational and Behavioral Statistics, 42(3), 308-340. \doi{10.3102/1076998617694879}
#' @export
#' @importFrom stats as.formula binomial coef fitted lm model.matrix pnorm
#' @importFrom lme4 VarCorr fixef glmer ranef
#' @importFrom statmod gauss.quad.prob
#' @examples 
#' data(sim)
#'
#' vartest.msmediate(data = sim, y = "y", treatment = "tr", mediator = "me", X = c("x1", "x2", "x3", "x4"), site = "site", npermute = 2)
#'
vartest.msmediate = function(data, y, treatment, mediator, X, site, npermute = 200) {
  data = as.data.frame(data)
  colnames(data)[which(colnames(data) == treatment)] = "tr"
  colnames(data)[which(colnames(data) == mediator)] = "me"
  colnames(data)[which(colnames(data) == site)] = "site"
  colnames(data)[which(colnames(data) == y)] = "y"
  
  # # Factorize categorical covariates (with fewer than 10 categories)
  # for(i in 1:length(X)){
  #   if(length(unique(data[, X[i]])) > 2 & length(unique(data[, X[i]])) < 10){
  #     data[, X[i]] = as.factor(data[, X[i]])
  #   }
  # }
  # covariates = model.matrix(as.formula(paste("~", paste(X, collapse = "+"))), data)[, -1]
  # data = data[, -which(colnames(data) %in% X)]
  # data = cbind(data, covariates)
  # X = colnames(covariates)
  
  est = function(data, y, treatment, mediator, X, site, permutation = F){
    ranX = 1
    data = data[order(data$site), ]
    data1 = data[data$tr == 1, ]
    data0 = data[data$tr == 0, ]
    
    l_tr = suppressWarnings({glmer(as.formula(paste("me", "~", paste(X, collapse="+"), "+(", ranX, "|site)")), data = data1, family = binomial, nAGQ = 10)})
    data$p1[data$tr == 1] = fitted(l_tr)
    l_ctrl = suppressWarnings({glmer(as.formula(paste("me", "~", paste(X, collapse="+"), "+(", ranX, "|site)")), data = data0, family = binomial, nAGQ = 10)})
    data$p0[data$tr == 0] = fitted(l_ctrl)
    
    predict.lmer = function(l, data, X, ranX) {
      nj = NULL
      for (j in as.numeric(rownames(ranef(l)$site))) {
        nj = c(nj, sum(data$site == j))
      }
      ranef = NULL
      for (k in 1:dim(ranef(l)$site)[2]) {
        ranef = cbind(ranef, rep(ranef(l)$site[, k], nj))
      }
      pred_logit = X %*% fixef(l) + apply(ranX * ranef, 1, sum)
      pred = exp(pred_logit)/(exp(pred_logit) + 1)
    }
    data$p0[data$tr == 1] = predict.lmer(l_ctrl, data1, model.matrix(l_tr), rep(1, nrow(data1)))
    data$p1[data$tr == 0] = predict.lmer(l_tr, data0, model.matrix(l_ctrl), rep(1, nrow(data0)))
    
    data$rmpw[data$tr == 1 & data$me == 1] = (data$p0/data$p1)[data$tr == 1 & data$me == 1]
    data$rmpw[data$tr == 1 & data$me == 0] = ((1 - data$p0)/(1 - data$p1))[data$tr == 1 & data$me == 0]
    data$rmpw[data$tr == 0] = 1
    
    data = data[order(data$site),]
    data0 = data[data$tr == 0, ]
    data1 = data[data$tr == 1, ]
    
    f = function(x, alpha, s, F, B){
      return(1 / (1 + exp(-(x %*% alpha + s * F * B)))) # When the estimated variance of both two or one of the two random effects is 1
    }
    
    nnodes = 10
    temp = gauss.quad.prob(nnodes, "normal") # Here, random effects are assumed normal
    nodes = temp$nodes
    weights = temp$weights
    # When there is one random effect, i.e. when the estimated variance of one random effect is 0
    B_1 = t(as.matrix(nodes))
    A_1 = weights
    
    x = model.matrix(as.formula(paste("me", "~", paste(X, collapse="+"))), data = data) # Here we have the same covariates in the models for both the treatment and control group
    h1_fun = function(data0, data){
      if(data0$tr[1] == 0){
        l0 = l_ctrl
      } else{
        l0 = l_tr
      }        
      alpha0 = fixef(l0)
      if(round(VarCorr(l0)$site[1], 5) == 0){
        F0 = 0
        v0 = 0
        s0 = matrix(0, nrow(x), 1)
        J0 = 0
        A0 = 1
        B0 = as.matrix(0)
      } else {
        F0 = sqrt(VarCorr(l0)$site[1])
        v0 = ranef(l0)$site
        s0 = matrix(1, nrow(x), 1)
        J0 = 1
        A0 = A_1
        B0 = B_1
      }
      p0 = NULL
      for(q in 1:length(A0)){
        p0 = cbind(p0, f(x = x, alpha = alpha0, s = s0, F = F0, B = B0[, q]))
      }
      p0 = as.matrix(p0)
      
      h_alpha0 = NULL
      h_F0 = NULL
      for(j in unique(data$site)){
        l0j = 0
        h_alpha0j = 0
        h_F0j = 0
        for(q in 1:length(A0)){
          p0q = p0[, q]
          if(data0$tr[1] == 0){
            f0jq = 1
            for(i in which(data$site == j)){
              f0jq = f0jq * (p0q[i]^data$me[i] * (1 - p0q[i])^(1 - data$me[i]))^(1 - data$tr[i])
            }
            l0j = l0j + f0jq * A0[q]
            h_alpha0j = h_alpha0j + as.numeric(f0jq * A0[q] * (1 - data$tr[data$site == j]) * (data$me[data$site == j] - p0q[data$site == j])) * x[data$site == j,]      
            h_F0j =  h_F0j + as.numeric(f0jq * A0[q] * (1 - data$tr[data$site == j]) * (data$me[data$site == j] - p0q[data$site == j])) * kronecker(t(B0[, q]), s0[data$site == j,]) %*% t(J0)    
          } else {
            f0jq = 1
            for(i in which(data$site == j)){
              f0jq = f0jq * (p0q[i]^data$me[i] * (1 - p0q[i])^(1 - data$me[i]))^(data$tr[i])
            }
            l0j = l0j + f0jq * A0[q]
            h_alpha0j = h_alpha0j + as.numeric(f0jq * A0[q] * data$tr[data$site == j] * (data$me[data$site == j] - p0q[data$site == j])) * x[data$site == j,]              
            h_F0j =  h_F0j + as.numeric(f0jq * A0[q] * data$tr[data$site == j] * (data$me[data$site == j] - p0q[data$site == j])) * kronecker(t(B0[, q]), s0[data$site == j,]) %*% t(J0)    
          }        
        }
        h_alpha0 = rbind(h_alpha0, h_alpha0j/l0j)
        h_F0 = rbind(h_F0, h_F0j/l0j)
      }
      
      if(round(VarCorr(l0)$site[1], 5) == 0){
        h1_0 = h_alpha0
      } else {
        h1_0 = cbind(h_alpha0, h_F0)
      }
      
      if(round(VarCorr(l0)$site[1], 5) == 0){
        theta0 = matrix(0, 1, length(unique(data0$site)))
      } else {
        theta0 = solve(F0) %*% t(v0)
      }
      
      return(list(h1 = h1_0, theta = theta0, s = s0, J = J0))
    }
    
    h1_list0 = h1_fun(data0, data)
    h1_list1 = h1_fun(data1, data)
    h1_0 = h1_list0$h1 
    h1_1 = h1_list1$h1
    h1 = cbind(h1_0, h1_1)
    p = ncol(h1)
    theta0 = h1_list0$theta
    theta1 = h1_list1$theta
    s0 = h1_list0$s
    s1 = h1_list1$s
    J0 = h1_list0$J
    J1 = h1_list1$J
    
    N = nrow(data)
    J = length(unique(data$site))
    h2 = matrix(0, N, J * 3)
    mu = NULL
    nj = NULL
    for (j in unique(data$site)) {
      dataj = data[data$site == j, ]
      nj = c(nj, nrow(dataj))
      mu_0 = mean(dataj$y[dataj$tr == 0])
      mu_counter = sum(dataj$y * dataj$rmpw * dataj$tr)/sum(dataj$rmpw * dataj$tr)
      mu_1 = mean(dataj$y[dataj$tr == 1])
      h_0 = (data$y - mu_0) * (1 - data$tr) * (data$site == j)
      h_counter = (data$y - mu_counter) * data$rmpw * data$tr * (data$site == j)
      h_1 = (data$y - mu_1) * data$tr * (data$site == j)
      h2[, (3 * (which(unique(data$site) == j) - 1) + 1):(3 * which(unique(data$site) == j))] = cbind(h_0, 
                                                                                                      h_counter, h_1)
      mu = c(mu, mu_0, mu_counter, mu_1)
    }
    
    h = cbind(h1, h2)  # Stack the moment functions from both steps
    H = 1/N * t(as.matrix(h)) %*% as.matrix(h)  # See Equation S13
    ## The derivation of R (See Equation S14)
    bdiag = function(A, B){
      C = matrix(0, nrow = nrow(A) + nrow(B), ncol = ncol(A) + ncol(B))
      C[1:nrow(A),1:ncol(A)] = A
      C[(1 + nrow(A)):nrow(C),(1 + ncol(A)):ncol(C)] = B
      return(C)
    }
    R11 = as.matrix(bdiag(-1/N * t(h1_0) %*% as.matrix(h1_0), -1/N * t(h1_1) %*% as.matrix(h1_1)))
    
    R22_diag = NULL
    for (j in unique(data$site)) {
      R22_diag = c(R22_diag, mean(-(1 - data$tr) * (data$site == j)), mean(-data$rmpw * data$tr * (data$site == j)), mean(-data$tr * (data$site == j)))
    }
    R22 = diag(R22_diag)
    
    p0 = data$p0
    p1 = data$p1
    R21 = NULL
    for (j in 1:J) {
      dh_counter_dalpha0 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (data$me/p1 - (1 - data$me)/(1 - p1)) *   
                                   p0 * (1 - p0) * x, 2, mean)
      dh_counter_dF0 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (data$me/p1 - (1 - data$me)/(1 - p1)) * p0 
                             * (1 - p0) * kronecker(t(theta0[, j]), s0) %*% t(J0), 2, mean)
      dh_counter_dalpha1 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (-data$me * p0/p1 * (1 - p1) + (1 - data$me) 
                                                                    * (1 - p0)/(1 - p1) * p1) * x, 2, mean)
      dh_counter_dF1 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (-data$me * p0/p1 * (1 - p1) + (1 - data$me) 
                                                                * (1 - p0)/(1 - p1) * p1) * kronecker(t(theta1[, j]), s1) %*% t(J1), 2, mean)
      if (ncol(h1_0) - ncol(x) == 0 & ncol(h1_1) - ncol(x) == 0) {
        R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dalpha1), rep(0, p)))
      } else if (ncol(h1_0) - ncol(x) == 0) {
        R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dalpha1, dh_counter_dF1), 
                               rep(0, p)))
      } else if (ncol(h1_1) - ncol(x) == 0) {
        R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dF0, dh_counter_dalpha1), 
                               rep(0, p)))
      } else {
        R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dF0, dh_counter_dalpha1, 
                                            dh_counter_dF1), rep(0, p)))
      }
    }
    
    R12 = matrix(0, p, 3 * J)
    
    R = matrix(0, p + 3 * J, p + 3 * J)
    R[1:p, 1:p] = R11
    R[1:p, (p + 1):(p + 3 * J)] = R12
    R[(p + 1):(p + 3 * J), 1:p] = R21
    R[(p + 1):(p + 3 * J), (p + 1):(p + 3 * J)] = R22
    
    V_all = 1/N * solve(R) %*% H %*% t(solve(R))  # See Appendix A
    V_mu = V_all[(ncol(V_all) - 3 * J + 1):ncol(V_all), (ncol(V_all) - 3 * J + 1):ncol(V_all)]  
    
    Phi = matrix(c(-1, 0, 1, -1, 0, 1), 2, 3)
    for (j in 1:(J - 1)) {
      Phi = as.matrix(bdiag(Phi, matrix(c(-1, 0, 1, -1, 0, 1), 2, 3)))
    }
    beta = Phi %*% mu # Estimated site-specific direct and indirect effects
    V_beta = Phi %*% V_mu %*% t(Phi)
    dej = beta[seq(1, 2 * J, 2)]
    iej = beta[seq(2, 2 * J, 2)]
    
    Psi = cbind(rep(c(1, 0), J), rep(c(0, 1), J))
    gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% beta
    de = gamma[1]
    ie = gamma[2]
    
    G = 0
    V = 0
    for(j in 1:J){
      G = G + (beta[(2 * j - 1): (2 * j)] - gamma) %*% t(beta[(2 * j - 1): (2 * j)] - gamma)    
      Vj = V_beta[(2 * j - 1):(2 * j), (2 * j - 1):(2 * j)]
      V = V + Vj
    }
    tau = 1/(J - 1) * (G - V + 1/J * t(Psi) %*% V_beta %*% Psi)
    
    if(permutation == F){
      V_gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% (V_beta + kronecker(diag(rep(1,J)), tau)) %*% Psi %*%   
        solve(t(Psi) %*% Psi)
      V_de = V_gamma[1, 1]
      V_ie = V_gamma[2, 2]   
      
      se_de = sqrt(V_de)
      se_ie = sqrt(V_ie)
      
      z_de = de/se_de
      z_ie = ie/se_ie
      
      p_de = (1 - pnorm(abs(z_de))) * 2
      p_ie = (1 - pnorm(abs(z_ie))) * 2
      
      est = cbind(round(c(de, ie), 5), round(c(se_de, se_ie), 5), round(c(z_de, z_ie), 3), round(c(p_de, p_ie), 3))
      
      sig = NULL
      sig[est[, 4] <= 0.001] = "**"
      sig[est[, 4] > 0.001 & est[, 4] <= 0.01] = "*"
      sig[est[, 4] > 0.01 & est[, 4] <= 0.05] = "."
      sig[est[, 4] > 0.05] = ""
      
      est = cbind(est, sig)
      est[est[, 4] < 0.001, 4] = "<0.001"
      est = as.data.frame(est)
      colnames(est) = c("Estimate", "Std.Error", "t value", "Pr(>|t|)", "")
      rownames(est) = c("Natrual Direct Effect", "Natrual Indirect Effect")
      
      tau.ori = tau
      tau.sub = tau
      
      if("TRUE" %in% (diag(tau.ori) < 0)){
        tau[which(diag(tau < 0)), ] = 0
        tau[, which(diag(tau < 0))] = 0
        tau.sub = tau.ori[-which(diag(tau.ori < 0)), ]
        tau.sub = tau.sub[, -which(diag(tau.sub < 0))]    
      }
      
      cor.sub = cov2cor(tau.sub)
      if(max(na.omit(cor.sub)) > 1)
        cor.sub = suppressWarnings(cor.smooth(cor.sub))
      cor = cor.sub
      
      if("TRUE" %in% (diag(tau.ori) < 0)){
        cor = tau.ori
        cor[which(diag(tau.ori < 0)), ] = NaN
        cor[, which(diag(tau.ori < 0))] = NaN
        cor[which(diag(tau.ori >= 0)), which(diag(tau.ori >= 0))] = cor.sub
      }
      
      cor_de_ie = cor[2, 1]
      
      var = cbind(round(c(tau[1, 1], tau[2, 2]), 3), round(c(sqrt(tau[1, 1]), sqrt(tau[2, 2])), 3), c("", round(cor_de_ie, 3)))
      var = as.data.frame(var)
      colnames(var) = c("Variance", "Std.Dev.", "Corr")
      rownames(var) = c("Natrual Direct Effect", "Natrual Indirect Effect")
    }
    
    chisq_de = 0
    chisq_ie = 0
    for (j in 1:J) {
      Vj = V_beta[(2 * (j - 1) + 1):(2 * j), (2 * (j - 1) + 1):(2 * j)]
      chisq_de = chisq_de + (dej[j] - de)^2/Vj[1, 1]
      chisq_ie = chisq_ie + (iej[j] - ie)^2/Vj[2, 2]
    }
    
    if(permutation == F)
      result = list(Random_effects = var, Fixed_effects = est, chisq = c(chisq_de = chisq_de, chisq_ie = chisq_ie))
    
    if(permutation == T)
      result = list(chisq = c(chisq_de = chisq_de, chisq_ie = chisq_ie))
    
    return(result)
  }
  
  result = suppressWarnings({est(data, y, treatment, mediator, X, site, permutation = F)})
  
  siteID = data$site
  chisq_de_pm = NULL
  chisq_ie_pm = NULL
  for (i in 1:npermute) {
    try({
      data$site = sample(siteID, length(siteID), replace = F)
      est_list_pm = suppressWarnings({est(data, y, treatment, mediator, X, site, permutation = T)$chisq})
      chisq_de_pm = c(chisq_de_pm, est_list_pm["chisq_de"])
      chisq_ie_pm = c(chisq_ie_pm, est_list_pm["chisq_ie"])
    }, silent = T)
  }
  
  pvalue_var_de = sum(chisq_de_pm >= result$chisq["chisq_de"])/length(chisq_de_pm)
  pvalue_var_ie = sum(chisq_ie_pm >= result$chisq["chisq_ie"])/length(chisq_ie_pm)
  
  pvalue = round(c(pvalue_var_de, pvalue_var_ie), 5)
  
  sig = NULL
  sig[pvalue <= 0.001] = "**"
  sig[pvalue > 0.001 & pvalue <= 0.01] = "*"
  sig[pvalue > 0.01 & pvalue <= 0.05] = "."
  sig[pvalue > 0.05] = ""
  result$Random_effects = cbind(result$Random_effects, pvalue, sig)
  colnames(result$Random_effects)[4] = "p-value"
  colnames(result$Random_effects)[5] = ""
  
  return(list(Random_effects = result$Random_effects, Fixed_effects = result$Fixed_effects))
}

#' Variance testing for multisite causal mediation analysis in the presence of complex sample and survey designs and non-random nonresponse
#'
#' This function performs hypothesis testing for the between-site variances of natural direct effect, natural indirect effect, pure indirect effect, and treatment-by-mediator interaction effect in the presence of complex sample and survey designs and non-random nonresponse, besides providing the same output as given by the function msmediate.weights().
#' 
#' @param data The data set for analysis.
#' @param y The name of the outcome variable (string).
#' @param treatment The name of the treatment variable (string).
#' @param mediator The name of the mediator variable (string).
#' @param response The name of the response variable (string), which is equal to 1 if the individual responded and 0 otherwise.
#' @param XR1 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XR0 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM1 A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM0 A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param site The variable name for the site ID (string).
#' @param sample.weight The variable name for the sample weight given by design (string).
#' @param npermute The number of permutations for the permutation test. The default value is 200. It may take a long time, depending on the sample size and the length of X.
#' @return A list contains the hypothesis testing results of the between-site variance of the causal effects, besides the same output as given by the function msmediate().
#' @author Xu Qin, Guanglei Hong, Jonah Deutsch, and Edward Bein
#' @references Qin, X., Hong, G., Deutsch, J., & Bein, E. (2019). Multisite causal mediation analysis in the presence of complex sample and survey designs and non-random non-response. Journal of the Royal Statistical Society: Series A (Statistics in Society), 182(4), 1343-1370. \doi{10.1111/rssa.12446}
#' @export
#' @importFrom stats as.formula binomial coef fitted lm model.matrix pnorm na.omit cov2cor
#' @importFrom lme4 VarCorr fixef glmer ranef glmerControl
#' @importFrom statmod gauss.quad.prob
#' @importFrom psych cor.smooth
#' @importFrom MASS ginv
#' @examples 
#' data(sim.weights)
#'
#' vartest.msmediate.weights(data = sim.weights, y = "y", treatment = "tr", mediator = "me", response = "R", XR1 = c("x1", "x2", "x3"), XR0 = c("x1", "x2", "x3"), XM1 = c("x1", "x2", "x3"), XM0 = c("x1", "x2", "x3"), site = "site", sample.weight = "WD", npermute = 2)
#'
vartest.msmediate.weights = function(data, y, treatment, mediator, response, XR1, XR0, XM1, XM0, site, sample.weight, npermute = 200) {
  data = as.data.frame(data)
  colnames(data)[which(colnames(data) == y)] = "y"
  colnames(data)[which(colnames(data) == treatment)] = "tr"
  colnames(data)[which(colnames(data) == mediator)] = "M"
  colnames(data)[which(colnames(data) == response)] = "R"
  colnames(data)[which(colnames(data) == site)] = "site"
  colnames(data)[which(colnames(data) == sample.weight)] = "WD"
  
  # # Factorize categorical covariates (with fewer than 10 categories)
  # transform = function(X){
  #   for(i in 1:length(X)){
  #     if(length(unique(data[, X[i]])) > 2 & length(unique(data[, X[i]])) < 10){
  #       data[, X[i]] = as.factor(data[, X[i]])
  #     }
  #   }
  #   covariates = model.matrix(as.formula(paste("~", paste(X, collapse = "+"))), data)
  #   X = colnames(covariates)
  #   return(list(covariates = covariates[, -1], X = X[-1]))
  # }
  # transform.XR1 = transform(XR1)
  # transform.XR0 = transform(XR0)
  # transform.XM1 = transform(XM1)
  # transform.XM0 = transform(XM0)
  # data = data[, -which(colnames(data) %in% unique(c(XR1, XR0, XM1, XM0)))]
  # XR1 = transform.XR1$X
  # XR0 = transform.XR0$X
  # XM1 = transform.XM1$X
  # XM0 = transform.XM0$X
  # covariates = cbind(transform.XR1$covariates, transform.XR0$covariates, transform.XM1$covariates, transform.XM0$covariates)
  # colnames(covariates) = c(XR1, XR0, XM1, XM0)
  # data = cbind(data, covariates)
  # data = data[, colnames(unique(as.matrix(data), MARGIN = 2))] 
   
  est = function(data, y, treatment, mediator, response, XR1, XR0, XM1, XM0, site, sample.weight, permutation = F){
    data = data[order(data$site), ]
    ######## Separate the data set into two, one for the treatment group and the other for the control group.
    data1 = data[data$tr == 1, ]
    data0 = data[data$tr == 0, ]
    ######## Nonresponse Weight Estimation in Step 1
    #### Fit multilevel logistic regressions of the response indicator
    ## Model fitted to the treatment group:
    lR1 = glmer(as.formula(paste("R", "~", paste(XR1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    ## Model fitted to the control group:
    lR0 = glmer(as.formula(paste("R", "~", paste(XR0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    #### Predict the response probabilities
    data$pR[data$tr == 1 & data$R == 1] = fitted(lR1)[data1$R == 1]
    data$pR[data$tr == 0 & data$R == 1] = fitted(lR0)[data0$R == 1]
    #### Numerator of the nonresponse weight
    lR.nu1 = glmer(R ~ 1 + (1|site), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    lR.nu0 = glmer(R ~ 1 + (1|site), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    data$pR.nu[data$tr == 1 & data$R == 1] = fitted(lR.nu1)[data1$R == 1]
    data$pR.nu[data$tr == 0 & data$R == 1] = fitted(lR.nu0)[data0$R == 1]  
    #### Construct the nonresponse weight
    data$WR[data$R == 1] = data$pR.nu[data$R == 1]/data$pR[data$R == 1]
    data$WR[data$R == 0] = 1 # Otherwise G is not invertible due to NA's. Because WR will not be used for nonrespondents, this will not affect the final results.
    data$pR[data$R == 0] = 1
    data$pR.nu[data$R == 0] = 1
    
    #### Moment functions for nonresponse weight estimation in step 1
    ## Generate G-H quadrature points and weights
    nnodes = 10
    temp = gauss.quad.prob(nnodes, "normal") # Here, random effects are assumed normal
    nodes = temp$nodes
    weights = temp$weights
    B_1 = t(as.matrix(nodes)) # When there is one random effect, i.e. when the estimated variance of one random effect is 0
    A_1 = weights # When there is one random effect, i.e. when the estimated variance of one random effect is 0
    ## Moment function
    data$R.nu = data$R # This is for the use in fjq, h_pij, h_sigmaj
    h1_fun = function(model, group, data){ # model = "M" for mediator model; model = "R" for response model; group = 0 for control group; group = 1 for treatment group.
      l = get(paste0("l", model, group))
      if(model == "R.nu"){
        x = as.matrix(cbind(rep(1, nrow(data))))
      } else {
        x = as.matrix(cbind(rep(1, nrow(data)), data[, get(paste0("X", model, group))]))   
      }	  
      pi = fixef(l)
      if(round(VarCorr(l)$site[1], 5) == 0){
        sigma = 0
        theta = matrix(0, length(unique(data$site)), 1)
        A = 1
        B = as.matrix(0)
      } else {
        sigma = sqrt(VarCorr(l)$site[1])
        theta = as.matrix(ranef(l)$site)/sigma
        A = A_1
        B = B_1
      }
      p = NULL
      for(q in 1:length(A)){
        p = cbind(p, 1 / (1 + exp(-(x %*% pi + sigma * B[, q]))))
      }
      p = as.matrix(p)
      
      h_pi = NULL
      h_sigma = NULL
      for(j in unique(data$site)){
        data$theta[data$site == j] = theta[as.numeric(rownames(theta)) == j, ]
        lj = 0
        h_pij = 0
        h_sigmaj = 0
        for(q in 1:length(A)){
          pq = p[, q]
          fjq = 1
          for(i in which(data$site == j)){
            fjq = fjq * (pq[i]^data[i, model] * (1 - pq[i])^(1 - data[i, model]))^(data$tr[i] == group)
          }
          lj = lj + fjq * A[q]
          h_pij = h_pij + fjq * A[q] * (data$tr[data$site == j] == group) * (data[data$site == j, model] - pq[data$site == j]) * as.matrix(x[data$site == j,])      
          h_sigmaj =  h_sigmaj + as.matrix(fjq * A[q] * (data$tr[data$site == j] == group) * (data[data$site == j, model] - pq[data$site == j])) * B[, q]
        }
        h_pi = rbind(h_pi, h_pij/lj)
        h_sigma = rbind(h_sigma, h_sigmaj/lj)
      }
      
      if(round(VarCorr(l)$site[1], 5) == 0){
        h1 = h_pi
        V = x
      } else {
        h1 = cbind(h_pi, h_sigma)
        V = cbind(x, data$theta)
      }
      
      return(list(h1 = h1, V = V, sigma = sigma))
    }
    
    h1_R1_list = h1_fun(model = "R", group = 1, data = data) #denominator for treatment group 
    h1_R0_list = h1_fun(model = "R", group = 0, data = data) #denominator for control group
    h1_R1 = h1_R1_list$h1
    h1_R0 = h1_R0_list$h1
    V_R1 = h1_R1_list$V
    V_R0 = h1_R0_list$V
    sigma_R1 = h1_R1_list$sigma
    sigma_R0 = h1_R0_list$sigma
    
    h1_R.nu1_list = h1_fun(model = "R.nu", group = 1, data = data) #numerator for treatment group 
    h1_R.nu0_list = h1_fun(model = "R.nu", group = 0, data = data) #numerator for control group
    h1_R.nu1 = h1_R.nu1_list$h1
    h1_R.nu0 = h1_R.nu0_list$h1
    V_R.nu1 = h1_R.nu1_list$V
    V_R.nu0 = h1_R.nu0_list$V
    sigma_R.nu1 = h1_R.nu1_list$sigma
    sigma_R.nu0 = h1_R.nu0_list$sigma
    
    ######## RMPW Estimation in Step 1
    #### Fit multilevel logistic regressions of the mediator
    ## Sample: respondents in the full sample
    data1 = data[data$tr == 1 & data$R == 1, ]
    data0 = data[data$tr == 0 & data$R == 1, ]
    ## Model fitted to the treatment group:
    lM1 = glmer(as.formula(paste("M", "~", paste(XM1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    ## Model fitted to the control group:
    lM0 = glmer(as.formula(paste("M", "~", paste(XM0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    
    XM1 = colnames(model.matrix(lM1))[-1]
    XM0 = colnames(model.matrix(lM0))[-1]
    
    #### Predict the mediator probabilities
    predict.lmer = function(l, data, X, ranX) {
      nj = NULL
      for (j in as.numeric(rownames(ranef(l)$site))) {
        nj = c(nj, sum(data$site == j))
      }
      ranef = NULL
      for (k in 1:dim(ranef(l)$site)[2]) {
        ranef = cbind(ranef, rep(ranef(l)$site[, k], nj))
      }
      X = cbind(rep(1, nrow(data)), data[, X])
      pred_logit = as.matrix(X) %*% fixef(l) + apply(ranX * ranef, 1, sum)
      pred = exp(pred_logit)/(exp(pred_logit) + 1)
    }
    data$pM1[data$tr == 1 & data$R == 1] = fitted(lM1)
    data$pM0[data$tr == 0 & data$R == 1] = fitted(lM0)
    data$pM1[data$tr == 0 & data$R == 1] = predict.lmer(lM1, data0, XM1, 1) 
    data$pM0[data$tr == 1 & data$R == 1] = predict.lmer(lM0, data1, XM0, 1) 
    
    #### Construct the RMPW weight
    data$rmpw[data$tr == 1 & data$M == 1 & data$R == 1] = data$pM0[data$tr == 1 & data$M == 1 & data$R == 1]/data$pM1[data$tr == 1 & data$M == 1 & data$R == 1]
    data$rmpw[data$tr == 1 & data$M == 0 & data$R == 1] = (1 - data$pM0[data$tr == 1 & data$M == 0 & data$R == 1])/(1 - data$pM1[data$tr == 1 & data$M == 0 & data$R == 1])
    data$rmpw[data$tr == 0 & data$M == 1 & data$R == 1] = data$pM1[data$tr == 0 & data$M == 1 & data$R == 1]/data$pM0[data$tr == 0 & data$M == 1 & data$R == 1]
    data$rmpw[data$tr == 0 & data$M == 0 & data$R == 1] = (1 - data$pM1[data$tr == 0 & data$M == 0 & data$R == 1])/(1 - data$pM0[data$tr == 0 & data$M == 0 & data$R == 1])
    data$rmpw[data$R == 0] = 0
    
    #### Moment functions for RMPW weight estimation in step 1
    h1_M1_list = h1_fun(model = "M", group = 1, data = data[data$R == 1, ]) #Mediator model for treatment group 
    h1_M0_list = h1_fun(model = "M", group = 0, data = data[data$R == 1, ]) #Mediator model for control group
    h1_M1 = matrix(0, nrow(data), ncol(h1_M1_list$h1))
    h1_M0 = matrix(0, nrow(data), ncol(h1_M0_list$h1))
    h1_M1[data$R == 1, ] = h1_M1_list$h1 #Mediator is modeled based on respondents. Therefore, the moment function of the mediator model for nonrespondents is 0.
    h1_M0[data$R == 1, ] = h1_M0_list$h1
    V_M1 = matrix(0, nrow(data), ncol(h1_M1_list$V))
    V_M0 = matrix(0, nrow(data), ncol(h1_M0_list$V))
    V_M1[data$R == 1, ] = h1_M1_list$V
    V_M0[data$R == 1, ] = h1_M0_list$V
    sigma_M1 = h1_M1_list$sigma
    sigma_M0 = h1_M0_list$sigma
    
    h1 = cbind(h1_R0, h1_R1, h1_R.nu0, h1_R.nu1, h1_M0, h1_M1)
    
    ######## Site-Specific Mean Potential Outcome Estimation in Step 2
    N = nrow(data)
    J = length(unique(data$site))
    h2 = matrix(0, N, J * 4)
    mu = NULL
    data$y[data$R == 0] = 0
    data$WD[data$R == 0] = 0
    for (j in unique(data$site)) {
      mu_0 = sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$y)/sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR)
      mu_counter0 = sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw * data$y)/sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw)
      mu_counter1 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw)
      mu_1 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR)
      
      h_0 = data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * (data$y - mu_0)
      h_counter0 = data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw * (data$y - mu_counter0)
      h_counter1 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw * (data$y - mu_counter1)
      h_1 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * (data$y - mu_1)
      
      mu = c(mu, mu_0, mu_counter0, mu_counter1, mu_1)
      h2[, (4 * (which(unique(data$site) == j) - 1) + 1):(4 * which(unique(data$site) == j))] = cbind(h_0, h_counter0, h_counter1, h_1)
    }
    
    if(sum(is.na(mu)) > 0){
      site.omit = which(is.na(mu))[4]/4
      J = J - sum(is.na(mu))/4
      h2 = h2[, -which(is.na(mu))]
      mu = na.omit(mu)
    }
    
    ######## Asymptotic Sampling Variance of the Two-Step Estimators
    #### Stack the moment functions from both steps
    h = cbind(h1, h2)
    H = 1/N * t(as.matrix(h)) %*% as.matrix(h)
    bdiag = function(A, B){
      C = matrix(0, nrow = nrow(A) + nrow(B), ncol = ncol(A) + ncol(B))
      C[1:nrow(A),1:ncol(A)] = A
      C[(1 + nrow(A)):nrow(C),(1 + ncol(A)):ncol(C)] = B
      return(C)
    }
    G11 = as.matrix(bdiag(bdiag(bdiag(-1/N * t(h1_R0) %*% as.matrix(h1_R0), -1/N * t(h1_R1) %*% as.matrix(h1_R1)), bdiag(-1/N * t(h1_R.nu0) %*% as.matrix(h1_R.nu0), -1/N * t(h1_R.nu1) %*% as.matrix(h1_R.nu1))), bdiag(-1/N * t(h1_M0) %*% as.matrix(h1_M0), -1/N * t(h1_M1) %*% as.matrix(h1_M1))))
    G22_diag = NULL
    if(length(unique(data$site)) != J){
      for (j in unique(data$site)[-site.omit]) {
        G22_diag = c(G22_diag, mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR),
                     mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR))
      }
    } else {
      for (j in unique(data$site)) {
        G22_diag = c(G22_diag, mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR),
                     mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$rmpw),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR))
      }
    }
    G22 = diag(G22_diag)
    G21 = matrix(0, 4 * J, ncol(G11))
    for(j in 1:J){
      G21[4 * j - 3, 1:ncol(h1_R0)] = apply(h2[, (4 * j - 3)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R0, 2, mean)
      G21[4 * j - 2, 1:ncol(h1_R0)] = apply(h2[, (4 * j - 2)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R0, 2, mean)
      G21[4 * j - 3, (ncol(h1_R0) + ncol(h1_R1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0))] = apply(h2[, (4 * j - 3)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu0, 2, mean)
      G21[4 * j - 2, (ncol(h1_R0) + ncol(h1_R1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0))] = apply(h2[, (4 * j - 2)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu0, 2, mean)
      G21[4 * j - 1, (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (4 * j - 1)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
      G21[4 * j, (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (4 * j)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
      G21[4 * j - 1, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (4 * j - 1)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
      G21[4 * j, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (4 * j)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
      G21[4 * j - 2, (ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + 1):(ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_M0))] = apply((h2[, (4 * j - 2)]/data$rmpw * (-data$M * data$pM1/data$pM0 * (1 - data$pM0) + (1 - data$M) * (1 - data$pM1)/(1 - data$pM0) * data$pM0) * V_M0)[data$R==1, ], 2, sum)/nrow(data)
      G21[4 * j - 2, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1)+ ncol(h1_M0) + 1):ncol(G11)] = apply((h2[, (4 * j - 2)]/data$rmpw * (data$M/data$pM0 - (1 - data$M)/(1 - data$pM0)) * data$pM1 * (1 - data$pM1) * V_M1)[data$R==1, ], 2, sum)/nrow(data)
      G21[4 * j - 1, (ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + 1):(ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_M0))] = apply((h2[, (4 * j - 1)]/data$rmpw * (data$M/data$pM1 - (1 - data$M)/(1 - data$pM1)) * data$pM0 * (1 - data$pM0) * V_M0)[data$R==1, ], 2, sum)/nrow(data)
      G21[4 * j - 1, (ncol(h1_R0) + ncol(h1_R1)+ ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_M0) + 1):ncol(G11)] = apply((h2[, (4 * j - 1)]/data$rmpw * (-data$M * data$pM0/data$pM1 * (1 - data$pM1) + (1 - data$M) * (1 - data$pM0)/(1 - data$pM1) * data$pM1) * V_M1)[data$R==1, ], 2, sum)/nrow(data)
    }
    G12 = matrix(0, nrow(G11), ncol(G22))
    G= cbind(rbind(G11, G21), rbind(G12, G22))
    V = 1/N * ginv(G) %*% H %*% t(ginv(G))
    V_mu = V[(ncol(V) - 4 * J + 1):ncol(V), (ncol(V) - 4 * J + 1):ncol(V)]  
    
    ######## Population Average Effect Estimation
    Phi = matrix(c(-1,-1,0,-1,1,0,0,0,1,-1,0,1,-1,0,-1,1,0,1,0,1),5, 4)
    beta = kronecker(diag(1, J), Phi) %*% mu
    V_beta = kronecker(diag(1, J), Phi) %*% V_mu %*% kronecker(diag(1, J), t(Phi))
    
    Psi = kronecker(rep(1, J), diag(1, 5))
    gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% beta
    
    if(permutation == F){
      ######## Between-Site Variance Estimation
      B = 0
      W = 1/(J * (J - 1)) * t(Psi) %*% V_beta %*% Psi
      for(j in 1:J){
        B = B + 1/(J - 1) * (beta[(5 * j - 4):(5 * j)] - gamma) %*% t((beta[(5 * j - 4):(5 * j)] - gamma))
        W = W - 1/(J - 1) * V_beta[(5 * j - 4):(5 * j), (5 * j - 4):(5 * j)]
      }
      tau = B + W
      
      ######## Standard Error of Population Average Effects
      V_gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% (V_beta + kronecker(diag(1, J), tau)) %*% Psi %*% solve(t(Psi) %*% Psi)
      SE_gamma = sqrt(diag(V_gamma))
      t_gamma = gamma/SE_gamma
      p_gamma = (1 - pnorm(abs(t_gamma))) * 2
      est_gamma = cbind(gamma, SE_gamma, t_gamma, p_gamma)
      est_gamma = as.data.frame(round(est_gamma, 3))
      
      sig = NULL
      sig[est_gamma[, 4] <= 0.001] = "**"
      sig[est_gamma[, 4] > 0.001 & est_gamma[, 4] <= 0.01] = "*"
      sig[est_gamma[, 4] > 0.01 & est_gamma[, 4] <= 0.05] = "."
      sig[est_gamma[, 4] > 0.05] = ""
      est_gamma = cbind(est_gamma, sig)
      est_gamma[est_gamma[, 4] < 0.001, 4] = "<0.001"
      est_gamma = as.data.frame(est_gamma)
      colnames(est_gamma) = c("Estimate", "Std.Error", "t value", "Pr(>|t|)", "")
      rownames(est_gamma) = c("ITT Effect on Outcome","Natural Direct Effect", "Natural Indirect Effect", "Pure Indirect Effect", "T-M Interaction Effect")
      
      ######## Between-Site Variance Output
      tau.ori = tau
      tau.sub = tau
      
      if("TRUE" %in% (diag(tau.ori) < 0)){
        tau[which(diag(tau < 0)), ] = 0
        tau[, which(diag(tau < 0))] = 0
        tau.sub = tau.ori[-which(diag(tau.ori < 0)), ]
        tau.sub = tau.sub[, -which(diag(tau.sub < 0))]    
      }
      
      cor.sub = cov2cor(tau.sub)
      if(max(na.omit(cor.sub)) > 1)
        cor.sub = suppressWarnings(cor.smooth(cor.sub))
      cor = cor.sub
      
      if("TRUE" %in% (diag(tau.ori) < 0)){
        cor = tau.ori
        cor[which(diag(tau.ori < 0)), ] = NaN
        cor[, which(diag(tau.ori < 0))] = NaN
        cor[which(diag(tau.ori >= 0)), which(diag(tau.ori >= 0))] = cor.sub
      }
      
      var = round(cbind(diag(tau), sqrt(diag(tau)), cor[, 1:4]), 3)
      var[1, 3:6] = ""
      var[2, 4:6] = ""
      var[3, 5:6] = ""
      var[4, 6] = ""
      var = as.data.frame(var)
      colnames(var) = c("Variance", "Std.Dev.", "Corr", "", "", "")
      rownames(var) = c("ITT Effect on Outcome","Natural Direct Effect", "Natural Indirect Effect", "Pure Indirect Effect", "T-M Interaction Effect")
    }   
    
    ######## Chisq for Between-Site Variance Testing
    chisq = 0
    for (j in 1:J) {
      Vj = V_beta[(5 * j - 4):(5 * j), (5 * j - 4):(5 * j)]
      chisq = chisq + (beta[(5 * j - 4):(5 * j)] - gamma)^2/diag(Vj)
    }
    
    if(permutation == F)
      return(list(Random_effects = var, Fixed_effects = est_gamma, chisq = chisq))
    if(permutation == T)
      return(list(chisq = chisq))
  }
  
  result = suppressWarnings({est(data, y, treatment, mediator, response, XR1, XR0, XM1, XM0, site, sample.weight, permutation = F)})
  
  siteID = data$site
  data_pm = data
  chisq_pm = NULL
  for (i in 1:npermute) {
    try({
      data_pm$site = sample(siteID, length(siteID), replace = F)
      chisq_pm = cbind(chisq_pm, suppressWarnings({est(data_pm, y, treatment, mediator, response, XR1, XR0, XM1, XM0, site, sample.weight, permutation = T)})$chisq)
    }, silent = T)
  }
  pvalue = round(apply((chisq_pm - matrix(rep(result$chisq, ncol(chisq_pm)), 5, ncol(chisq_pm))) >= 0, 1, mean), 3)
  
  sig = NULL
  sig[pvalue <= 0.001] = "**"
  sig[pvalue > 0.001 & pvalue <= 0.01] = "*"
  sig[pvalue > 0.01 & pvalue <= 0.05] = "."
  sig[pvalue > 0.05] = ""
  result$Random_effects = cbind(result$Random_effects, pvalue, sig)
  colnames(result$Random_effects)[4:8] = ""
  colnames(result$Random_effects)[7] = "p-value"
  
  return(list(Random_effects = result$Random_effects, Fixed_effects = result$Fixed_effects))
}

#' Complex multisite causal mediation analysis with two concurrent mediators in the presence of complex sample and survey designs and non-random nonresponse
#' 
#' This function is used to estimate both the population average and between-site variance of a natural indirect transmitted through each mediator and a natural direct effect when two concurrent (conditionally independent) mediators are involved in the mediation mechanism. 
#' It incorporates a sample weight to adjust for complex sample and survey designs and employs an estimated nonresponse weight to account for non-random nonresponse.
#'
#' @param data The data set for analysis.
#' @param y The name of the outcome variable (string).
#' @param treatment The name of the treatment variable (string).
#' @param mediator A vecor of two mediator variable names (string).
#' @param response The name of the response variable (string), which is equal to 1 if the individual responded and 0 otherwise.
#' @param XR1 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XR0 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM11 A vector of variable names (string) of pretreatment covariates in the propensity score model for the first mediator under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM10 A vector of variable names (string) of pretreatment covariates in the propensity score model for the first mediator under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM21 A vector of variable names (string) of pretreatment covariates in the propensity score model for the second mediator under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM20 A vector of variable names (string) of pretreatment covariates in the propensity score model for the second mediator under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param site The variable name for the site ID (string).
#' @param sample.weight The variable name for the sample weight given by design (string).
#' @return A list contains the estimates of the between-site variances of indirect effect via M1 given M2(0) ("I.M1(0)"), indirect effect via M2 given M1(1) ("I.M2(1)"), direct effect ("D"), indirect effect via M1 given M2(1) ("I.M1(1)"), indirect effect via M2 given M1(0) ("I.M2(0)"), interaction effect between M1 and M2 ("I.M1*M2"), and the correlations between the effects across sites ($Random_effects), and the population average effect estimates along with their hypothesis testing results ($Fixed_effects).
#' @author Xu Qin, Jonah Deutsch, and Guanglei Hong
#' @references Qin, X., Deutsch, J, & Hong, G. (2021). Unpacking Complex Mediation Mechanisms and Their Heterogeneity between Sites in A Job Corps Evaluation. The Journal of Policy Analysis and Management, 40(1), 158-190. \doi{10.1002/pam.22268}
#' @export
#' @importFrom stats as.formula binomial coef fitted lm model.matrix pnorm na.omit cov2cor
#' @importFrom lme4 VarCorr fixef glmer ranef glmerControl
#' @importFrom statmod gauss.quad.prob
#' @importFrom psych cor.smooth
#' @importFrom MASS ginv
#' @examples 
#' data(sim.weights)
#'
#' msmediate.concurrent(data = sim.weights, y = "y", treatment = "tr", mediator = c("me", "me2"), response = "R", XR1 = c("x1", "x2", "x3"), XR0 = c("x1", "x2", "x3"), XM11 = c("x1", "x2", "x3"), XM10 = c("x1", "x2", "x3"), XM21 = c("x1", "x2", "x3"), XM20 = c("x1", "x2", "x3"), site = "site", sample.weight = "WD")
#'
msmediate.concurrent = function(data, y, treatment, mediator, response, XR1, XR0, XM11, XM10, XM21, XM20, site, sample.weight){
  data = as.data.frame(data)
  colnames(data)[which(colnames(data) == y)] = "y"
  colnames(data)[which(colnames(data) == treatment)] = "tr"
  colnames(data)[which(colnames(data) == mediator[1])] = "MV"
  colnames(data)[which(colnames(data) == mediator[2])] = "ME"
  colnames(data)[which(colnames(data) == response)] = "R"
  XMV1 = XM11 
  XMV0 = XM10
  XME1 = XM21 
  XME0 = XM20
  colnames(data)[which(colnames(data) == site)] = "site"
  colnames(data)[which(colnames(data) == sample.weight)] = "WD"
  data = data[order(data$site), ]
  
  ######## Separate the data set into two, one for the Job Corps group and the other for the control group.
  data1 = data[data$tr == 1, ]
  data0 = data[data$tr == 0, ]
  ######## Nonresponse Weight Estimation in Step 1
  #### Fit multilevel logistic regressions of the response indicator
  ## Sample: full sample
  ## Model fitted to the Job Corps group:
  lR1 = glmer(as.formula(paste("R", "~", paste(XR1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  ## Model fitted to the control group:
  lR0 = glmer(as.formula(paste("R", "~", paste(XR0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  #### Predict the response probabilities
  data$pR[data$tr == 1 & data$R == 1] = fitted(lR1)[data1$R == 1]
  # data$pR[data$tr == 1 & data$R == 0] = 1 - fitted(lR1)[data1$R == 0] #This is for balance checking, while it is not used in the data analysis, which is focused on R = 1.
  data$pR[data$tr == 0 & data$R == 1] = fitted(lR0)[data0$R == 1]
  # data$pR[data$tr == 0 & data$R == 0] = 1 - fitted(lR0)[data0$R == 0] #This is for balance checking, while it is not used in the data analysis, which is focused on R = 1.
  #### Numerator
  lR.nu1 = glmer(R ~ 1 + (1|site), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  lR.nu0 = glmer(R ~ 1 + (1|site), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  data$pR.nu[data$tr == 1 & data$R == 1] = fitted(lR.nu1)[data1$R == 1]
  data$pR.nu[data$tr == 0 & data$R == 1] = fitted(lR.nu0)[data0$R == 1]
  
  #### Construct the nonresponse weight
  data$WR[data$R == 1] = data$pR.nu[data$R == 1]/data$pR[data$R == 1]
  data$WR[data$R == 0] = 1 # Otherwise G is not invertible due to NA's. Because WR will not be used for nonrespondents, this will not affect the final results.
  data$pR[data$R == 0] = 1
  data$pR.nu[data$R == 0] = 1
  
  #### Moment functions for nonresponse weight estimation in step 1
  ## Generate G-H quadrature points and weights
  nnodes = 10
  temp = gauss.quad.prob(nnodes, "normal") # Here, random effects are assumed normal
  nodes = temp$nodes
  weights = temp$weights
  # When there is one random effect, i.e. when the estimated variance of one random effect is 0
  B_1 = t(as.matrix(nodes))
  A_1 = weights
  ## Moment function
  data$R.nu = data$R # This is for the use in fjq, h_pij, h_sigmaj
  h1_fun = function(model, group, data){ # model = "M" for mediator model; model = "R" for response model; group = 0 for control group; group = 1 for treatment group.
    l = get(paste0("l", model, group))
    if(model == "R.nu"){
      x = as.matrix(cbind(rep(1, nrow(data))))
    } else {
      x = as.matrix(cbind(rep(1, nrow(data)), data[, get(paste0("X", model, group))]))
    }     
    pi = fixef(l)
    if(round(VarCorr(l)$site[1], 5) == 0){
      sigma = 0
      theta = matrix(0, length(unique(data$site)), 1)
      A = 1
      B = as.matrix(0)
    } else {
      sigma = sqrt(VarCorr(l)$site[1])
      theta = as.matrix(ranef(l)$site)/sigma
      A = A_1
      B = B_1
    }
    p = NULL
    for(q in 1:length(A)){
      p = cbind(p, 1 / (1 + exp(-(x %*% pi + sigma * B[, q]))))
    }
    p = as.matrix(p)
    
    h_pi = NULL
    h_sigma = NULL
    for(j in unique(data$site)){
      data$theta[data$site == j] = theta[as.numeric(rownames(theta)) == j, ]
      lj = 0
      h_pij = 0
      h_sigmaj = 0
      for(q in 1:length(A)){
        pq = p[, q]
        fjq = 1
        for(i in which(data$site == j)){
          fjq = fjq * (pq[i]^data[i, model] * (1 - pq[i])^(1 - data[i, model]))^(data$tr[i] == group)
        }
        lj = lj + fjq * A[q]
        h_pij = h_pij + fjq * A[q] * (data$tr[data$site == j] == group) * (data[data$site == j, model] - pq[data$site == j]) * as.matrix(x[data$site == j,])      
        h_sigmaj =  h_sigmaj + as.matrix(fjq * A[q] * (data$tr[data$site == j] == group) * (data[data$site == j, model] - pq[data$site == j])) * B[, q]
      }
      h_pi = rbind(h_pi, h_pij/lj)
      h_sigma = rbind(h_sigma, h_sigmaj/lj)
    }
    
    if(round(VarCorr(l)$site[1], 5) == 0){
      h1 = h_pi
      V = x
    } else {
      h1 = cbind(h_pi, h_sigma)
      V = cbind(x, data$theta)
    }
    
    return(list(h1 = h1, V = V, sigma = sigma))
  }
  
  h1_R1_list = h1_fun(model = "R", group = 1, data = data) #denominator for treatment group 
  h1_R0_list = h1_fun(model = "R", group = 0, data = data) #denominator for control group
  h1_R1 = h1_R1_list$h1
  h1_R0 = h1_R0_list$h1
  V_R1 = h1_R1_list$V
  V_R0 = h1_R0_list$V
  sigma_R1 = h1_R1_list$sigma
  sigma_R0 = h1_R0_list$sigma
  
  h1_R.nu1_list = h1_fun(model = "R.nu", group = 1, data = data) #numerator for treatment group 
  h1_R.nu0_list = h1_fun(model = "R.nu", group = 0, data = data) #numerator for control group
  h1_R.nu1 = h1_R.nu1_list$h1
  h1_R.nu0 = h1_R.nu0_list$h1
  V_R.nu1 = h1_R.nu1_list$V
  V_R.nu0 = h1_R.nu0_list$V
  sigma_R.nu1 = h1_R.nu1_list$sigma
  sigma_R.nu0 = h1_R.nu0_list$sigma
  
  ######## RMPW Estimation in Step 1
  #### Fit multilevel logistic regressions of the mediator
  ## Sample: respondents in the 48-month interview
  data1 = data[data$tr == 1 & data$R == 1, ]
  data0 = data[data$tr == 0 & data$R == 1, ]
  
  ## Model fitted to each mediator in each of two different treatment groups:
  lMV1 = glmer(as.formula(paste("MV ~", paste(XMV1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  lMV0 = glmer(as.formula(paste("MV ~", paste(XMV0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  lME1 = glmer(as.formula(paste("ME ~", paste(XME1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  lME0 = glmer(as.formula(paste("ME ~", paste(XME0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  
  XMV1 = colnames(model.matrix(lMV1))[-1]
  XMV0 = colnames(model.matrix(lMV0))[-1]
  XME1 = colnames(model.matrix(lME1))[-1]
  XME0 = colnames(model.matrix(lME0))[-1]
  
  rmpw = function(lM1, lM0, M, XM1, XM0){
    #### Predict the mediator probabilities
    predict.lmer = function(l, data, X, ranX) {
      nj = NULL
      for (j in as.numeric(rownames(ranef(l)$site))) {
        nj = c(nj, sum(data$site == j))
      }
      ranef = NULL
      for (k in 1:dim(ranef(l)$site)[2]) {
        ranef = cbind(ranef, rep(ranef(l)$site[, k], nj))
      }
      X = cbind(rep(1, nrow(data)), data[, X])
      pred_logit = as.matrix(X) %*% fixef(l) + apply(ranX * ranef, 1, sum)
      pred = exp(pred_logit)/(exp(pred_logit) + 1)
    }
    data$pM1[data$tr == 1 & data$R == 1] = fitted(lM1)
    data$pM0[data$tr == 0 & data$R == 1] = fitted(lM0)
    data$pM1[data$tr == 0 & data$R == 1] = predict.lmer(lM1, data0, XM1, 1) 
    data$pM0[data$tr == 1 & data$R == 1] = predict.lmer(lM0, data1, XM0, 1) 
    #### Construct the RMPW weight
    data$rmpw[data$tr == 1 & data[, M] == 1 & data$R == 1] = data$pM0[data$tr == 1 & data[, M] == 1 & data$R == 1]/data$pM1[data$tr == 1 & data[, M] == 1 & data$R == 1]
    data$rmpw[data$tr == 1 & data[, M] == 0 & data$R == 1] = (1 - data$pM0[data$tr == 1 & data[, M] == 0 & data$R == 1])/(1 - data$pM1[data$tr == 1 & data[, M] == 0 & data$R == 1])
    data$rmpw[data$tr == 0 & data[, M] == 1 & data$R == 1] = data$pM1[data$tr == 0 & data[, M] == 1 & data$R == 1]/data$pM0[data$tr == 0 & data[, M] == 1 & data$R == 1]
    data$rmpw[data$tr == 0 & data[, M] == 0 & data$R == 1] = (1 - data$pM1[data$tr == 0 & data[, M] == 0 & data$R == 1])/(1 - data$pM0[data$tr == 0 & data[, M] == 0 & data$R == 1])
    data$rmpw[data$R == 0] = 0
    
    colnames(data)[which(colnames(data) == "pM1")] = paste0("p", M, "1")
    colnames(data)[which(colnames(data) == "pM0")] = paste0("p", M, "0")
    colnames(data)[which(colnames(data) == "rmpw")] = paste0("rmpw.", M)
    
    return(data)
  }
  data = rmpw(lMV1, lMV0, "MV", XMV1, XMV0)
  data = rmpw(lME1, lME0, "ME", XME1, XME0)
  
  #### Moment functions for RMPW weight estimation in step 1
  h1_MV1_list = h1_fun(model = "MV", group = 1, data = data[data$R == 1, ]) #Mediator model for treatment group 
  h1_MV0_list = h1_fun(model = "MV", group = 0, data = data[data$R == 1, ]) #Mediator model for control group
  h1_MV1 = matrix(0, nrow(data), ncol(h1_MV1_list$h1))
  h1_MV0 = matrix(0, nrow(data), ncol(h1_MV0_list$h1))
  h1_MV1[data$R == 1, ] = h1_MV1_list$h1 #Mediator is modeled based on respondents. Therefore, the moment function of the mediator model for nonrespondents is 0.
  h1_MV0[data$R == 1, ] = h1_MV0_list$h1
  V_MV1 = matrix(0, nrow(data), ncol(h1_MV1_list$V))
  V_MV0 = matrix(0, nrow(data), ncol(h1_MV0_list$V))
  V_MV1[data$R == 1, ] = h1_MV1_list$V
  V_MV0[data$R == 1, ] = h1_MV0_list$V
  sigma_MV1 = h1_MV1_list$sigma
  sigma_MV0 = h1_MV0_list$sigma
  
  h1_ME1_list = h1_fun(model = "ME", group = 1, data = data[data$R == 1, ]) #Mediator model for treatment group 
  h1_ME0_list = h1_fun(model = "ME", group = 0, data = data[data$R == 1, ]) #Mediator model for control group
  h1_ME1 = matrix(0, nrow(data), ncol(h1_ME1_list$h1))
  h1_ME0 = matrix(0, nrow(data), ncol(h1_ME0_list$h1))
  h1_ME1[data$R == 1, ] = h1_ME1_list$h1 #Mediator is modeled based on respondents. Therefore, the moment function of the mediator model for nonrespondents is 0.
  h1_ME0[data$R == 1, ] = h1_ME0_list$h1
  V_ME1 = matrix(0, nrow(data), ncol(h1_ME1_list$V))
  V_ME0 = matrix(0, nrow(data), ncol(h1_ME0_list$V))
  V_ME1[data$R == 1, ] = h1_ME1_list$V
  V_ME0[data$R == 1, ] = h1_ME0_list$V
  sigma_ME1 = h1_ME1_list$sigma
  sigma_ME0 = h1_ME0_list$sigma
  
  h1 = cbind(h1_R0, h1_R1, h1_R.nu0, h1_R.nu1, h1_MV0, h1_MV1, h1_ME0, h1_ME1)
  
  ######## Site-Specific Mean Potential Outcome Estimation in Step 2
  N = nrow(data)
  J = length(unique(data$site))
  h2 = matrix(0, N, J * 5)
  mu = NULL
  data$y[data$R == 0] = 0
  data$WD[data$R == 0] = 0
  for (j in unique(data$site)) {
    mu0.0.0 = sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$y)/sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR)
    mu1.1.1 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR)
    mu1.0.1 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV)
    mu1.1.0 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME)
    mu1.0.0 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME)
    
    h0.0.0 = data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * (data$y - mu0.0.0)
    h1.1.1 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * (data$y - mu1.1.1)
    h1.0.1 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * (data$y - mu1.0.1)
    h1.1.0 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME * (data$y - mu1.1.0)
    h1.0.0 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME * (data$y - mu1.0.0)
    
    mu = c(mu, mu0.0.0, mu1.1.1, mu1.0.1, mu1.1.0, mu1.0.0)
    h2[, (5 * (which(unique(data$site) == j) - 1) + 1):(5 * which(unique(data$site) == j))] = cbind(h0.0.0, h1.1.1, h1.0.1, h1.1.0, h1.0.0)
  }
  
  if(sum(is.na(mu)) > 0){
    site.omit = which(is.na(mu))[5]/5
    J = J - sum(is.na(mu))/5
    h2 = h2[, -which(is.na(mu))]
    mu = na.omit(mu)
  }
  
  ######## Asymptotic Sampling Variance of the Two-Step Estimators
  #### Stack the moment functions from both steps
  h = cbind(h1, h2)
  H = 1/N * t(as.matrix(h)) %*% as.matrix(h)
  bdiag = function(A, B){
    C = matrix(0, nrow = nrow(A) + nrow(B), ncol = ncol(A) + ncol(B))
    C[1:nrow(A),1:ncol(A)] = A
    C[(1 + nrow(A)):nrow(C),(1 + ncol(A)):ncol(C)] = B
    return(C)
  }
  G11 = as.matrix(bdiag(bdiag(bdiag(-1/N * t(h1_R0) %*% as.matrix(h1_R0), -1/N * t(h1_R1) %*% as.matrix(h1_R1)), bdiag(-1/N * t(h1_R.nu0) %*% as.matrix(h1_R.nu0), -1/N * t(h1_R.nu1) %*% as.matrix(h1_R.nu1))), bdiag(bdiag(-1/N * t(h1_MV0) %*% as.matrix(h1_MV0), -1/N * t(h1_MV1) %*% as.matrix(h1_MV1)), bdiag(-1/N * t(h1_ME0) %*% as.matrix(h1_ME0), -1/N * t(h1_ME1) %*% as.matrix(h1_ME1)))))
  G22_diag = NULL
  if(length(unique(data$site)) != J){
    for (j in unique(data$site)[-site.omit]) {
      G22_diag = c(G22_diag, mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME)
      )
    }
  } else {
    for (j in unique(data$site)) {
      G22_diag = c(G22_diag, mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME),
                   mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME)
      )
    }
  }
  G22 = diag(G22_diag)
  G21 = matrix(0, 5 * J, ncol(G11))
  for(j in 1:J){
    G21[5 * j - 4, 1:ncol(h1_R0)] = apply(h2[, (5 * j - 4)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R0, 2, mean)
    G21[5 * j - 4, (ncol(h1_R0) + ncol(h1_R1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0))] = apply(h2[, (5 * j - 4)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu0, 2, mean)
    G21[(5 * j - 3), (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (5 * j - 3)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
    G21[(5 * j - 3), (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (5 * j - 3)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
    
    G21[(5 * j - 2), (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (5 * j - 2)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
    G21[(5 * j - 1), (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (5 * j - 1)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
    G21[(5 * j), (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (5 * j)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
    
    G21[(5 * j - 2), (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (5 * j - 2)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
    G21[(5 * j - 1), (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (5 * j - 1)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
    G21[(5 * j), (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (5 * j)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
    
    G21[5 * j - 2, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0))] = apply((h2[, (5 * j - 2)]/data$rmpw.MV * (data$MV/data$pMV1 - (1 - data$MV)/(1 - data$pMV1)) * data$pMV0 * (1 - data$pMV0) * V_MV0)[data$R==1, ], 2, sum)/nrow(data)
    G21[5 * j - 2, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0) + ncol(h1_MV1))] = apply((h2[, (5 * j - 2)]/data$rmpw.MV * (-data$MV * data$pMV0/data$pMV1 * (1 - data$pMV1) + (1 - data$MV) * (1 - data$pMV0)/(1 - data$pMV1) * data$pMV1) * V_MV1)[data$R==1, ], 2, sum)/nrow(data)
    G21[5 * j - 1, (ncol(G11) - ncol(h1_ME0) - ncol(h1_ME1) + 1):(ncol(G11) - ncol(h1_ME1))] = apply((h2[, (5 * j - 1)]/data$rmpw.ME * (data$ME/data$pME1 - (1 - data$ME)/(1 - data$pME1)) * data$pME0 * (1 - data$pME0) * V_ME0)[data$R==1, ], 2, sum)/nrow(data)
    G21[5 * j - 1, (ncol(G11) - ncol(h1_ME1) + 1):ncol(G11)] = apply((h2[, (5 * j - 1)]/data$rmpw.ME * (-data$ME * data$pME0/data$pME1 * (1 - data$pME1) + (1 - data$ME) * (1 - data$pME0)/(1 - data$pME1) * data$pME1) * V_ME1)[data$R==1, ], 2, sum)/nrow(data)
    G21[5 * j, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0))] = apply((h2[, (5 * j)]/data$rmpw.MV * (data$MV/data$pMV1 - (1 - data$MV)/(1 - data$pMV1)) * data$pMV0 * (1 - data$pMV0) * V_MV0)[data$R==1, ], 2, sum)/nrow(data)
    G21[5 * j, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0) + ncol(h1_MV1))] = apply((h2[, (5 * j)]/data$rmpw.MV * (-data$MV * data$pMV0/data$pMV1 * (1 - data$pMV1) + (1 - data$MV) * (1 - data$pMV0)/(1 - data$pMV1) * data$pMV1) * V_MV1)[data$R==1, ], 2, sum)/nrow(data)
    G21[5 * j, (ncol(G11) - ncol(h1_ME0) - ncol(h1_ME1) + 1):(ncol(G11) - ncol(h1_ME1))] = apply((h2[, (5 * j)]/data$rmpw.ME * (data$ME/data$pME1 - (1 - data$ME)/(1 - data$pME1)) * data$pME0 * (1 - data$pME0) * V_ME0)[data$R==1, ], 2, sum)/nrow(data)
    G21[5 * j, (ncol(G11) - ncol(h1_ME1) + 1):ncol(G11)] = apply((h2[, (5 * j)]/data$rmpw.ME * (-data$ME * data$pME0/data$pME1 * (1 - data$pME1) + (1 - data$ME) * (1 - data$pME0)/(1 - data$pME1) * data$pME1) * V_ME1)[data$R==1, ], 2, sum)/nrow(data)
  }
  G12 = matrix(0, nrow(G11), ncol(G22))
  G= cbind(rbind(G11, G21), rbind(G12, G22))
  V = 1/N * ginv(G) %*% H %*% t(ginv(G))
  V_mu = V[(ncol(V) - 5 * J + 1):ncol(V), (ncol(V) - 5 * J + 1):ncol(V)]  
  
  ######## Population Average Effect Estimation
  Phi = matrix(c(0, 0, -1, 0, 0, 0, -1,
                 0, 1, 0, 1, 0, 1, 1,
                 0, 0, 0, -1, 1, -1, 0, 
                 1, -1, 0, 0, 0, -1, 0,
                 -1, 0, 1, 0, -1, 1, 0), 7, 5)
  beta = kronecker(diag(1, J), Phi) %*% mu
  V_beta = kronecker(diag(1, J), Phi) %*% V_mu %*% kronecker(diag(1, J), t(Phi))
  
  Psi = kronecker(rep(1, J), diag(1, 7))
  gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% beta
  
  ######## Between-Site Variance Estimation
  B = 0
  W = 1/(J * (J - 1)) * t(Psi) %*% V_beta %*% Psi
  for(j in 1:J){
    B = B + 1/(J - 1) * (beta[(7 * j - 6):(7 * j)] - gamma) %*% t((beta[(7 * j - 6):(7 * j)] - gamma))
    W = W - 1/(J - 1) * V_beta[(7 * j - 6):(7 * j), (7 * j - 6):(7 * j)]
  }
  tau = B + W
  
  ######## Standard Error of Population Average Effects
  V_gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% (V_beta + kronecker(diag(1, J), tau)) %*% Psi %*% solve(t(Psi) %*% Psi)
  SE_gamma = sqrt(diag(V_gamma))
  t_gamma = gamma/SE_gamma
  p_gamma = (1 - pnorm(abs(t_gamma))) * 2
  est_gamma = cbind(gamma, SE_gamma, t_gamma, p_gamma)
  est_gamma = as.data.frame(round(est_gamma, 3))
  colnames(est_gamma) = c("Estimate", "Std.Error", "t value", "Pr(>|t|)")
  rownames(est_gamma) = c("I.M1(0)", "I.M2(1)", "D", "I.M1(1)", "I.M2(0)", "I.M1*M2", "ITT")
  
  ######## Between-Site Variance Output
  tau.ori = tau
  tau.sub = tau
  
  if("TRUE" %in% (diag(tau.ori) < 0)){
    tau[which(diag(tau < 0)), ] = 0
    tau[, which(diag(tau < 0))] = 0
    tau.sub = tau.ori[-which(diag(tau.ori < 0)), ]
    tau.sub = tau.sub[, -which(diag(tau.sub < 0))]    
  }
  
  cor.sub = cov2cor(tau.sub)
  if(max(na.omit(cor.sub)) > 1)
    cor.sub = suppressWarnings(cor.smooth(cor.sub))
  cor = cor.sub
  
  if("TRUE" %in% (diag(tau.ori) < 0)){
    cor = tau.ori
    cor[which(diag(tau.ori < 0)), ] = NaN
    cor[, which(diag(tau.ori < 0))] = NaN
    cor[which(diag(tau.ori >= 0)), which(diag(tau.ori >= 0))] = cor.sub
  }
  
  var = round(cbind(diag(tau), sqrt(diag(tau)), cor[, 1:6]), 3)
  var[1, 3:8] = ""
  var[2, 4:8] = ""
  var[3, 5:8] = ""
  var[4, 6:8] = ""
  var[5, 7:8] = ""
  var[6, 8] = ""
  var = as.data.frame(var)
  colnames(var) = c("Variance", "Std.Dev.", "Corr", "", "", "", "", "")
  rownames(var) = c("I.M1(0)", "I.M2(1)", "D", "I.M1(1)", "I.M2(0)", "I.M1*M2", "ITT")
  
  return(list(Random_effects = var, Fixed_effects = est_gamma))
}

#' Variance testing for complex multisite causal mediation analysis with two concurrent mediators in the presence of complex sample and survey designs and non-random nonresponse
#'
#' This function performs hypothesis testing for the between-site variances besides providing the same output as given by the function msmediate.concurrent().
#' 
#' @param data The data set for analysis.
#' @param y The name of the outcome variable (string).
#' @param treatment The name of the treatment variable (string).
#' @param mediator A vecor of two mediator variable names (string).
#' @param response The name of the response variable (string), which is equal to 1 if the individual responded and 0 otherwise.
#' @param XR1 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XR0 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM11 A vector of variable names (string) of pretreatment covariates in the propensity score model for the first mediator under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM10 A vector of variable names (string) of pretreatment covariates in the propensity score model for the first mediator under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM21 A vector of variable names (string) of pretreatment covariates in the propensity score model for the second mediator under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM20 A vector of variable names (string) of pretreatment covariates in the propensity score model for the second mediator under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param site The variable name for the site ID (string).
#' @param sample.weight The variable name for the sample weight given by design (string).
#' @param npermute The number of permutations for the permutation test. The default value is 200. It may take a long time, depending on the sample size and the length of X.
#' @return A list contains the hypothesis testing results of the between-site variance of the causal effects, besides the same output as given by the function msmediate().
#' @author Xu Qin, Jonah Deutsch, and Guanglei Hong
#' @references Qin, X., Deutsch, J, & Hong, G. (2021). Unpacking Complex Mediation Mechanisms and Their Heterogeneity between Sites in A Job Corps Evaluation. The Journal of Policy Analysis and Management, 40(1), 158-190. \doi{10.1002/pam.22268}
#' @export
#' @importFrom stats as.formula binomial coef fitted lm model.matrix pnorm na.omit cov2cor
#' @importFrom lme4 VarCorr fixef glmer ranef glmerControl 
#' @importFrom statmod gauss.quad.prob
#' @importFrom psych cor.smooth
#' @importFrom MASS ginv
#' @examples 
#' data(sim.weights)
#' set.seed(1)
#' vartest.msmediate.concurrent(data = sim.weights, y = "y", treatment = "tr", mediator = c("me", "me2"), response = "R", XR1 = c("x1", "x2", "x3"), XR0 = c("x1", "x2", "x3"), XM11 = c("x1", "x2", "x3"), XM10 = c("x1", "x2", "x3"), XM21 = c("x1", "x2", "x3"), XM20 = c("x1", "x2", "x3"), site = "site", sample.weight = "WD", npermute = 2)
#'
vartest.msmediate.concurrent = function(data, y, treatment, mediator, response, XR1, XR0, XM11, XM10, XM21, XM20, site, sample.weight, npermute = 200) {
  data = as.data.frame(data)
  colnames(data)[which(colnames(data) == y)] = "y"
  colnames(data)[which(colnames(data) == treatment)] = "tr"
  colnames(data)[which(colnames(data) == mediator[1])] = "MV"
  colnames(data)[which(colnames(data) == mediator[2])] = "ME"
  colnames(data)[which(colnames(data) == response)] = "R"
  XMV1 = XM11 
  XMV0 = XM10
  XME1 = XM21 
  XME0 = XM20
  colnames(data)[which(colnames(data) == site)] = "site"
  colnames(data)[which(colnames(data) == sample.weight)] = "WD"
  
  est = function(data, y, treatment, mediator, response, XR1, XR0, XM11, XM10, XM21, XM20, site, sample.weight, permutation = F){
    data = data[order(data$site), ]
    
    ######## Separate the data set into two, one for the Job Corps group and the other for the control group.
    data1 = data[data$tr == 1, ]
    data0 = data[data$tr == 0, ]
    ######## Nonresponse Weight Estimation in Step 1
    #### Fit multilevel logistic regressions of the response indicator
    ## Sample: full sample
    ## Model fitted to the Job Corps group:
    lR1 = glmer(as.formula(paste("R", "~", paste(XR1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    ## Model fitted to the control group:
    lR0 = glmer(as.formula(paste("R", "~", paste(XR0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    #### Predict the response probabilities
    data$pR[data$tr == 1 & data$R == 1] = fitted(lR1)[data1$R == 1]
    # data$pR[data$tr == 1 & data$R == 0] = 1 - fitted(lR1)[data1$R == 0] #This is for balance checking, while it is not used in the data analysis, which is focused on R = 1.
    data$pR[data$tr == 0 & data$R == 1] = fitted(lR0)[data0$R == 1]
    # data$pR[data$tr == 0 & data$R == 0] = 1 - fitted(lR0)[data0$R == 0] #This is for balance checking, while it is not used in the data analysis, which is focused on R = 1.
    #### Numerator
    lR.nu1 = glmer(R ~ 1 + (1|site), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    lR.nu0 = glmer(R ~ 1 + (1|site), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    data$pR.nu[data$tr == 1 & data$R == 1] = fitted(lR.nu1)[data1$R == 1]
    data$pR.nu[data$tr == 0 & data$R == 1] = fitted(lR.nu0)[data0$R == 1]
    
    #### Construct the nonresponse weight
    data$WR[data$R == 1] = data$pR.nu[data$R == 1]/data$pR[data$R == 1]
    data$WR[data$R == 0] = 1 # Otherwise G is not invertible due to NA's. Because WR will not be used for nonrespondents, this will not affect the final results.
    data$pR[data$R == 0] = 1
    data$pR.nu[data$R == 0] = 1
    
    #### Moment functions for nonresponse weight estimation in step 1
    ## Generate G-H quadrature points and weights
    nnodes = 10
    temp = gauss.quad.prob(nnodes, "normal") # Here, random effects are assumed normal
    nodes = temp$nodes
    weights = temp$weights
    # When there is one random effect, i.e. when the estimated variance of one random effect is 0
    B_1 = t(as.matrix(nodes))
    A_1 = weights
    ## Moment function
    data$R.nu = data$R # This is for the use in fjq, h_pij, h_sigmaj
    h1_fun = function(model, group, data){ # model = "M" for mediator model; model = "R" for response model; group = 0 for control group; group = 1 for treatment group.
      l = get(paste0("l", model, group))
      if(model == "R.nu"){
        x = as.matrix(cbind(rep(1, nrow(data))))
      } else {
        x = as.matrix(cbind(rep(1, nrow(data)), data[, get(paste0("X", model, group))]))
      }     
      pi = fixef(l)
      if(round(VarCorr(l)$site[1], 5) == 0){
        sigma = 0
        theta = matrix(0, length(unique(data$site)), 1)
        A = 1
        B = as.matrix(0)
      } else {
        sigma = sqrt(VarCorr(l)$site[1])
        theta = as.matrix(ranef(l)$site)/sigma
        A = A_1
        B = B_1
      }
      p = NULL
      for(q in 1:length(A)){
        p = cbind(p, 1 / (1 + exp(-(x %*% pi + sigma * B[, q]))))
      }
      p = as.matrix(p)
      
      h_pi = NULL
      h_sigma = NULL
      for(j in unique(data$site)){
        data$theta[data$site == j] = theta[as.numeric(rownames(theta)) == j, ]
        lj = 0
        h_pij = 0
        h_sigmaj = 0
        for(q in 1:length(A)){
          pq = p[, q]
          fjq = 1
          for(i in which(data$site == j)){
            fjq = fjq * (pq[i]^data[i, model] * (1 - pq[i])^(1 - data[i, model]))^(data$tr[i] == group)
          }
          lj = lj + fjq * A[q]
          h_pij = h_pij + fjq * A[q] * (data$tr[data$site == j] == group) * (data[data$site == j, model] - pq[data$site == j]) * as.matrix(x[data$site == j,])      
          h_sigmaj =  h_sigmaj + as.matrix(fjq * A[q] * (data$tr[data$site == j] == group) * (data[data$site == j, model] - pq[data$site == j])) * B[, q]
        }
        h_pi = rbind(h_pi, h_pij/lj)
        h_sigma = rbind(h_sigma, h_sigmaj/lj)
      }
      
      if(round(VarCorr(l)$site[1], 5) == 0){
        h1 = h_pi
        V = x
      } else {
        h1 = cbind(h_pi, h_sigma)
        V = cbind(x, data$theta)
      }
      
      return(list(h1 = h1, V = V, sigma = sigma))
    }
    
    h1_R1_list = h1_fun(model = "R", group = 1, data = data) #denominator for treatment group 
    h1_R0_list = h1_fun(model = "R", group = 0, data = data) #denominator for control group
    h1_R1 = h1_R1_list$h1
    h1_R0 = h1_R0_list$h1
    V_R1 = h1_R1_list$V
    V_R0 = h1_R0_list$V
    sigma_R1 = h1_R1_list$sigma
    sigma_R0 = h1_R0_list$sigma
    
    h1_R.nu1_list = h1_fun(model = "R.nu", group = 1, data = data) #numerator for treatment group 
    h1_R.nu0_list = h1_fun(model = "R.nu", group = 0, data = data) #numerator for control group
    h1_R.nu1 = h1_R.nu1_list$h1
    h1_R.nu0 = h1_R.nu0_list$h1
    V_R.nu1 = h1_R.nu1_list$V
    V_R.nu0 = h1_R.nu0_list$V
    sigma_R.nu1 = h1_R.nu1_list$sigma
    sigma_R.nu0 = h1_R.nu0_list$sigma
    
    ######## RMPW Estimation in Step 1
    #### Fit multilevel logistic regressions of the mediator
    ## Sample: respondents in the 48-month interview
    data1 = data[data$tr == 1 & data$R == 1, ]
    data0 = data[data$tr == 0 & data$R == 1, ]
    
    ## Model fitted to each mediator in each of two different treatment groups:
    lMV1 = glmer(as.formula(paste("MV ~", paste(XMV1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    lMV0 = glmer(as.formula(paste("MV ~", paste(XMV0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    lME1 = glmer(as.formula(paste("ME ~", paste(XME1, collapse="+"), "+(1|site)")), data = data1, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    lME0 = glmer(as.formula(paste("ME ~", paste(XME0, collapse="+"), "+(1|site)")), data = data0, family = binomial, nAGQ = 1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    
    XMV1 = colnames(model.matrix(lMV1))[-1]
    XMV0 = colnames(model.matrix(lMV0))[-1]
    XME1 = colnames(model.matrix(lME1))[-1]
    XME0 = colnames(model.matrix(lME0))[-1]
    
    rmpw = function(lM1, lM0, M, XM1, XM0){
      #### Predict the mediator probabilities
      predict.lmer = function(l, data, X, ranX) {
        nj = NULL
        for (j in as.numeric(rownames(ranef(l)$site))) {
          nj = c(nj, sum(data$site == j))
        }
        ranef = NULL
        for (k in 1:dim(ranef(l)$site)[2]) {
          ranef = cbind(ranef, rep(ranef(l)$site[, k], nj))
        }
        X = cbind(rep(1, nrow(data)), data[, X])
        pred_logit = as.matrix(X) %*% fixef(l) + apply(ranX * ranef, 1, sum)
        pred = exp(pred_logit)/(exp(pred_logit) + 1)
      }
      data$pM1[data$tr == 1 & data$R == 1] = fitted(lM1)
      data$pM0[data$tr == 0 & data$R == 1] = fitted(lM0)
      data$pM1[data$tr == 0 & data$R == 1] = predict.lmer(lM1, data0, XM1, 1) 
      data$pM0[data$tr == 1 & data$R == 1] = predict.lmer(lM0, data1, XM0, 1) 
      #### Construct the RMPW weight
      data$rmpw[data$tr == 1 & data[, M] == 1 & data$R == 1] = data$pM0[data$tr == 1 & data[, M] == 1 & data$R == 1]/data$pM1[data$tr == 1 & data[, M] == 1 & data$R == 1]
      data$rmpw[data$tr == 1 & data[, M] == 0 & data$R == 1] = (1 - data$pM0[data$tr == 1 & data[, M] == 0 & data$R == 1])/(1 - data$pM1[data$tr == 1 & data[, M] == 0 & data$R == 1])
      data$rmpw[data$tr == 0 & data[, M] == 1 & data$R == 1] = data$pM1[data$tr == 0 & data[, M] == 1 & data$R == 1]/data$pM0[data$tr == 0 & data[, M] == 1 & data$R == 1]
      data$rmpw[data$tr == 0 & data[, M] == 0 & data$R == 1] = (1 - data$pM1[data$tr == 0 & data[, M] == 0 & data$R == 1])/(1 - data$pM0[data$tr == 0 & data[, M] == 0 & data$R == 1])
      data$rmpw[data$R == 0] = 0
      
      colnames(data)[which(colnames(data) == "pM1")] = paste0("p", M, "1")
      colnames(data)[which(colnames(data) == "pM0")] = paste0("p", M, "0")
      colnames(data)[which(colnames(data) == "rmpw")] = paste0("rmpw.", M)
      
      return(data)
    }
    data = rmpw(lMV1, lMV0, "MV", XMV1, XMV0)
    data = rmpw(lME1, lME0, "ME", XME1, XME0)
    
    #### Moment functions for RMPW weight estimation in step 1
    h1_MV1_list = h1_fun(model = "MV", group = 1, data = data[data$R == 1, ]) #Mediator model for treatment group 
    h1_MV0_list = h1_fun(model = "MV", group = 0, data = data[data$R == 1, ]) #Mediator model for control group
    h1_MV1 = matrix(0, nrow(data), ncol(h1_MV1_list$h1))
    h1_MV0 = matrix(0, nrow(data), ncol(h1_MV0_list$h1))
    h1_MV1[data$R == 1, ] = h1_MV1_list$h1 #Mediator is modeled based on respondents. Therefore, the moment function of the mediator model for nonrespondents is 0.
    h1_MV0[data$R == 1, ] = h1_MV0_list$h1
    V_MV1 = matrix(0, nrow(data), ncol(h1_MV1_list$V))
    V_MV0 = matrix(0, nrow(data), ncol(h1_MV0_list$V))
    V_MV1[data$R == 1, ] = h1_MV1_list$V
    V_MV0[data$R == 1, ] = h1_MV0_list$V
    sigma_MV1 = h1_MV1_list$sigma
    sigma_MV0 = h1_MV0_list$sigma
    
    h1_ME1_list = h1_fun(model = "ME", group = 1, data = data[data$R == 1, ]) #Mediator model for treatment group 
    h1_ME0_list = h1_fun(model = "ME", group = 0, data = data[data$R == 1, ]) #Mediator model for control group
    h1_ME1 = matrix(0, nrow(data), ncol(h1_ME1_list$h1))
    h1_ME0 = matrix(0, nrow(data), ncol(h1_ME0_list$h1))
    h1_ME1[data$R == 1, ] = h1_ME1_list$h1 #Mediator is modeled based on respondents. Therefore, the moment function of the mediator model for nonrespondents is 0.
    h1_ME0[data$R == 1, ] = h1_ME0_list$h1
    V_ME1 = matrix(0, nrow(data), ncol(h1_ME1_list$V))
    V_ME0 = matrix(0, nrow(data), ncol(h1_ME0_list$V))
    V_ME1[data$R == 1, ] = h1_ME1_list$V
    V_ME0[data$R == 1, ] = h1_ME0_list$V
    sigma_ME1 = h1_ME1_list$sigma
    sigma_ME0 = h1_ME0_list$sigma
    
    h1 = cbind(h1_R0, h1_R1, h1_R.nu0, h1_R.nu1, h1_MV0, h1_MV1, h1_ME0, h1_ME1)
    
    ######## Site-Specific Mean Potential Outcome Estimation in Step 2
    N = nrow(data)
    J = length(unique(data$site))
    h2 = matrix(0, N, J * 5)
    mu = NULL
    data$y[data$R == 0] = 0
    data$WD[data$R == 0] = 0
    for (j in unique(data$site)) {
      mu0.0.0 = sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * data$y)/sum(data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR)
      mu1.1.1 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR)
      mu1.0.1 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV)
      mu1.1.0 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME)
      mu1.0.0 = sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME * data$y)/sum(data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME)
      
      h0.0.0 = data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR * (data$y - mu0.0.0)
      h1.1.1 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * (data$y - mu1.1.1)
      h1.0.1 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * (data$y - mu1.0.1)
      h1.1.0 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME * (data$y - mu1.1.0)
      h1.0.0 = data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME * (data$y - mu1.0.0)
      
      mu = c(mu, mu0.0.0, mu1.1.1, mu1.0.1, mu1.1.0, mu1.0.0)
      h2[, (5 * (which(unique(data$site) == j) - 1) + 1):(5 * which(unique(data$site) == j))] = cbind(h0.0.0, h1.1.1, h1.0.1, h1.1.0, h1.0.0)
    }
    
    if(sum(is.na(mu)) > 0){
      site.omit = which(is.na(mu))[5]/5
      J = J - sum(is.na(mu))/5
      h2 = h2[, -which(is.na(mu))]
      mu = na.omit(mu)
    }
    
    ######## Asymptotic Sampling Variance of the Two-Step Estimators
    #### Stack the moment functions from both steps
    h = cbind(h1, h2)
    H = 1/N * t(as.matrix(h)) %*% as.matrix(h)
    bdiag = function(A, B){
      C = matrix(0, nrow = nrow(A) + nrow(B), ncol = ncol(A) + ncol(B))
      C[1:nrow(A),1:ncol(A)] = A
      C[(1 + nrow(A)):nrow(C),(1 + ncol(A)):ncol(C)] = B
      return(C)
    }
    G11 = as.matrix(bdiag(bdiag(bdiag(-1/N * t(h1_R0) %*% as.matrix(h1_R0), -1/N * t(h1_R1) %*% as.matrix(h1_R1)), bdiag(-1/N * t(h1_R.nu0) %*% as.matrix(h1_R.nu0), -1/N * t(h1_R.nu1) %*% as.matrix(h1_R.nu1))), bdiag(bdiag(-1/N * t(h1_MV0) %*% as.matrix(h1_MV0), -1/N * t(h1_MV1) %*% as.matrix(h1_MV1)), bdiag(-1/N * t(h1_ME0) %*% as.matrix(h1_ME0), -1/N * t(h1_ME1) %*% as.matrix(h1_ME1)))))
    G22_diag = NULL
    if(length(unique(data$site)) != J){
      for (j in unique(data$site)[-site.omit]) {
        G22_diag = c(G22_diag, mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME)
        )
      }
    } else {
      for (j in unique(data$site)) {
        G22_diag = c(G22_diag, mean(-data$R * (data$site == j) * (data$tr == 0) * data$WD * data$WR),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.ME),
                     mean(-data$R * (data$site == j) * (data$tr == 1) * data$WD * data$WR * data$rmpw.MV * data$rmpw.ME)
        )
      }
    }
    G22 = diag(G22_diag)
    G21 = matrix(0, 5 * J, ncol(G11))
    for(j in 1:J){
      G21[5 * j - 4, 1:ncol(h1_R0)] = apply(h2[, (5 * j - 4)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R0, 2, mean)
      G21[5 * j - 4, (ncol(h1_R0) + ncol(h1_R1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0))] = apply(h2[, (5 * j - 4)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu0, 2, mean)
      G21[(5 * j - 3), (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (5 * j - 3)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
      G21[(5 * j - 3), (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (5 * j - 3)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
      
      G21[(5 * j - 2), (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (5 * j - 2)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
      G21[(5 * j - 1), (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (5 * j - 1)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
      G21[(5 * j), (ncol(h1_R0) + 1):(ncol(h1_R0) + ncol(h1_R1))] = apply(h2[, (5 * j)]/data$WR * (-data$pR.nu * (1 - data$pR)/data$pR) * V_R1, 2, mean)
      
      G21[(5 * j - 2), (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (5 * j - 2)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
      G21[(5 * j - 1), (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (5 * j - 1)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
      G21[(5 * j), (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1))] = apply(h2[, (5 * j)]/data$WR * (data$pR.nu * (1 - data$pR.nu)/data$pR) * V_R.nu1, 2, mean)
      
      G21[5 * j - 2, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0))] = apply((h2[, (5 * j - 2)]/data$rmpw.MV * (data$MV/data$pMV1 - (1 - data$MV)/(1 - data$pMV1)) * data$pMV0 * (1 - data$pMV0) * V_MV0)[data$R==1, ], 2, sum)/nrow(data)
      G21[5 * j - 2, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0) + ncol(h1_MV1))] = apply((h2[, (5 * j - 2)]/data$rmpw.MV * (-data$MV * data$pMV0/data$pMV1 * (1 - data$pMV1) + (1 - data$MV) * (1 - data$pMV0)/(1 - data$pMV1) * data$pMV1) * V_MV1)[data$R==1, ], 2, sum)/nrow(data)
      G21[5 * j - 1, (ncol(G11) - ncol(h1_ME0) - ncol(h1_ME1) + 1):(ncol(G11) - ncol(h1_ME1))] = apply((h2[, (5 * j - 1)]/data$rmpw.ME * (data$ME/data$pME1 - (1 - data$ME)/(1 - data$pME1)) * data$pME0 * (1 - data$pME0) * V_ME0)[data$R==1, ], 2, sum)/nrow(data)
      G21[5 * j - 1, (ncol(G11) - ncol(h1_ME1) + 1):ncol(G11)] = apply((h2[, (5 * j - 1)]/data$rmpw.ME * (-data$ME * data$pME0/data$pME1 * (1 - data$pME1) + (1 - data$ME) * (1 - data$pME0)/(1 - data$pME1) * data$pME1) * V_ME1)[data$R==1, ], 2, sum)/nrow(data)
      G21[5 * j, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0))] = apply((h2[, (5 * j)]/data$rmpw.MV * (data$MV/data$pMV1 - (1 - data$MV)/(1 - data$pMV1)) * data$pMV0 * (1 - data$pMV0) * V_MV0)[data$R==1, ], 2, sum)/nrow(data)
      G21[5 * j, (ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0) + 1):(ncol(h1_R0) + ncol(h1_R1) + ncol(h1_R.nu0) + ncol(h1_R.nu1) + ncol(h1_MV0) + ncol(h1_MV1))] = apply((h2[, (5 * j)]/data$rmpw.MV * (-data$MV * data$pMV0/data$pMV1 * (1 - data$pMV1) + (1 - data$MV) * (1 - data$pMV0)/(1 - data$pMV1) * data$pMV1) * V_MV1)[data$R==1, ], 2, sum)/nrow(data)
      G21[5 * j, (ncol(G11) - ncol(h1_ME0) - ncol(h1_ME1) + 1):(ncol(G11) - ncol(h1_ME1))] = apply((h2[, (5 * j)]/data$rmpw.ME * (data$ME/data$pME1 - (1 - data$ME)/(1 - data$pME1)) * data$pME0 * (1 - data$pME0) * V_ME0)[data$R==1, ], 2, sum)/nrow(data)
      G21[5 * j, (ncol(G11) - ncol(h1_ME1) + 1):ncol(G11)] = apply((h2[, (5 * j)]/data$rmpw.ME * (-data$ME * data$pME0/data$pME1 * (1 - data$pME1) + (1 - data$ME) * (1 - data$pME0)/(1 - data$pME1) * data$pME1) * V_ME1)[data$R==1, ], 2, sum)/nrow(data)
    }
    G12 = matrix(0, nrow(G11), ncol(G22))
    G= cbind(rbind(G11, G21), rbind(G12, G22))
    V = 1/N * ginv(G) %*% H %*% t(ginv(G))
    V_mu = V[(ncol(V) - 5 * J + 1):ncol(V), (ncol(V) - 5 * J + 1):ncol(V)]  
    
    ######## Population Average Effect Estimation
    Phi = matrix(c(0, 0, -1, 0, 0, 0, -1,
                   0, 1, 0, 1, 0, 1, 1,
                   0, 0, 0, -1, 1, -1, 0, 
                   1, -1, 0, 0, 0, -1, 0,
                   -1, 0, 1, 0, -1, 1, 0), 7, 5)
    beta = kronecker(diag(1, J), Phi) %*% mu
    V_beta = kronecker(diag(1, J), Phi) %*% V_mu %*% kronecker(diag(1, J), t(Phi))
    
    Psi = kronecker(rep(1, J), diag(1, 7))
    gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% beta
    
    ######## Between-Site Variance Estimation
    B = 0
    W = 1/(J * (J - 1)) * t(Psi) %*% V_beta %*% Psi
    for(j in 1:J){
      B = B + 1/(J - 1) * (beta[(7 * j - 6):(7 * j)] - gamma) %*% t((beta[(7 * j - 6):(7 * j)] - gamma))
      W = W - 1/(J - 1) * V_beta[(7 * j - 6):(7 * j), (7 * j - 6):(7 * j)]
    }
    tau = B + W
    
    ######## Standard Error of Population Average Effects
    V_gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% (V_beta + kronecker(diag(1, J), tau)) %*% Psi %*% solve(t(Psi) %*% Psi)
    SE_gamma = sqrt(diag(V_gamma))
    t_gamma = gamma/SE_gamma
    p_gamma = (1 - pnorm(abs(t_gamma))) * 2
    est_gamma = cbind(gamma, SE_gamma, t_gamma, p_gamma)
    est_gamma = as.data.frame(round(est_gamma, 3))
    colnames(est_gamma) = c("Estimate", "Std.Error", "t value", "Pr(>|t|)")
    rownames(est_gamma) = c("I.M1(0)", "I.M2(1)", "D", "I.M1(1)", "I.M2(0)", "I.M1*M2", "ITT")
    
    ######## Between-Site Variance Output
    tau.ori = tau
    tau.sub = tau
    
    if("TRUE" %in% (diag(tau.ori) < 0)){
      tau[which(diag(tau < 0)), ] = 0
      tau[, which(diag(tau < 0))] = 0
      tau.sub = tau.ori[-which(diag(tau.ori < 0)), ]
      tau.sub = tau.sub[, -which(diag(tau.sub < 0))]    
    }
    
    cor.sub = cov2cor(tau.sub)
    if(max(na.omit(cor.sub)) > 1)
      cor.sub = suppressWarnings(cor.smooth(cor.sub))
    cor = cor.sub
    
    if("TRUE" %in% (diag(tau.ori) < 0)){
      cor = tau.ori
      cor[which(diag(tau.ori < 0)), ] = NaN
      cor[, which(diag(tau.ori < 0))] = NaN
      cor[which(diag(tau.ori >= 0)), which(diag(tau.ori >= 0))] = cor.sub
    }
    
    var = round(cbind(diag(tau), sqrt(diag(tau)), cor[, 1:6]), 3)
    var[1, 3:8] = ""
    var[2, 4:8] = ""
    var[3, 5:8] = ""
    var[4, 6:8] = ""
    var[5, 7:8] = ""
    var[6, 8] = ""
    var = as.data.frame(var)
    colnames(var) = c("Variance", "Std.Dev.", "Corr", "", "", "", "", "")
    rownames(var) = c("I.M1(0)", "I.M2(1)", "D", "I.M1(1)", "I.M2(0)", "I.M1*M2", "ITT")
    
    ######## Chisq for Between-Site Variance Testing
    chisq = 0
    for (j in 1:J) {
      Vj = V_beta[(7 * j - 6):(7 * j), (7 * j - 6):(7 * j)]
      chisq = chisq + (beta[(7 * j - 6):(7 * j)] - gamma)^2/diag(Vj)
    }
    
    if(permutation == F)
      return(list(Random_effects = var, Fixed_effects = est_gamma, chisq = chisq))
    if(permutation == T)
      return(list(chisq = chisq))
  }
  
  result = suppressWarnings({est(data, y, treatment, mediator, response, XR1, XR0, XM11, XM10, XM21, XM20, site, sample.weight, permutation = F)})
  
  siteID = data$site
  data_pm = data
  chisq_pm = NULL
  for (i in 1:npermute) {
    try({
      data_pm$site = sample(siteID, length(siteID), replace = F)
      chisq_pm = cbind(chisq_pm, suppressWarnings({est(data_pm, y, treatment, mediator, response, XR1, XR0, XM11, XM10, XM21, XM20, site, sample.weight, permutation = T)})$chisq)
    }, silent = T)
  }
  pvalue = round(apply((chisq_pm - matrix(rep(result$chisq, ncol(chisq_pm)), 7, ncol(chisq_pm))) >= 0, 1, mean), 3)
  
  sig = NULL
  sig[pvalue <= 0.001] = "**"
  sig[pvalue > 0.001 & pvalue <= 0.01] = "*"
  sig[pvalue > 0.01 & pvalue <= 0.05] = "."
  sig[pvalue > 0.05] = ""
  result$Random_effects = cbind(result$Random_effects, pvalue, sig)
  colnames(result$Random_effects)[4:8] = ""
  colnames(result$Random_effects)[7] = "p-value"
  
  return(list(Random_effects = result$Random_effects, Fixed_effects = result$Fixed_effects))
}
  
#' Balance checking for causal mediation analysis in multisite trials
#' 
#' This function is used to check if, within a treatment group, the estimated nonresponse weight balances the distribution of the observed covariates between the respondents and the nonrespondents, or if the estimated RMPW weight balances the distribution of the observed covariates between those whose mediator takes value 1 and those whose mediator takes value 0. 
#'
#' @param data The data set for analysis.
#' @param y The name of the outcome variable (string).
#' @param treatment The name of the treatment variable (string).
#' @param mediator The name of the mediator variable (string).
#' @param response The name of the response variable (string), which is equal to 1 if the individual responded and 0 otherwise.
#' @param XR1 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XR0 A vector of variable names (string) of pretreatment covariates in the propensity score model for the response under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM1 A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the treatment condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param XM0 A vector of variable names (string) of pretreatment covariates in the propensity score model for the mediator under the control condition. For now, the multilevel propensity score model only allows for one random intercept.
#' @param X A vector of variable names (string) of all the pretreatment covariates to be checked balance for.
#' @param site The variable name for the site ID (string).
#' @return A list of tables containing the balance checking results for the response before weighting ($balance.R$balance1 under the treatment condition and $balance.R$balance0 under the control condition) and after weighting ($balance.R$balance1.adj under the treatment condition and $balance.R$balance0.adj under the control condition); and the balance checking results for the mediator before weighting ($balance.M$balance1 under the treatment condition and $balance.M$balance0 under the control condition) and after weighting ($balance.M$balance1.adj under the treatment condition and $balance.M$balance0.adj under the control condition). It also contains a set of balance checking plots corresponding to the tables.
#' \item{balance.mean}{Population average of standardized bias. The standardized bias is calculated by dividing the unweighted (before weighting) or weighted (after weighting) mean difference between response or mediator levels in each covariate by the standard deviation of the covariate }
#' \item{balance.sd}{Between-site standard deviation of standardized bias.}
#' \item{balance.lower}{Lower bound of the 95\% plausible value range of the site-specific standardized bias.}
#' \item{balance.upper}{Upper bound of the 95\% plausible value range of the site-specific standardized bias.}
#' @author Xu Qin, Guanglei Hong, Jonah Deutsch, and Edward Bein
#' @references Qin, X., Hong, G., Deutsch, J., & Bein, E. (2019). Multisite causal mediation analysis in the presence of complex sample and survey designs and non-random non-response. Journal of the Royal Statistical Society: Series A (Statistics in Society), 182(4), 1343-1370. \doi{10.1111/rssa.12446}
#' @export
#' @importFrom stats as.formula binomial coef fitted lm model.matrix pnorm sd
#' @importFrom lme4 VarCorr fixef glmer lmer ranef glmerControl
#' @importFrom ggplot2 ggplot aes geom_errorbarh geom_point labs theme element_blank geom_vline xlim
#' @examples 
#' data(sim.weights)
#'
#' balance(data = sim.weights, y = "y", treatment = "tr", mediator = "me", response = "R", XR1 = c("x1", "x2", "x3"), XR0 = c("x1", "x2", "x3"), XM1 = c("x1", "x2", "x3"), XM0 = c("x1", "x2", "x3"), X = c("x1", "x2", "x3"), site = "site")
                     
balance = function(data, y, treatment, mediator, response, XR1, XR0, XM1, XM0, X, site){
  data = as.data.frame(data)
  colnames(data)[which(colnames(data) == y)] = "y"
  colnames(data)[which(colnames(data) == treatment)] = "tr"
  colnames(data)[which(colnames(data) == mediator)] = "M"
  colnames(data)[which(colnames(data) == response)] = "R"
  colnames(data)[which(colnames(data) == site)] = "site"
  data = data[order(data$site), ]
  # # Factorize categorical covariates (with fewer than 10 categories)
  # transform = function(X){
  #   for(i in 1:length(X)){
  #     if(length(unique(data[, X[i]])) > 2 & length(unique(data[, X[i]])) < 10){
  #       data[, X[i]] = as.factor(data[, X[i]])
  #     }
  #   }
  #   covariates = model.matrix(as.formula(paste("~", paste(X, collapse = "+"))), data)
  #   X = colnames(covariates)
  #   return(list(covariates = covariates[, -1], X = X[-1]))
  # }
  # transform.XR1 = transform(XR1)
  # transform.XR0 = transform(XR0)
  # transform.XM1 = transform(XM1)
  # transform.XM0 = transform(XM0)
  # data = data[, -which(colnames(data) %in% unique(c(XR1, XR0, XM1, XM0)))]
  # XR1 = transform.XR1$X
  # XR0 = transform.XR0$X
  # XM1 = transform.XM1$X
  # XM0 = transform.XM0$X
  # covariates = cbind(transform.XR1$covariates, transform.XR0$covariates, transform.XM1$covariates, transform.XM0$covariates)
  # colnames(covariates) = c(XR1, XR0, XM1, XM0)
  # data = cbind(data, covariates)
  # data = data[, colnames(unique(as.matrix(data), MARGIN = 2))] 
  data1 = data[data$tr == 1, ]
  data0 = data[data$tr == 0, ]
  
  #### Balance checking
  check = function(var, X1, X0, X){
    if(var == "M"){
      data = data[data$R == 1, ]
      data1 = data[data$tr == 1, ]
      data0 = data[data$tr == 0, ]
    }

    ## Model fitted to the treatment group:
    l1 = glmer(as.formula(paste(var, "~", paste(X1, collapse="+"), "+(1|site)")), data = data1, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    ## Model fitted to the control group:
    l0 = glmer(as.formula(paste(var, "~", paste(X0, collapse="+"), "+(1|site)")), data = data0, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    #### Predict the response probabilities
    data$p[data$tr == 1 & data[, var] == 1] = fitted(l1)[data1[, var] == 1]
    data$p[data$tr == 1 & data[, var] == 0] = 1 - fitted(l1)[data1[, var] == 0]
    data$p[data$tr == 0 & data[, var] == 1] = fitted(l0)[data0[, var] == 1]
    data$p[data$tr == 0 & data[, var] == 0] = 1 - fitted(l0)[data0[, var] == 0]
    ## Construct the IPTW weight
    l1 = glmer(as.formula(paste(var, "~ 1 + (1|site)")), data = data1, family = binomial)
    l0 = glmer(as.formula(paste(var, "~ 1 + (1|site)")), data = data0, family = binomial)
    data$iptw[data$tr == 1 & data[, var] == 1] = fitted(l1)[data1[, var] == 1]/data$p[data$tr == 1 & data[, var] == 1]
    data$iptw[data$tr == 1 & data[, var] == 0] = (1 - fitted(l1)[data1[, var] == 0])/data$p[data$tr == 1 & data[, var] == 0]
    data$iptw[data$tr == 0 & data[, var] == 1] = fitted(l0)[data0[, var] == 1]/data$p[data$tr == 0 & data[, var] == 1]
    data$iptw[data$tr == 0 & data[, var] == 0] = (1 - fitted(l0)[data0[, var] == 0])/data$p[data$tr == 0 & data[, var] == 0]
    ## Add logit of propensity scores to the vector of covariates
    data$logitp[data$tr == 1 & data[, var] == 1] = log(data$p[data$tr == 1 & data[, var] == 1]/(1 - data$p[data$tr == 1 & data[, var] == 1]))
    data$logitp[data$tr == 1 & data[, var] == 0] = log((1 - data$p[data$tr == 1 & data[, var] == 0])/data$p[data$tr == 1 & data[, var] == 0]) 
    data$logitp[data$tr == 0 & data[, var] == 1] = log(data$p[data$tr == 0 & data[, var] == 1]/(1 - data$p[data$tr == 0 & data[, var] == 1]))
    data$logitp[data$tr == 0 & data[, var] == 0] = log((1 - data$p[data$tr == 0 & data[, var] == 0])/data$p[data$tr == 0 & data[, var] == 0]) 
    data1 = data[data$tr == 1, ]
    data0 = data[data$tr == 0, ]
    X = c(X, "logitp")
    
    #### Calculate population average and between-site standard deviation of the standardized bias
    bias = function(X, weights = NULL, data){
      balance.mean = NULL
      balance.sd = NULL
      balance.lower = NULL
      balance.upper = NULL
      for(k in 1:length(X)){
        if(length(unique(data[, X[k]])) == 2){
          l.balance = suppressWarnings({glmer(as.formula(paste(X[k], "~", var, "+ (", var, "|site)")), data = data, weights = weights, family = binomial, nAGQ = 1)})
          r0j = fixef(l.balance)[1] + ranef(l.balance)$site[, 1]
          r1j = fixef(l.balance)[2] + ranef(l.balance)$site[, 2]
          mean0j = (1 + exp(-r0j))^(-1)
          mean1j = (1 + exp(-r0j - r1j))^(-1)
          mean.difference = mean(mean1j - mean0j)
          sd.difference = sd(mean1j - mean0j)
        } else {
          l.balance = suppressWarnings({lmer(as.formula(paste(X[k], "~", var, "+ (", var, "|site)")), data = data, weights = weights)})
          mean.difference = fixef(l.balance)[2]
          sd.difference = as.data.frame(VarCorr(l.balance))$sdcor[2]
        }
        balance.mean = c(balance.mean, mean.difference/sd(data[ , X[k]]))
        balance.sd = c(balance.sd, sd.difference/sd(data[ , X[k]]))
        balance.lower = c(balance.lower, (mean.difference - 1.96 * sd.difference)/sd(data[ , X[k]]))
        balance.upper = c(balance.upper, (mean.difference + 1.96 * sd.difference)/sd(data[ , X[k]]))
      }
      results = cbind(balance.mean = balance.mean, balance.sd = balance.sd, balance.lower = balance.lower, balance.upper = balance.upper)
      rownames(results) = X
      return(results)
    }
    
    #### Report the population average and between-site standard deviation of the standardized bias before and after weighting for each pretreatment covariate under each treatment condition
    return(
      list(balance1 = round(bias(X, weights = NULL, data1), 3),
           balance1.adj = round(bias(X, weights = data1$iptw, data1), 3),
           balance0 = round(bias(X, weights = NULL, data0), 3),
           balance0.adj = round(bias(X, weights = data0$iptw, data0), 3))
    )
  }
  
  balance.R = check("R", XR1, XR0, X)
  balance.M = check("M", XM1, XM0, X)
  
  balance.plot = function(balance, var, treatment, adj){
    meandiff = balance[, 1]
    min = balance[, 3]
    max = balance[, 4]
    covariates = rownames(balance)
    covariates = factor(covariates, levels = rownames(balance)[order(abs(meandiff))])
    balance = cbind.data.frame(covariates = covariates, meandiff = meandiff, min = min, max = max)
    
    ggplot(balance, aes(meandiff, covariates)) +
      geom_errorbarh(aes(y = covariates, x = meandiff, xmin = min, xmax = max), data = balance, col="grey", size=1.2) +
      geom_point(size=3, shape=21, col="grey", fill = "grey") +
      labs(x = "Standardized Mean Differences", y = "", title = paste("Imbalance Between", var, "Levels", adj, "Weighting in", treatment, "Group")) +
      theme(legend.title=element_blank()) +
      geom_vline(xintercept=c(-0.25,-0.1,0.1,0.25), linetype="longdash", col = c("red", "blue", "blue", "red")) +
      xlim(-2.5, 2.5)
  }
  
  plot.balance.R.1 = balance.plot(balance.R$balance1, "Response", "Treatment", "Before")
  plot.balance.R.1.adj = balance.plot(balance.R$balance1.adj, "Response", "Treatment", "After")
  plot.balance.R.0 = balance.plot(balance.R$balance0, "Response", "Control", "Before")
  plot.balance.R.0.adj = balance.plot(balance.R$balance0.adj, "Response", "Control", "After")
  plot.balance.M.1 = balance.plot(balance.M$balance1, "Mediator", "Treatment", "Before")
  plot.balance.M.1.adj = balance.plot(balance.M$balance1.adj, "Mediator", "Treatment", "After")
  plot.balance.M.0 = balance.plot(balance.M$balance0, "Mediator", "Control", "Before")
  plot.balance.M.0.adj = balance.plot(balance.M$balance0.adj, "Mediator", "Control", "After")
  
  return(
    list(balance.R = balance.R,
         balance.M = balance.M,
         plot.balance.R.1 = plot.balance.R.1,
         plot.balance.R.1.adj = plot.balance.R.1.adj,
         plot.balance.R.0 = plot.balance.R.0,
         plot.balance.R.0.adj = plot.balance.R.0.adj,
         plot.balance.M.1 = plot.balance.M.1,
         plot.balance.M.1.adj = plot.balance.M.1.adj,
         plot.balance.M.0 = plot.balance.M.0,
         plot.balance.M.0.adj = plot.balance.M.0.adj)
  )
}